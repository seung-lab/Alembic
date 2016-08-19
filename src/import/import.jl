global PREVIOUS_OVERVIEW_RESOLUTION = 86/3000
global PREVIOUS_IMPORT_BB = ImageRegistration.BoundingBox(10110,19850,39000,28000)
global OVERVIEW_RESOLUTION = 95.3/3840 # 3.58/225.0 
global OVERVIEW_IMPORT_BB = ImageRegistration.BoundingBox(290,570,1110,800)
global IMPORT_BB = snap_bb(scale_bb(OVERVIEW_IMPORT_BB, 1/OVERVIEW_RESOLUTION))

function get_src_dir(z_index)
	LOADFILE = readdlm("/media/tmacrina/667FB0797A5072D7/3D_align/mosaiced_images_160729_google_cloud_upload_seung_import.csv",',')
	i = findfirst(i -> i == z_index, LOADFILE[:,1])
	if i > 0
		return joinpath("/media/tmacrina/", LOADFILE[i, 2])
	else
		return nothing
	end
end

function load_metadata(z_index)
	dir = get_src_dir(z_index)
	files = readdir(dir)
	meta_fn = files[findfirst(i -> contains(i, "meta"), files)]
	return JSON.parsefile(joinpath(dir, meta_fn), dicttype=Dict, use_mmap=true)
end

function plot_meta_mean(z_index)
	meta = load_metadata(z_index)
	means = [i["mean"] for i in meta[2]["data"]]
	return plt[:scatter](1:length(means), means)
end

function load_contrast_histogram(z_index; reset=false)
	auto_fn = get_path("contrast_stretch", premontaged(1, z_index))
	if !isfile(auto_fn) || reset
		calculate_contrast_histogram(z_index)
	end
	return readdlm(auto_fn)
end

function load_import_table(z_index; reset=false)
	import_fn = get_path("import", premontaged(1, z_index))
	if !isfile(import_fn) || reset
		dir = get_src_dir(z_index)
		files = readdir(dir)
		trakem_fn = files[findfirst(i -> contains(i, "trackem"), files)]
		initialize_import_table(joinpath(dir, trakem_fn), z_index)
	end
	import_table = readdlm(import_fn)
	update_import_include!(z_index, import_table)
	return import_table
end

function initialize_import_table(import_src_fn, z_index)
	println("Initializing import table for section $z_index from $import_src_fn")
	import_table = readdlm(import_src_fn)

	# pull out row & col from filename
	path_no_ext = [split(i, ".")[1] for i in import_table[:,1]]
	waf = [1 for i in path_no_ext]
	col = [parse(split(i, "_")[end-1])+1 for i in path_no_ext]
	row = [parse(split(i, "_")[end])+1 for i in path_no_ext]
	sec = [z_index for i in path_no_ext]
	height = [0 for i in path_no_ext]
	width = [0 for i in path_no_ext]
	include = [true for i in path_no_ext]
	roi = [true for i in path_no_ext]

	# include full path to the file
	dir = get_src_dir(z_index)
	import_table[:,1] = [joinpath(dir, fn) for fn in import_table[:,1]]
	# dst_path = map(get_path, zip(waf, sec, row, col))
	import_table = hcat(import_table, [waf sec row col height width include roi])
	save_import_table(z_index, import_table)
end

function save_import_table(z_index, import_table)
	import_dst_fn = get_path("import", premontaged(1, z_index))
	println("Writing import table for section $z_index to $import_dst_fn")
	writedlm(import_dst_fn, import_table)
end

function bias_correction(img, bias, bias_mean)
	# return (img ./ bias) * bias_mean
	return img - bias + bias_mean
end

function contrast_stretch(img, minintensity, maxintensity)
	return min(1.0, max(0.0, (img-minintensity) / (maxintensity-minintensity)))
end

function auto_contrast_stretch(img, bins=20)
	hist = nquantile(img[:], bins)
	minintensity = hist[2]
	maxintensity = hist[end-1]
	return min(1.0, max(0.0, (img-minintensity) / (maxintensity-minintensity)))
end

function convert_float64_to_uint8(img)
	return Array{UInt8, 2}(round(UInt8, img * 255))
end

# function generate_clahe(clip_limit=3.0, img_size=(3840,3840), block_size=128)
# 	@pyimport cv2
# 	tile_grid_size = round(Int64, img_size[1] / block_size), round(Int64, img_size[2] / block_size)
# 	return cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=tile_grid_size)
# end

function apply_clahe(img, clahe)
	return clahe[:apply](img)
end

function adjust_tile(z_index, i)
	import_table = load_import_table(z_index)
	bias = load_bias_image(z_index)
	bias_mean = mean(bias)
	dtype = Float64
	# contrast_hist = load_contrast_histogram(z_index, reset=reset)
	# minval = contrast_hist[2]
	# maxval = contrast_hist[length(contrast_hist)-1]
	maxval = dtype(typemax(UInt16))

	src_fn = get_src_fn(import_table,i)
	if isfile(src_fn)
		img = get_image_disk(src_fn, dtype)
		img = bias_correction(img, bias, bias_mean)
		img = contrast_stretch(img)
		img = convert_float64_to_uint8(img)
		return img
	end
	return nothing
end

function load_bias_image(z_index; reset=false)
	bias_fn = get_path("contrast_bias", premontaged(1,z_index))
	if !isfile(bias_fn) || reset
		calculate_section_bias_field(z_index)
	end
	return get_image_disk(bias_fn, Float64)
end

function calculate_section_bias_field(z_index, N=100)
	println("Calculating contrast bias field for $z_index using $N samples")
	indices = []
	tiles = Array{Array{Float64,2}, 1}()

	# import_table = map(load_import_table, z_indices)
	# import_table = vcat(import_table...)
	import_table = load_import_table(z_index)
	import_indices = get_included_indices(import_table)
	samples = rand(import_indices, N)

	bias = @parallel (+) for i=samples
		println(i)
		Images.imfilter_gaussian_no_nans!(get_image_disk(get_src_fn(import_table, i), Float64), [10,10]) / N
	end

	bias[bias .== 0] = 1
	bias_fn = get_path("contrast_bias", premontaged(1,z_index))
	f = h5open(bias_fn, "w")
	f["img"] = bias
	close(f)
	println("Contrast bias image written to $bias_fn")
end

function calculate_contrast_histogram(z_index, N=100, bins=20)
	println("Calculating overall contrast histogram for section $(z_index) using $N samples")
	import_table = load_import_table(z_index)	
	bias = load_bias_image(z_index)
	bias_mean = mean(bias)
	dtype = Float64

	import_indices = get_included_indices(import_table)
	samples = rand(import_indices, N)

	function sample_image(i)
	# for i in samples
		println("Loading sample $i / $(length(samples))")
		img = get_image_disk(get_src_fn(import_table, i), dtype)
		smp = bias_correction(img, bias, bias_mean)[3840*1900:3840*1920]
		# push!(intensities, smp)
		return smp
	end

	intensities = pmap(sample_image, samples)
	intensities = vcat(intensities...)
	distr = nquantile(intensities, bins)
	stretch_fn = get_path("contrast_stretch", premontaged(1, z_index))
	writedlm(stretch_fn, distr)
	println("Contrast histogram written to $stretch_fn")
	return intensities
end

function import_overview(z_index)
	dir = get_src_dir(z_index)
	if dir != nothing
		files = readdir(dir)
		k = findfirst(i -> contains(i, "montage"), files)
		index = (1,z_index,OVERVIEW_INDEX,OVERVIEW_INDEX)
		if k > 0
			src_fn = joinpath(dir, files[k])
			dst_fn = string(get_path(index)[1:end-3], ".tif")
			run(`cp $src_fn $dst_fn`)
		end
	end
end

function import_tiles(z_index; reset=false)
	# try
	import_table = load_import_table(z_index)
	N = size(import_table, 1)

	import_indices = get_included_indices(import_table)
	if !reset
		import_indices = get_unwritten_included_indices(import_table)
	end
	n = length(import_indices)
	# contrast_clusters = get_contrast_clusters(import_table)
	bias = load_bias_image(z_index, reset=reset)
	bias_mean = mean(bias)
	# contrast_hist = load_contrast_histogram(z_index, reset=reset)
	# minintensity = contrast_hist[2]
	# maxintensity = contrast_hist[end-1]
	dtype = Float64

	function import_tile(i)
		src_fn = get_src_fn(import_table, i)
		dst_fn = get_dst_fn(import_table, i)
		index = get_import_index(import_table, i)
		println("Importing $index, # $i / $n")
		println("\tsrc: $src_fn")

		img = get_image_disk(src_fn, dtype)
		sz = size(img)
		img = bias_correction(img, bias, bias_mean)
		# img = contrast_stretch(img, minintensity, maxintensity)
		img = auto_contrast_stretch(img)
		img = convert_float64_to_uint8(img)
		# resin = is_adjusted_resin(img, 20, z_index)
		println("\tdst: $dst_fn")
		f = h5open(dst_fn, "w")
		chunksize = 1000
		@time f["img", "chunk", (chunksize,chunksize)] = img
		close(f)
		return i, sz...
	end

	results = pmap(import_tile, import_indices)

	indices = [i[1] for i in results]
	order = sortperm(indices)
	height = [i[2] for i in results][order]
	width = [i[3] for i in results][order]
	# resin = [i[4] for i in results]

	println("Appending image size to import table")
	import_table[import_indices,9] = height
	import_table[import_indices,10] = width

	import_dst_fn = get_path("import", premontaged(1, z_index))
	println("Writing import table for section $z_index to $import_dst_fn")
	writedlm(import_dst_fn, import_table)
	initialize_offsets(z_index)
	# end
end

function downsample_tiles(z_index, thumbnail_scale=0.25)
	import_table = load_import_table(z_index)
	N = size(import_table, 1)
	existing_indices = create_src_isfile_mask(import_table)
	indices = collect(1:N)[existing_indices]
	n = length(indices)
	dtype = Float64

	function downsample_tile(i)
		src_fn = get_src_fn(import_table, i)
		index = get_import_index(import_table, i)
		println("Downsampling $index, # $i / $n")
		img = get_image_disk(src_fn, dtype)
		thumbnail, _ = imscale(img, thumbnail_scale)
	    write_thumbnail(thumbnail, index, thumbnail_scale)
	    return true
	end

	results = pmap(downsample_tile, indices)
end

function concat_thumbs(z_index, M=10)
	import_table = load_import_table(z_index)
	N = size(import_table, 1)
	existing_indices = create_src_isfile_mask(import_table)
	sample = rand(collect(1:N)[existing_indices], M)
	indices = map(get_import_index, repeated(import_table), sample)
	n = length(indices)

	thumb_array = zeros(Float64, 960, 960, n)

	@sync @parallel for (i, index) in enumerate(indices)
		@async thumb_array[:,:,i] = get_image_disk(get_path("thumbnail", index), Float64)
	end

	return thumb_array
end

function get_src_fn(import_table, i)
	return import_table[i,1]
end

function get_dst_fn(import_table, i)
	index = get_import_index(import_table, i)
	return get_path(index)
end

function get_src_filenames(import_table)
	return [get_src_fn(import_table, i) for i in 1:size(import_table,1)]
end

function get_dst_filenames(import_table)
	return [get_dst_fn(import_table, i) for i in 1:size(import_table,1)]
end

function get_row_numbers(import_table)
	return import_table[:,7]
end

function get_col_numbers(import_table)
	return import_table[:,8]
end

function get_import_index(import_table, i)
	return (import_table[i,5:8]...)
end

function get_import_indices(import_table::Array{Any,2})
	return [get_import_index(import_table, i) for i in 1:size(import_table,1)]
end

# function get_import_indices(z_index::Int64)
# 	import_table = load_import_table(z_index)
# 	return get_import_indices(import_table)
# end

function get_import_offset(import_table, i)
	return [import_table[i,3], import_table[i,2]]
end

function update_import_offset!(import_table, i, new_offset)
	import_table[i,3], import_table[i,2] = new_offset
end

function get_import_offsets(import_table)
	return [get_import_offset(import_table, i) for i in 1:size(import_table,1)]
end

function create_dst_isfile_mask(import_table)
	return convert(BitArray, map(isfile, get_dst_filenames(import_table)))
end

function create_src_isfile_mask(import_table)
	return convert(BitArray, map(isfile, get_src_filenames(import_table)))
end

function get_import_sizes(import_table)
	return [(import_table[i,9:10]...) for i in 1:size(import_table,1)]
end

# function get_contrast_clusters(import_table)
# 	return unique(import_table[:,12])
# end

# function get_contrast_id(import_table, i)
# 	return import_table[i,12]
# end

# function create_contrast_cluster_mask(import_table, i)
# 	return convert(BitArray, import_table[:,12] .== i)
# end

function initialize_offsets(z_index)
	import_table = load_import_table(z_index)
	imported_tiles = create_dst_isfile_mask(import_table)
	indices = get_import_indices(import_table)
	offsets = get_import_offsets(import_table)
	# sizes = get_import_sizes(import_table)
	sizes = [(3840,3840) for i in 1:size(import_table,1)]
	update_offsets(indices[imported_tiles], offsets[imported_tiles], sizes[imported_tiles])
end

function setup_registry_after_python_import(dir)
	indices = []
	offsets = []
	sizes = []
	for fn in readdir(dir)
		if fn[end-2:end] == ".h5"
			sz = h5read(joinpath(dir, fn), "size")
			index = parse_name(fn)
			println("$sz, $index")
			push!(indices, index)
			push!(offsets, [0,0])
			push!(sizes, sz)
		end
	end
	update_offsets(indices, offsets, sizes)
end

function get_bbs(import_table, sz=[3840,3840])
	offsets = get_import_offsets(import_table)
	return [ImageRegistration.BoundingBox(offset..., sz...) for offset in offsets]
end

function get_tform_bbs(import_table, tform=eye(3), sz=[3840,3840])
	bbs = get_bbs(import_table, sz)
	return [tform_bb(bb, tform) for bb in bbs]
end

function get_bbs_pts(bbs)
	return [[bb_to_pts(bb)] for bb in bbs]
end

function tform_bbs_pts(bbs_pts, tform)
	tform_pts = [[pts ones(size(pts,1),1)]*tform for pts in bbs_pts]
	return [pts[:, 1:2] for pts in tform_pts]
end

function snap_bbs_pts(bbs_pts)
	return [round(Int64, pts) for pts in bbs_pts]
end

function get_tform(index::Index, scale=1.0)
	s = make_scale_matrix(scale)
	return s*load("cumulative_transform", index)*s^-1
end

function view_import_bb(z_index::Int64, tform=eye(3))
	import_table = load_import_table(z_index)
	polys = make_import_outline(import_table, tform)
	# polys = [poly[1:4,:] for poly in polys]
	polys = tform_bbs_pts(polys, make_scale_matrix(0.05))
	indices = get_import_indices(import_table)
	view_polys(polys, indices)
end

function view_roi_previous(z_index::Int64)
	roi = PREVIOUS_IMPORT_BB
	tform = get_tform(overview(1,z_index), PREVIOUS_OVERVIEW_RESOLUTION)
	return view_roi(z_index, roi=roi, tform=tform)
end

function view_roi(z_index::Int64; roi::ImageRegistration.BoundingBox=IMPORT_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	scale = 0.05
	import_table = load_import_table(z_index)
	polys = make_import_outline(import_table, tform)
	# polys = [poly[1:4,:] for poly in polys]
	polys = tform_bbs_pts(polys, make_scale_matrix(scale))
	indices = get_import_indices(import_table)
	view_polys(polys, indices, bb_to_pts(scale_bb(roi, scale)))
end

function make_import_outline(import_table, tform=eye(3))
	bbs = get_bbs(import_table)
	pts = get_bbs_pts(bbs)
	return tform_bbs_pts(pts, tform)
end

function get_roi_mask(import_table, roi::ImageRegistration.BoundingBox=IMPORT_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	# possible speed up: run bounding box intersection first
	# all_bbs = get_tform_bbs(z_index)
	# roi_intersects = convert(BitArray, [intersects(bb, roi) for bb in all_bbs])
	# bbs = all_bbs[roi_intersects]
	# indices = all_indices[roi_intersects]
	tiles = make_import_outline(import_table, tform)
	return convert(BitArray, [poly_intersects(tile, bb_to_pts(roi)) for tile in tiles])
end

function get_tile_indices(z_index::Int64, roi::ImageRegistration.BoundingBox=IMPORT_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	import_table = load_import_table(z_index)
	all_indices = get_import_indices(import_table)
	roi_mask = get_roi_mask(import_table, roi, tform)
	return all_indices[roi_mask]
end

function update_import_include!(z_index, import_table, roi::ImageRegistration.BoundingBox=IMPORT_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	include_mask = create_flagged_mask(import_table)
	roi_mask = get_roi_mask(import_table, roi, tform)
	exisiting_mask = create_src_isfile_mask(import_table)
	import_table[:,12] = include_mask & roi_mask & exisiting_mask
	save_import_table(z_index, import_table)
end

function get_included_indices(import_table)
	indices =  collect(1:size(import_table,1))
	mask = create_included_mask(import_table)
	return indices[mask]
end

function get_included_indices(import_table)
	indices =  collect(1:size(import_table,1))
	mask = create_included_mask(import_table)
	return indices[mask]
end

function get_unwritten_included_indices(import_table)
	indices =  collect(1:size(import_table,1))
	included = create_included_mask(import_table)
	unwritten = create_dst_isfile_mask(import_table)
	mask = included & unwritten
	return indices[mask]
end

function create_included_mask(import_table; reset=true)
	return convert(BitArray, import_table[:,12])
end

function create_flagged_mask(import_table; reset=true)
	return convert(BitArray, import_table[:,11])
end

function reset_flagged_mask!(import_table)
	import_table[:,11] = trues(length(import_table[:,11]))
end

function rename_files(subdir)
	dir_path = OVERVIEW_DIR_PATH
	files = readdir(dir_path)
	for f in files
		if f[end-1:end] == "h5"
			src_fn = joinpath(dir_path, f)
			index = parse_name(f)
			dst_fn = get_path(overview(index))
			run(`mv $src_fn $dst_fn`)
		end
	end
end

function need_to_reimage_for_bb(z_indices, bb::ImageRegistration.BoundingBox)
	outside_roi = Dict()
	for z_index in z_indices
		try
			bb_tiles = get_tile_indices(z_index, bb)
			# sd = outside_import_roi(z_index, bb)
			outside_roi[z_index] = length(bb_tiles)
		end
	end
	fn = joinpath(homedir(), "Desktop", "tiles_within_section.txt")
	writedlm(fn, outside_roi)
	return outside_roi
end

function outside_import_roi(z_index::Int64)
	import_table = load_import_table(z_index)
	roi_tiles = get_tile_indices(z_index)
	all_tiles = get_import_indices(import_table)
	return setdiff(all_tiles, roi_tiles)
end

"""
Assumes there are no gaps in the tiles for a given row
"""
function shift_row_offsets!(import_table, row_num, shift_in_tiles)
	row = get_row_numbers(import_table)
	indices = collect(1:size(import_table,1))[row .== row_num]
	col = get_col_numbers(import_table)[row .== row_num]
	order = shift_in_tiles > 0 ? sortperm(col) : sortperm(col, rev=true)
	indices = indices[order]
	shift_in_tiles = abs(shift_in_tiles)
	for k in shift_in_tiles+1:length(col)
		offset = get_import_offset(import_table, indices[k])
		update_import_offset!(import_table, indices[k-shift_in_tiles], offset)
	end
	delta_offset = shift_in_tiles*(get_import_offset(import_table, indices[2]) - get_import_offset(import_table, indices[1]))
	for k in length(col)-shift_in_tiles+1:length(col)
		offset = get_import_offset(import_table, indices[k]) + delta_offset
		update_import_offset!(import_table, indices[k], offset)
	end
end

function shift_rows!(z_index, row_range, shift_in_tiles)
	import_table = load_import_table(z_index)
	for row in row_range
		shift_row_offsets!(import_table, row, shift_in_tiles)
	end
	update_import_include!(z_index, import_table, reset=true)
	save_import_table(z_index, import_table)
end

function get_overview_size(index::Index)
	path = get_path(overview(index))
	fid = h5open(path, "r")
	dset = fid["img"]
	return [size(dset)...]
end

function get_overview_resolution(index::Index, tile_size=[3840,3840])
	overview_size = get_overview_size(index)
	z_index = index[2]
	import_table = load_import_table(z_index)
	overview_origin = get_overview_origin(import_table)
	last_tile_offset = get_last_offset(import_table) + tile_size
	resolutions = tile_size ./ (last_tile_offset - overview_origin) .* overview_size
	return resolutions ./ tile_size
end

function get_width_in_tiles(import_table)
	return maximum(import_table[:,8])
end

function get_height_in_tiles(import_table)
	return maximum(import_table[:,7])
end

function get_overview_origin(z_index::Int64)
	return get_overview_origin(load_import_table(z_index))
end

function get_overview_origin(import_table)
	return [minimum(import_table[:,3]), minimum(import_table[:,2])]
end

function get_last_offset(import_table)
	return [maximum(import_table[:,3]), maximum(import_table[:,2])]
end

function get_montage_original_offset(index::Index)
	z_index = index[2]
	import_table = load_import_table(z_index)
	tiles = get_index_range(premontaged(index), premontaged(index))
	rc = hcat([[t[3:4]...] for t in tiles]...)'
	min_rc = minimum(rc[:,1]), minimum(rc[:,2])
	indices = get_import_indices(import_table)
	k = findfirst(i->i == (index[1:2]..., min_rc...), indices)
	global_offset = get_import_offset(import_table, k)
	overview_origin = get_overview_origin(import_table)
	return global_offset - overview_origin
end

function get_montage_tform_from_overview(index::Index)
	src_index = index
	dst_index = get_preceding(index)
	tform = load("relative_transform", overview(src_index))
	moving_offset = get_montage_original_offset(src_index)
	fixed_offset = get_montage_original_offset(dst_index)
	moving_t = make_translation_matrix(moving_offset)
	fixed_t = make_translation_matrix(fixed_offset)
	src_scale = make_scale_matrix(get_overview_resolution(src_index)[1])
	dst_scale = make_scale_matrix(get_overview_resolution(dst_index)[1])
	return moving_t*src_scale*tform*dst_scale^-1*fixed_t^-1
end

function initialize_montage_registry_from_overview_tform(index::Index)
	tform = get_montage_tform_from_overview(index)
	update_registry(index, tform)
end

function initialize_montage_registry_from_overview_tform(firstindex::Index, lastindex::Index)
	for index in get_index_range(firstindex, lastindex)
		initialize_montage_registry_from_overview_tform(index)
	end
end

function reset_montage_registry_tform(firstindex::Index, lastindex::Index)
	for index in get_index_range(firstindex, lastindex)
		update_registry(index, rotation=0, offset=[0,0])
	end
end