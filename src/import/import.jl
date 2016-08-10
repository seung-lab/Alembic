#global LOADFILE = readdlm("/media/tmacrina/667FB0797A5072D7/3D_align/mosaiced_images_160729_google_cloud_upload_seung_import.csv",',')
global OVERVIEW_RESOLUTION = 86/3000
global IMPORT_ROI_BB = BoundingBox(10110,19850,39000,28000)

function get_src_dir(z_index)
	i = findfirst(i -> i == z_index, LOADFILE[:,1])
	if i > 0
		return joinpath("/media/tmacrina/", LOADFILE[i, 2])
	else
		return nothing
	end
end

function initialize_import_table(import_src_fn, z_index)
	println("Initializing import table for section $z_index from $import_src_fn")
	import_table = readdlm(import_src_fn)
	import_dst_fn = get_path("import", premontaged(1, z_index))

	# pull out row & col from filename
	path_no_ext = [split(i, ".")[1] for i in import_table[:,1]]
	waf = [1 for i in path_no_ext]
	col = [parse(split(i, "_")[end-1])+1 for i in path_no_ext]
	row = [parse(split(i, "_")[end])+1 for i in path_no_ext]
	sec = [z_index for i in path_no_ext]
	# include full path to the file
	dir = get_src_dir(z_index)
	import_table[:,1] = [joinpath(dir, fn) for fn in import_table[:,1]]
	file_exists = map(isfile, import_table[:,1])
	dst_path = map(get_path, zip(waf, sec, row, col))
	import_table = hcat(import_table, [waf sec row col dst_path file_exists])
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

# function contrast_stretch(img, bins=20)
# 	hist = nquantile(img[:], bins)
# 	minintensity = hist[2]
# 	maxintensity = hist[end-1]
# 	return min(1.0, max(0.0, (img-minintensity) / (maxintensity-minintensity)))
# end

function convert_float64_to_uint8(img)
	return Array{UInt8, 2}(round(UInt8, img * 255))
end

function generate_clahe(clip_limit=3.0, img_size=(3840,3840), block_size=128)
	@pyimport cv2
	tile_grid_size = round(Int64, img_size[1] / block_size), round(Int64, img_size[2] / block_size)
	return cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=tile_grid_size)
end

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

function calculate_section_bias_field(z_index, N=100)
	println("Calculating contrast bias field for $z_index using $N samples")
	indices = []
	tiles = Array{Array{Float64,2}, 1}()

	# import_table = map(load_import_table, z_indices)
	# import_table = vcat(import_table...)
	import_table = load_import_table(z_index)
	import_indices = get_roi_existing_indices(z_index)
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

	import_indices = get_roi_existing_indices(z_index)
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
	return readdlm(import_fn)
end

function load_bias_image(z_index; reset=false)
	bias_fn = get_path("contrast_bias", premontaged(1,z_index))
	if !isfile(bias_fn) || reset
		calculate_section_bias_field(z_index)
	end
	return get_image_disk(bias_fn, Float64)
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
	try
		import_table = load_import_table(z_index, reset=reset)
		N = size(import_table, 1)

		import_indices = get_roi_existing_indices(z_index)
		n = length(import_indices)

		bias = load_bias_image(z_index, reset=reset)
		bias_mean = mean(bias)
		contrast_hist = load_contrast_histogram(z_index, reset=reset)
		minintensity = contrast_hist[2]
		maxintensity = contrast_hist[end-1]
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
			img = contrast_stretch(img, minintensity, maxintensity)
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

		println("Appending image size & resin test results to import table")
		if size(import_table,2) < 11
			import_table = hcat(import_table, [zeros(Int64, N) zeros(Int64, N)])
		end
		import_table[import_indices,11] = height
		import_table[import_indices,12] = width

		import_dst_fn = get_path("import", premontaged(1, z_index))
		println("Writing import table for section $z_index to $import_dst_fn")
		writedlm(import_dst_fn, import_table)
		initialize_offsets(z_index)
	end
end

function downsample_tiles(z_index, thumbnail_scale=0.25)
	import_table = load_import_table(z_index)
	N = size(import_table, 1)
	existing_indices = get_import_src_isfile(import_table)
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
	existing_indices = get_import_src_isfile(import_table)
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
	return import_table[i,9]
end

function get_src_filenames(import_table)
	return [get_src_fn(import_table, i) for i in 1:size(import_table,1)]
end

function get_dst_filenames(import_table)
	return [get_dst_fn(import_table, i) for i in 1:size(import_table,1)]
end

function get_import_index(import_table, i)
	return (import_table[i,5:8]...)
end

function get_import_indices(import_table::Array{Any,2})
	return [get_import_index(import_table, i) for i in 1:size(import_table,1)]
end

function get_import_indices(z_index::Int64)
	import_table = load_import_table(z_index)
	return get_import_indices(import_table)
end

function get_import_offsets(import_table)
	return [[import_table[i,3], import_table[i,2]] for i in 1:size(import_table,1)]
end

function get_import_dst_isfile(import_table)
	return convert(BitArray, map(isfile, get_dst_filenames(import_table)))
end

function get_import_src_isfile(import_table)
	return convert(BitArray, map(isfile, get_src_filenames(import_table)))
end

function get_import_sizes(import_table)
	return [(import_table[i,11:12]...) for i in 1:size(import_table,1)]
end

function initialize_offsets(z_index)
	import_table = load_import_table(z_index)
	imported_tiles = get_import_dst_isfile(import_table)
	indices = get_import_indices(import_table)
	offsets = get_import_offsets(import_table)
	sizes = get_import_sizes(import_table)
	# sizes = [(3840,3840) for i in 1:size(import_table,1)]
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

function get_bbs(z_index::Int64, sz=[3840,3840])
	import_table = load_import_table(z_index)
	offsets = get_import_offsets(import_table)
	return [BoundingBox(offset..., sz...) for offset in offsets]
end

function get_tform_bbs(z_index::Int64, tform=eye(3), sz=[3840,3840])
	bbs = get_bbs(z_index, sz)
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
	polys = make_import_outline(z_index, tform)
	# polys = [poly[1:4,:] for poly in polys]
	polys = tform_bbs_pts(polys, make_scale_matrix(0.05))
	indices = get_import_indices(z_index)
	view_polys(polys, indices)
end

function view_roi(z_index::Int64, roi::BoundingBox=IMPORT_ROI_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	scale = 0.05
	polys = make_import_outline(z_index, tform)
	# polys = [poly[1:4,:] for poly in polys]
	polys = tform_bbs_pts(polys, make_scale_matrix(scale))
	indices = get_import_indices(z_index)
	view_polys(polys, indices, bb_to_pts(scale_bb(roi, scale)))
end

function make_import_outline(z_index::Int64, tform=eye(3))
	bbs = get_bbs(z_index)
	pts = get_bbs_pts(bbs)
	return tform_bbs_pts(pts, tform)
end

function get_roi_mask(z_index::Int64, roi::BoundingBox=IMPORT_ROI_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	# possible speed up: run bounding box intersection first
	# all_bbs = get_tform_bbs(z_index)
	# roi_intersects = convert(BitArray, [intersects(bb, roi) for bb in all_bbs])
	# bbs = all_bbs[roi_intersects]
	# indices = all_indices[roi_intersects]
	tiles = make_import_outline(z_index, tform)
	return convert(BitArray, [poly_intersects(tile, bb_to_pts(roi)) for tile in tiles])
end

function get_tile_indices(z_index::Int64, roi::BoundingBox=IMPORT_ROI_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	all_indices = get_import_indices(z_index)
	roi_mask = get_roi_mask(z_index, roi, tform)
	return all_indices[roi_mask]
end

function get_roi_existing_mask(z_index::Int64, roi::BoundingBox=IMPORT_ROI_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	import_table = load_import_table(z_index)
	return get_roi_mask(z_index, roi, tform) & get_import_src_isfile(import_table)
end

function get_roi_existing_indices(z_index::Int64, roi::BoundingBox=IMPORT_ROI_BB, tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	import_table = load_import_table(z_index)
	roi_mask = get_roi_mask(z_index, roi, tform)
	return collect(1:size(import_table,1))[roi_mask]
end

# function compile_adjusted_tile_based_stats(z_index, N=100)
# 	println("Calculating statistics on original tiles for section $(z_index) using $N samples")
# 	import_table = load_import_table(z_index)
# 	bins = 20

# 	bias = load_bias_image()
# 	bias_mean = mean(bias)
# 	dtype = Float64
# 	# contrast_hist = load_contrast_histogram(z_index, reset=reset)
# 	# minval = contrast_hist[2]
# 	# maxval = contrast_hist[length(contrast_hist)-1]
# 	maxval = dtype(typemax(UInt16))
# 	minintensity = 39000 / maxval
# 	maxintensity = 63000 / maxval
# 	clahe = generate_clahe()

# 	function calculate_stats(i)
# 		index = get_import_index(import_table, i)
# 		println("Calculating for $index")
# 		src_fn = get_src_fn(import_table,i)
# 		if isfile(src_fn)
# 			img = get_image_disk(src_fn, dtype)
# 			img = bias_correction(img, bias, bias_mean)
# 			img = contrast_stretch(img, minintensity, maxintensity)
# 			img = convert_float64_to_uint8(img)
# 			# img = apply_clahe(img, clahe)
# 			hist = nquantile(img[:], bins)
# 			mu = mean(img[:])
# 			stdev = std(img[:])
# 			krt = kurtosis(img[:])
# 			return [index[3:4]..., mu, stdev, krt, hist...]
# 		end
# 		return [index[3:4]..., 0, 0, 0, zeros(bins+1)...]
# 	end

# 	existing_indices = collect(1:size(import_table,1))[get_import_src_isfile(import_table)]
# 	samples = rand(existing_indices, N)
# 	stats = pmap(calculate_stats, samples) 
# 	stats = hcat(stats...)'
# 	stats_fn = get_path("stats", premontaged(1, z_index))
# 	writedlm(stats_fn, stats)
# 	return stats
# end

# function compile_raw_tile_based_stats(z_index, N=300)
# 	println("Calculating statistics on original tiles for section $(z_index) using $N samples")
# 	import_table = load_import_table(z_index)
# 	bins = 20

# 	function calculate_stats(i)
# 		index = get_import_index(import_table, i)
# 		println("Calculating for $index")
# 		fn = get_src_fn(import_table,i)
# 		if isfile(fn)
# 			img = get_image_disk(fn, Float64)
# 			hist = nquantile(img[:], bins)
# 			mu = mean(img[:])
# 			stdev = std(img[:])
# 			krt = kurtosis(img[:])
# 			return [index[3:4]..., mu, stdev, krt, hist...]
# 		end
# 		return [index[3:4]..., 0, 0, 0, zeros(bins+1)...]
# 	end

# 	existing_indices = collect(1:size(import_table,1))[get_import_src_isfile(import_table)]
# 	samples = rand(existing_indices, N)
# 	stats = pmap(calculate_stats, samples) 
# 	stats = hcat(stats...)'
# 	stats_fn = get_path("stats", premontaged(1, z_index))
# 	writedlm(stats_fn, stats)
# 	return stats
# end

# function is_adjusted_resin(img, n=20, z_index=0)
# 	threshold = 100
# 	if z_index == 3
# 		threshold = 50
# 	end
# 	dist = nquantile(img[:], n)
# 	return dist[2] > threshold
# end

# function is_original_resin(img, n=20, z_index=0)
# 	thresholds = [0.79, 0.74, 0.75, 0.69, 0.75, 0.75]
# 	dist = nquantile(img[:], n)
# 	return dist[2] > thresholds[z_index]
# end

# function load_raw_tiles(z_index, file_indices; include_resin=false)
# 	import_table = load_import_table(z_index)

# 	function load_image(i)
# 		index = get_import_index(import_table, i)
# 		println("Loading $index")
# 		fn = get_src_fn(import_table, i)
# 		if isfile(fn)
# 			img = get_image_disk(fn, Float64)
# 			return index, img, !is_original_resin(img, 20, z_index), true
# 		end
# 		return index, Array{Float64,2}(), false, false
# 	end

# 	results = map(load_image, file_indices)
# 	indices = [i[1] for i in results]
# 	tiles = [i[2] for i in results]
# 	has_tissue = [i[3] for i in results]
# 	tile_exists = [i[4] for i in results]
# 	mask = include_resin ? has_tissue : tile_exists
# 	return indices[mask], tiles[mask]
# end

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
