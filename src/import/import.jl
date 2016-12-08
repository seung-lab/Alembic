global OVERVIEW_RESOLUTION = 95.3/3840 # 3.58/225.0 
# global LOCAL_RAW_DIR = joinpath(homedir(), "raw")
global LOCAL_RAW_DIR = joinpath(BUCKET, DATASET, "raw")
global GCLOUD_RAW_DIR = "gs://243774_8973/"
global GCLOUD_BUCKET = "gs://seunglab_alembic/datasets/"

# function gsutil_download(remote_file::AbstractString, local_file::Union{AbstractString, IO, Void}=nothing)
#    download_cmd = `gsutil -m cp
#        "$(bucket.provider.prefix)/$(bucket.name)/$remote_file" -`
#    s3_output = Pipe()
#    # open the cmd in write mode. this automatically takes the 2nd arg
#    # (stdio) and uses it as redirection of STDOUT.
#    (s3_input, process) = open(download_cmd, "w", s3_output)
#    close(s3_output.in)
# end

function get_loadfile()
	loadfile_sub_path = joinpath(DATASET, "161201_aibs_import.csv")
	loadfile_localpath = joinpath(BUCKET, loadfile_sub_path)
	if !isfile(loadfile_localpath)
		loadfile_remotepath = joinpath(GCLOUD_BUCKET, loadfile_sub_path)
		Base.run(`sudo gsutil -m cp $loadfile_remotepath $loadfile_localpath`)
	end
	return readdlm(loadfile_localpath, ',')
end

function get_src_dir(z_index)
	# LOADFILE = readdlm("/media/tmacrina/667FB0797A5072D7/3D_align/mosaiced_images_160729_google_cloud_upload_seung_import.csv",',')
	# LOADFILE = readdlm(joinpath(homedir(), "seungmount/research/Julimaps/datasets/pinky/import_aibs.csv"),',')
	loadfile = get_loadfile()
	i = findfirst(i -> i == z_index, loadfile[:,1])
	if i > 0
		return loadfile[i,2]
	end
end

function make_local_raw_dir()
	if !isdir(LOCAL_RAW_DIR)
		mkdir(LOCAL_RAW_DIR)
	end
end

function remove_premontaged_files(z_index)
	import_table = load_import_table(z_index)
	localpaths = get_local_tile_imported_paths(import_table)
	for path in localpaths
		if isfile(path)
			println(`rm $path`)
			Base.run(`rm $path`)
		end
	end
end

function remove_local_raw_dir()
	if isdir(LOCAL_RAW_DIR)
		Base.run(`rm -rf $LOCAL_RAW_DIR`)
	end
end

function get_local_raw_path(z_index)
	src_dir = get_src_dir(z_index)
	path = joinpath(LOCAL_RAW_DIR, src_dir)
	return path
end

function get_remote_raw_path(z_index)
	src_dir = get_src_dir(z_index)
	return joinpath(GCLOUD_RAW_DIR, src_dir)
end

function get_trakem_file(z_index)
	println("Downloading trakem file for $z_index")
	remote_raw_path = get_remote_raw_path(z_index)
	src = joinpath(remote_raw_path, "_trackem_\*")
	local_raw_path = get_local_raw_path(z_index)
	dst = joinpath(local_raw_path, "trakem_import.txt")
	Base.run(`gsutil -m cp $src $dst`)
	return readdlm(dst, '\t')
end

function download_subdir_files(z_index)
	dir = OVERVIEW_DIR
	subdirs = [CORRESPONDENCE_DIR]	
	for subdir in subdirs
		path = joinpath(dir, subdir)
		localpath = joinpath(BUCKET, DATASET, path)
		remotepath = joinpath(GCLOUD_BUCKET, DATASET, path)
		Base.run(`gsutil -m rsync -r $remotepath $localpath`)
	end
end

function upload_subdir_files(z_index)
	subdirs = ["thumbnail", "outline", "import", "contrast_bias", "contrast_stretch"]	
	for subdir in subdirs
		path = truncate_path(get_path(subdir, premontaged(1,z_index)))
		localpath = joinpath(BUCKET, DATASET, path)
		remotepath = joinpath(GCLOUD_BUCKET, DATASET, path)
		Base.run(`gsutil -m cp -r $localpath $remotepath`)
	end
end

function sync_to_upload()
	dir = PREMONTAGED_DIR
	println("Syncing subdirs for $dir")
	localpath = joinpath(BUCKET, DATASET, dir)
	remotepath = joinpath(GCLOUD_BUCKET, DATASET, dir)
	Base.run(`gsutil -m rsync -r $localpath $remotepath`)
end

function sync_to_download(z_index)
	dir = PREMONTAGED_DIR
	println("Syncing subdirs for $dir")
	localpath = joinpath(BUCKET, DATASET, dir)
	remotepath = joinpath(GCLOUD_BUCKET, DATASET, dir, "(1,$z_index,*")
	Base.run(`gsutil -m cp -r $remotepath $localpath`)
end

function download_raw_tiles(z_index; roi_only=true, overwrite=false)
	println("Downloading raw tiles for $z_index")
	import_table = load_import_table(z_index)
	remotepaths = get_remote_tile_raw_paths(import_table)
	localpaths = get_local_tile_raw_paths(import_table)
	if roi_only
		import_indices = get_included_indices(import_table)
		remotepaths = remotepaths[import_indices]
		localpaths = localpaths[import_indices]
	end
	if !overwrite
		tiles_do_not_exist_locally = [!isfile(f) for f in localpaths]
		remotepaths = remotepaths[tiles_do_not_exist_locally]
		localpaths = localpaths[tiles_do_not_exist_locally]
	end
	download_cmds = [`gsutil cp $src $dst` for (src, dst) in zip(remotepaths, localpaths)]
	pmap(Base.run, download_cmds)
end

function sync_premontage_registry()
	localpath = get_registry_path(premontaged(1,1))
	remotepath = joinpath(GCLOUD_BUCKET, DATASET, PREMONTAGED_DIR)
	Base.run(`gsutil cp $localpath $remotepath`)
end

function upload_imported_tiles(z_index)
	println("Uploading imported tiles for $z_index")
	import_table = load_import_table(z_index)
	localpaths = get_local_tile_imported_paths(import_table)
	remotepaths = get_remote_tile_imported_paths(import_table)
	tiles_exist_locally = [isfile(f) for f in localpaths]
	localpaths = localpaths[tiles_exist_locally]
	remotepaths = remotepaths[tiles_exist_locally]
	download_cmds = [`gsutil cp $src $dst` for (src, dst) in zip(localpaths, remotepaths)]
	pmap(Base.run, download_cmds)
end

function load_metadata(z_index)
	dir, _ = get_src_dir(z_index)
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
		trakem_table = get_trakem_file(z_index)
		initialize_import_table(trakem_table, z_index)
	end
	import_table = readdlm(import_fn, '\t')
	update_import_include!(z_index, import_table)
	return import_table
end

function initialize_import_table(trakem_table, z_index)
	println("Initializing import table for section $z_index")
	import_table = trakem_table

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

	dir = get_src_dir(z_index)
	import_table[:,1] = [joinpath(dir, fn) for fn in import_table[:,1]]
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
	maxval = 1.0
	minval = 1.0/255
	return min(maxval, max(minval, (img-minintensity) / (maxintensity-minintensity)))
end

function auto_contrast_stretch(img, bins=20)
	hist = nquantile(img[:], bins)
	minintensity = hist[2]
	maxintensity = hist[end-1]
	return contrast_stretch(img, minintensity, maxintensity)
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

	src_fn = get_local_tile_raw_path(import_table,i)
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

function clean_zeros(img, fill=0.5)
	return img[img .== 0] = fill
end

function find_zeros_coords(img)
	bl = []
	for i=1:size(img,1)
		for j=1:size(img,2)
			if img[i,j] == 0
				push!(bl, (i,j))
			end
		end
	end
	return bl
end

"""
AIBS tiles in pinky have known black pixel at (254,1123)
"""
function clean_aibs_image(img)
	sample = img[253:255, 1122:1124][:]
	fill = sum(sample) / (length(sample) - 1)
	img[254,1123] = fill
	return img
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
		img = get_image_disk(get_local_tile_raw_path(import_table, i), Float64)
		img = clean_aibs_image(img)
		Images.imfilter_gaussian_no_nans!(img, [10,10]) / N
	end

	# bias[bias .== 0] = 1 # not necessary with clean_aibs_image
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
		img = get_image_disk(get_local_tile_raw_path(import_table, i), dtype)
		img = clean_aibs_image(img)
		smp = bias_correction(img, bias, bias_mean)[rand(1:length(img), 14745)] #[3840*1900:3840*1920] # sample the middle of the tile to be faster
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

function import_overview_gcloud(z_index)
	dir = get_src_dir(z_index)
	if dir != nothing
		index = (1,z_index,OVERVIEW_INDEX,OVERVIEW_INDEX)
		dst_fn = string(get_path(index)[1:end-3], ".tif")
		src_fn = joinpath(dir, "_montage\*")
		f = "sudo gsutil -m cp $src_fn $dst_fn"
		return f
	end
	return nothing
end

function import_overview_gcloud(firstz, lastz)
	paths = []
	for z = firstz:lastz
		a = import_overview_gcloud(z)
		if a != nothing
			push!(paths, a)
		end
	end
	return paths
end

function import_overview(z_index)
	dir, _ = get_src_dir(z_index)
	if dir != nothing
		files = readdir(dir)
		k = findfirst(i -> contains(i, "montage"), files)
		index = (1,z_index,OVERVIEW_INDEX,OVERVIEW_INDEX)
		if k > 0
			src_fn = joinpath(dir, files[k])
			dst_fn = string(get_path(index)[1:end-3], ".tif")
			f = `cp $src_fn $dst_fn`
			println(f)
			Base.run(f)
		end
	end
end

function save_import_tile(fn, img)
	println("dst: $fn")
	f = h5open(fn, "w")
	chunksize = 500
	@time f["img", "chunk", (chunksize,chunksize)] = img
	close(f)	
end

function fix_contrast(src_index::Index, ref_index::Index)
	z_index = src_index[2]
	import_table = load_import_table(z_index)
	indices = get_import_indices(import_table)
	src_k = findfirst(i -> i == src_index, indices)
	ref_k = findfirst(i -> i == ref_index, indices)
	import_indices = get_included_indices(import_table)
	assert(src_k in import_indices)
	assert(ref_k in import_indices)
	println("Fixing $src_index contrast using $ref_index")
	bias = load_bias_image(z_index)
	bias_mean = mean(bias)

	ref_fn = get_local_tile_raw_path(import_table, ref_k)
	src_fn = get_local_tile_raw_path(import_table, src_k)

	dtype = Float64
	ref_img = get_image_disk(ref_fn, dtype)
	ref_img = bias_correction(ref_img, bias, bias_mean)
	src_img = get_image_disk(src_fn, dtype)
	src_img = bias_correction(src_img, bias, bias_mean)

	bins = 20
	hist = nquantile(ref_img[:], bins)
	minintensity = hist[2]
	maxintensity = hist[end-1]
	src_img = min(1.0, max(0.0, (src_img-minintensity) / (maxintensity-minintensity)))
	src_img = convert_float64_to_uint8(src_img)

	dst_fn = get_local_tile_dst_path(import_table, src_k)
	save_import_tile(dst_fn, src_img)
	return src_img
end

function premontage_cluster(z_range::UnitRange{Int64})
	loadfile = get_loadfile()
	gentrify_list = loadfile[:,1]
	pr = []
	for z in z_range
		try
			if z in gentrify_list
				gentrify_tiles(z)
			else 
				premontage_cluster(z)
			end
			push!(pr, [z,1])
		catch
			push!(pr, [z,0])
		end
		writedlm(joinpath(homedir(), "import_issues.txt"), pr)
	end
	sync_to_upload()
end

function premontage_cluster(z_index::Int64)
	sync_to_download(z_index)
	initialize_offsets(z_index)
	premontage(premontaged(1,z_index))
	remove_premontaged_files(z_index)
end

function gentrify_tiles(z_range::UnitRange{Int64})
	pr = []
	for z in z_range
		try 
			gentrify_tiles(z)
			push!(pr, [z,1])
		catch
			push!(pr, [z,0])
		end
		writedlm(joinpath(homedir(), "import_issues.txt"), pr)
	end
end

function gentrify_tiles(z_index::Int64)
	download_subdir_files(z_index)
	make_local_raw_dir()
	download_raw_tiles(z_index)
	import_tiles(z_index)
	premontage(premontaged(1,z_index))
	sync_to_upload()
	remove_local_raw_dir()
	remove_premontaged_files(z_index)
end

function import_tiles(z_index; from_current=false, reset=false, overwrite_offsets=true)
	thumbnail_scale = OVERVIEW_RESOLUTION
	import_table = load_import_table(z_index)
	N = size(import_table, 1)

	import_indices = get_included_indices(import_table)
	if from_current
		import_indices = get_unwritten_included_indices(import_table)
	end
	n = length(import_indices)
	bins = 20
	# contrast_clusters = get_contrast_clusters(import_table)
	bias = load_bias_image(z_index, reset=reset)
	bias_mean = mean(bias)
	contrast_hist = load_contrast_histogram(z_index, reset=reset)
	minintensity = contrast_hist[2]
	maxintensity = contrast_hist[end-1]
	dtype = Float64

	function import_tile(i)
		src_fn = get_local_tile_raw_path(import_table, i)
		dst_fn = get_local_tile_dst_path(import_table, i)
		index = get_import_index(import_table, i)
		println("Importing $index, # $i / $n")
		println("\tsrc: $src_fn")

		img = get_image_disk(src_fn, dtype)
		sz = size(img)
		offset = get_import_offset(import_table, i)
		img = bias_correction(img, bias, bias_mean)
		hist = nquantile(img[:], bins)
		# if hist[2] < minintensity_threshold
		# 	img = contrast_stretch(img, hist[2], hist[end-1])
		# else
		img = contrast_stretch(img, minintensity, maxintensity)
		# end
		img = convert_float64_to_uint8(img)
		# resin = is_adjusted_resin(img, 20, z_index)
		save_import_tile(dst_fn, img)
		return i, sz..., offset, img
	end

	results = pmap(import_tile, import_indices)

	offset = [i[4] for i in results]
	tiles = [i[5] for i in results]
	println("Creating imported stage stitched thumbnail")
	ss, ss_offset = merge_images(tiles, offset)
	thumbnail, thumbnail_offset = imscale(ss, thumbnail_scale)
    write_thumbnail(thumbnail, premontaged(1,z_index), thumbnail_scale)

	if overwrite_offsets
		indices = [i[1] for i in results]
		order = sortperm(indices)
		height = [i[2] for i in results][order]
		width = [i[3] for i in results][order]

		println("Appending image size to import table")
		import_table[import_indices,9] = height
		import_table[import_indices,10] = width

		import_dst_fn = get_path("import", premontaged(1, z_index))
		println("Writing import table for section $z_index to $import_dst_fn")
		writedlm(import_dst_fn, import_table)
		initialize_offsets(z_index)
	end
end

function downsample_tiles(z_index, thumbnail_scale=0.25)
	import_table = load_import_table(z_index)
	N = size(import_table, 1)
	existing_indices = create_src_isfile_mask(import_table)
	indices = collect(1:N)[existing_indices]
	n = length(indices)
	dtype = Float64

	function downsample_tile(i)
		src_fn = get_local_tile_raw_path(import_table, i)
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

function get_remote_tile_raw_path(import_table, i)
	return joinpath(GCLOUD_RAW_DIR, import_table[i,1])
end

function get_local_tile_raw_path(import_table, i)
	return joinpath(LOCAL_RAW_DIR, import_table[i,1])
end

function get_remote_tile_imported_path(import_table, i)
	index = get_import_index(import_table, i)
	return joinpath(GCLOUD_BUCKET, get_path(index))
end

function get_local_tile_dst_path(import_table, i)
	index = get_import_index(import_table, i)
	return get_path(index)
end

function get_local_tile_raw_paths(import_table)
	return [get_local_tile_raw_path(import_table, i) for i in 1:size(import_table,1)]
end

function get_remote_tile_imported_paths(import_table)
	return [get_remote_tile_imported_path(import_table, i) for i in 1:size(import_table,1)]
end

function get_local_tile_imported_paths(import_table)
	return [get_local_tile_dst_path(import_table, i) for i in 1:size(import_table,1)]
end

function get_remote_tile_raw_paths(import_table)
	return [get_remote_tile_raw_path(import_table, i) for i in 1:size(import_table,1)]
end

function get_dst_filenames(import_table)
	return [get_local_tile_dst_path(import_table, i) for i in 1:size(import_table,1)]
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
	return convert(BitArray, map(isfile, get_local_tile_raw_paths(import_table)))
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

function initialize_offsets(z_index; )
	import_table = load_import_table(z_index)
	included_tiles = get_included_indices(import_table)
	indices = get_import_indices(import_table)
	offsets = get_import_offsets(import_table)
	# sizes = get_import_sizes(import_table)
	sizes = [(3840,3840) for i in 1:size(import_table,1)]
	update_offsets(indices[included_tiles], offsets[included_tiles], sizes[included_tiles])
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
	return [bb_to_pts(bb) for bb in bbs]
end

function tform_bbs_pts(bbs_pts, tform)
	tform_pts = [[pts ones(size(pts,1),1)]*tform for pts in bbs_pts]
	return [pts[:, 1:2] for pts in tform_pts]
end

function snap_bbs_pts(bbs_pts)
	return [round(Int64, pts) for pts in bbs_pts]
end

function get_tform(index::Index, scale=1.0)
	tform = eye(3)
	tform_path = get_path("cumulative_transform", index)
	if isfile(tform_path)
		tform = load("cumulative_transform", index)
	end
	s = make_scale_matrix(scale)
	return s*tform*s^-1
end

function view_import_bb(z_index::Int64, tform=eye(3))
	import_table = load_import_table(z_index)
	polys = make_import_outline(import_table, tform)
	# polys = [poly[1:4,:] for poly in polys]
	polys = tform_bbs_pts(polys, make_scale_matrix(0.05))
	indices = get_import_indices(import_table)
	view_polys(polys, indices)
end

function get_roi(z_index::Int64; scale=1.0)
	return load("correspondence", overview(1,z_index)) * scale
end

function view_roi(z_index::Int64; tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	s = 0.05
	import_table = load_import_table(z_index)
	polys = make_import_outline(import_table, tform)
	# polys = [poly[1:4,:] for poly in polys]
	polys = tform_bbs_pts(polys, make_scale_matrix(s))
	indices = get_import_indices(import_table)
	roi = get_roi(z_index, scale=s/OVERVIEW_RESOLUTION)
	view_polys(polys, indices, roi)
end

function make_import_outline(import_table, tform=eye(3))
	bbs = get_bbs(import_table)
	pts = get_bbs_pts(bbs)
	return tform_bbs_pts(pts, tform)
end

function get_roi_mask(import_table, roi, tform)
	# possible speed up: run bounding box intersection first
	# all_bbs = get_tform_bbs(z_index)
	# roi_intersects = convert(BitArray, [intersects(bb, roi) for bb in all_bbs])
	# bbs = all_bbs[roi_intersects]
	# indices = all_indices[roi_intersects]
	tiles = make_import_outline(import_table, tform)
	return convert(BitArray, [poly_intersects(tile, roi) for tile in tiles])
end

function get_tile_indices(z_index::Int64, roi=get_roi(z_index, scale=1/OVERVIEW_RESOLUTION), tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	import_table = load_import_table(z_index)
	all_indices = get_import_indices(import_table)
	roi_mask = get_roi_mask(import_table, roi, tform)
	return all_indices[roi_mask]
end

function update_import_include!(z_index, import_table, roi=get_roi(z_index, scale=1/OVERVIEW_RESOLUTION), tform=get_tform(overview(1,z_index), OVERVIEW_RESOLUTION))
	include_mask = create_flagged_mask(import_table)
	roi_mask = get_roi_mask(import_table, roi, tform)
	# exisiting_mask = create_src_isfile_mask(import_table)
	import_table[:,12] = include_mask & roi_mask # & exisiting_mask
	# save_import_table(z_index, import_table)
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
	unwritten = !create_dst_isfile_mask(import_table)
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
			Base.run(`mv $src_fn $dst_fn`)
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

function expunge_imported_tiles_outside_roi(z_index::Int64)
	outside = outside_import_roi(z_index)
	exists = convert(BitArray, [isfile(get_path(index)) for index in outside])
	expunge_tiles(outside[exists])
end

function expunge_tiles_not_included(z_index::Int64)
	import_table = load_import_table(z_index)
	indices = get_import_indices(import_table)
	not_included = !create_included_mask(import_table)
	for index in indices[not_included]
		if isfile(get_path(index))
			expunge_tile(index)
		end
	end
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
	println("$index: $tform")
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

function copy_overviews_to_montages(z_range)
	for z = z_range
		src_index = overview(1,z)
		dst_index = montaged(1,z)
		src_fn = get_path(src_index)
		dst_fn = get_path(dst_index)
		if isfile(src_fn)
			f = `cp $src_fn $dst_fn`
			# println(f)
			Base.run(f)
			# println(dst_index, get_overview_size(dst_index))
			update_registry(dst_index, offset=[0,0], image_size=get_overview_size(dst_index))
		end
	end
end

function copy_montage_import_folders_to_overview_correspondences(z_range)
	for z in z_range
		println("Copying montaged(1,$z) import polygon")
		pts = load("import", montaged(1,z)) 
		correspondences_path = get_path("correspondence", overview(1,z))
		writedlm(correspondences_path, pts)
	end
end

function copy_between_projects(z_range, stage, dir_name, src_dataset, dst_dataset)
	gp = get_path(dir_name, eval(stage)(1,z_range[1]))
	datasets_path = join(split(gp, "/")[1:end-4], "/")
	src_path = joinpath(datasets_path, src_dataset)
	dst_path = joinpath(datasets_path, dst_dataset)
	for z in z_range
		src_subdir = get_path(dir_name, eval(stage)(1,z))
		src_subdir = join(split(src_subdir, "/")[end-2:end], "/")
		src = joinpath(src_path, src_subdir)
		dst_subdir = get_path(dir_name, overview(1,z))
		dst_subdir = join(split(dst_subdir, "/")[end-2:end], "/")
		dst = joinpath(dst_path, dst_subdir)
		if isfile(src)
			f = `cp $src $dst`
			println(f)
			Base.run(f)
		else
			println("Does not exist: $src")
		end
	end
end

function copy_pinky_overviews_from_gcloud(fn=joinpath(homedir(), "seungmount/research/Alembic/datasets/pinky_overviews/161206_import_reimaged.csv"))
	import_table = readdlm(fn, ',')
	src_dir = GCLOUD_RAW_DIR
	dst_dir = joinpath(homedir(), "seungmount/research/Alembic/datasets/pinky_overviews/0_overview/second_reimaged_sections")
	for i = 1:size(import_table,1)
		z = import_table[i,1]
		src_fn = joinpath(src_dir, import_table[i,2], "_montage\*")
		dst_fn = joinpath(dst_dir, "1,$(z)_overview.tif")
		f = `sudo gsutil -m cp $src_fn $dst_fn`
		println(f)
		Base.run(f)
	end
end

function resave_overview(z_index)
	index = (1,z_index,OVERVIEW_INDEX,OVERVIEW_INDEX)
	dtype = UInt8
	src_fn = string(get_path(index)[1:end-3], ".tif")
	img = get_image_disk(src_fn, dtype)
	dst_fn = get_path(index)
	save_import_tile(dst_fn, img)
end

function download_one_tile_aibs_alignment(r, c, z_range=3279:3528)
	GCLOUD_RAW_DIR = "gs://aibs_alignment/20161111_output_TS5"
	dst_dir = "/usr/people/tmacrina/seungmount/research/Alembic/datasets/AIBS_pinky_test/"
	dst_dir = joinpath(dst_dir, "$(r)_$(c)")
	if !isdir(dst_dir)
		mkdir(dst_dir)
	end
	for z in z_range
		src = joinpath(GCLOUD_RAW_DIR, "$z", "0", "$r", "$(c).png")
		dst = joinpath(dst_dir, "$(r)_$(c)_$(z).png")
		f = `sudo gsutil -m cp -r $src $dst`
		println(f)
		Base.run(f)
	end
end
