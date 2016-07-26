function get_src_dir(z_index)
	return joinpath("/media/tmacrina/Data/AIBS_actual_trial/Dropbox/", "$z_index")
end

function initialize_import_table(import_src_fn, z_index)
	println("Initializing import table for section $z_index from $import_src_fn")
	import_table = readdlm(import_src_fn)
	import_dst_fn = get_path("import", premontaged(1, z_index+1))

	# pull out row & col from filename
	path_no_ext = [split(i, ".")[1] for i in import_table[:,1]]
	waf = [1 for i in path_no_ext]
	col = [parse(split(i, "_")[9])+1 for i in path_no_ext]
	row = [parse(split(i, "_")[10])+1 for i in path_no_ext]
	sec = [z_index+1 for i in path_no_ext]
	# include full path to the file
	dir = get_src_dir(z_index)
	import_table[:,1] = [joinpath(dir, fn) for fn in import_table[:,1]]
	file_exists = map(isfile, import_table[:,1])
	dst_path = map(get_path, zip(waf, sec, row, col))
	import_table = hcat(import_table, [waf sec row col dst_path file_exists])
	println("Writing import table for section $z_index to $import_dst_fn")
	writedlm(import_dst_fn, import_table)
end

function compile_tile_based_stats(z_index, N=60)
	println("Calculating statistics on original tiles for section $(z_index+1) using $N samples")
	import_table = load_import_table(z_index)

	function calculate_stats(i)
		index = get_import_index(import_table, i)
		println("Calculating for $index")
		fn = get_src_fn(import_table,i)
		if isfile(fn)
			img = get_image_disk(fn, Float64)
			hist = nquantile(img[:], 20)
			krt = kurtosis(img[:])
			return [index[3:4]..., krt, hist...]
		end
		return [index[3:4]..., 0, zeros(21)...]
	end

	samples = rand(1:size(import_table,1), N)
	stats = pmap(calculate_stats, samples) 
	stats = hcat(stats...)'
	stats_fn = get_path("stats", premontaged(1, z_index+1))
	writedlm(stats_fn, stats)
	return stats
end

function is_resin(img, n=20, z_index=0)
	thresholds = [0.79, 0.74, 0.75]
	dist = nquantile(img[:], n)
	return dist[2] > thresholds[z_index+1]
end

function load_raw_tiles(z_index, file_indices; include_resin=false)
	import_table = load_import_table(z_index)

	function load_image(i)
		index = get_import_index(import_table, i)
		println("Loading $index")
		fn = get_src_fn(import_table, i)
		if isfile(fn)
			img = get_image_disk(fn, Float64)
			return index, img, !is_resin(img, 20, z_index), true
		end
		return index, Array{Float64,2}(), false, false
	end

	results = pmap(load_image, file_indices)
	indices = [i[1] for i in results]
	tiles = [i[2] for i in results]
	has_tissue = [i[3] for i in results]
	tile_exists = [i[4] for i in results]
	mask = include_resin ? has_tissue : tile_exists
	return indices[mask], tiles[mask]
end

function calculate_bias_field(z_index, N=100)
	println("Calculating contrast bias field for section $(z_index+1) using $N samples")
	import_table = load_import_table(z_index)	
	samples = rand(1:size(import_table,1), N)
	indices, tiles = load_raw_tiles(z_index, samples, include_resin=false)
	bias = sum(tiles) /length(tiles)
	bias[bias .== 0] = 1e-5
	bias_fn = get_path("contrast_bias", premontaged(1, z_index+1))
	f = h5open(bias_fn, "w")
	f["img"] = bias
	close(f)
	println("Contrast bias image written to $bias_fn")
end

function calculate_contrast_histogram(z_index, N=60)
	println("Calculating overall contrast histogram for section $(z_index+1) using $N samples")
	import_table = load_import_table(z_index)	
	bias = load_bias_image(z_index)
	bias_mean = mean(bias)
	samples = rand(1:size(import_table,1), N)
	indices, tiles = load_raw_tiles(z_index, samples, include_resin=false)

	function adjust_images(img)
		img = (img ./ bias) * bias_mean 
		return img[:]
	end

	intensities = pmap(adjust_images, tiles)
	intensities = vcat(intensities...)
	distr = nquantile(intensities, 20)
	auto_fn = get_path("contrast_auto", premontaged(1, z_index+1))
	writedlm(auto_fn, distr)
	println("Contrast histogram written to $auto_fn")
end

function load_import_table(z_index; reset=false)
	import_fn = get_path("import", premontaged(1, z_index+1))
	if !isfile(import_fn) || reset
		dir = get_src_dir(z_index)
		files = readdir(dir)
		trakem_fn = files[findfirst(i -> contains(i, "trackem"), files)]
		initialize_import_table(joinpath(dir, trakem_fn), z_index)
	end
	return readdlm(import_fn)
end

function load_bias_image(z_index; reset=false)
	bias_fn = get_path("contrast_bias", premontaged(1, z_index+1))
	if !isfile(bias_fn) || reset
		calculate_bias_field(z_index)
	end
	return h5read(bias_fn, "img")
end

function load_contrast_histogram(z_index; reset=false)
	auto_fn = get_path("contrast_auto", premontaged(1, z_index+1))
	if !isfile(auto_fn) || reset
		calculate_contrast_histogram(z_index)
	end
	return readdlm(auto_fn)
end

function import_tiles(z_index; reset=false)
	import_table = load_import_table(z_index, reset=reset)
	N = size(import_table, 1)
	bias = load_bias_image(z_index, reset=reset)
	bias_mean = mean(bias)
	contrast_hist = load_contrast_histogram(z_index, reset=reset)
	minval = contrast_hist[2]
	maxval = contrast_hist[length(contrast_hist)-1]

	function bias_correction(img)
		return (img ./ bias) * bias_mean
	end

	function auto_adjust(img) #, n=20)
		return min(1, max(0, (img-minval) / (maxval-minval)))
	end

	function convert_float64_to_uint8(img)
		return Array{UInt8, 2}(round(UInt8, img * 255))
	end

	function import_tile(i)
		src_fn = get_src_fn(import_table, i)
		dst_fn = get_dst_fn(import_table, i)
		index = get_import_index(import_table, i)
		println("Importing $index, # $i / $N")
		println("\tsrc: $src_fn")

		img = get_image_disk(src_fn, Float64)
		sz = size(img)
		resin = is_resin(img, 20, z_index)
		img = bias_correction(img)
		img = auto_adjust(img)
		img = convert_float64_to_uint8(img)
		if !resin
			println("\tdst: $dst_fn")
			f = h5open(dst_fn, "w")
			chunksize = 1000
			@time f["img", "chunk", (chunksize,chunksize)] = img
			close(f)
		end
		return i, sz..., resin
	end

	tile_indices = 1:N
	tile_exists = get_import_src_isfile(import_table)
	existing_tile_indices = tile_indices[tile_exists]
	results = pmap(import_tile, existing_tile_indices)

	indices = [i[1] for i in results]
	height = [i[2] for i in results]
	width = [i[3] for i in results]
	resin = [i[4] for i in results]

	println("Appending image size & resin test results to import table")
	import_table = hcat(import_table, [zeros(Int64, N) zeros(Int64, N) falses(N)])
	import_table[existing_tile_indices,11] = height
	import_table[existing_tile_indices,12] = width
	import_table[existing_tile_indices,13] = resin

	import_dst_fn = get_path("import", premontaged(1, z_index+1))
	println("Writing import table for section $z_index to $import_dst_fn")
	writedlm(import_dst_fn, import_table)
	initialize_offsets(z_index)
end

function get_src_fn(import_table, i)
	return import_table[i,1]
end

function get_dst_fn(import_table, i)
	return import_table[i,9]
end

function get_dst_filenames(import_table)
	return [get_dst_fn(import_table, i) for i in 1:size(import_table,1)]
end

function get_import_index(import_table, i)
	return (import_table[i,5:8]...)
end

function get_import_indices(import_table)
	return [get_import_index(import_table, i) for i in 1:size(import_table,1)]
end

function get_import_offsets(import_table)
	return [[import_table[i,3], import_table[i,2]] for i in 1:size(import_table,1)]
end

function get_import_dst_isfile(import_table)
	return convert(BitArray, map(isfile, get_dst_filenames(import_table)))
end

function get_import_src_isfile(import_table)
	return convert(BitArray, import_table[:,10])
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
