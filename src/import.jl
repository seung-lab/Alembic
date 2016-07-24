function calculate_sample_stats(k, N=60)
	
	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	function calculate_stats(i, dir, meta, N)
		println("# $i / ", N)
		fn = joinpath(dir, meta[i,1])
		row, col = meta[i,6], meta[i,7]
		if isfile(fn)
			img = read(fn)
			hist = nquantile(img[:], 20)
			krt = kurtosis(img[:])
			return [row, col, krt, hist...]
		end
		return [row, col, 0, zeros(21)...]
	end

	dir = joinpath("/media/tmacrina/Data/AIBS_actual_trial/Dropbox (MIT)/", "$k")
	files = readdir(dir)
	j = findfirst(i -> contains(i, "trackem"), files)
	stats = []
	if j > 0
		trakem_fn = files[j]
		suffix = string("_sec", @sprintf("%05d", k))
		meta_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))
		if !isfile(meta_fn)
			build_AIBS_actual_meta(joinpath(dir, trakem_fn), k)
		end
		meta = readdlm(meta_fn)

		samples = rand(1:size(meta,1), N)
		stats = pmap(calculate_stats, samples, repeated(dir), repeated(meta), repeated(N)) 
		stats = hcat(stats...)'
		stats_fn = joinpath(RAW_DIR, string("stats", suffix, ".txt"))
		writedlm(stats_fn, stats)
		return stats
	end
end

function is_AIBS_actual_resin(img)
	dist = nquantile(img[:], 20)
	return dist[2] > 0.75
end

function build_AIBS_actual_bias(meta, dir, suffix="")
	
	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	function load_images(i, dir, meta)
		println("Compiling meta $i")
		fn = joinpath(dir, meta[i,1])
		if isfile(fn)
			img = read(fn)
			return img, !is_AIBS_actual_resin(img)
		end
		return Array{Float64,2}(), false
	end

	println("Constructing bias adjustment image")
	samples = rand(1:size(meta,1), 50)
	results = pmap(load_images, samples, repeated(dir), repeated(meta))
	imgs = [i[1] for i in results]
	tissue = [i[2] for i in results]
	bias = sum(imgs[tissue .== true]) / sum(tissue)
	bias[bias .== 0] = 1e-5
	f = h5open(joinpath(RAW_DIR, string("bias", suffix, ".h5")), "w")
	f["img"] = bias
	close(f)
end

function build_AIBS_actual_adjustment(meta, dir, suffix, n=20)
	println("Constructing bias adjusted histogram")

	adj_fn = joinpath(RAW_DIR, string("adj", suffix, ".txt"))
	intensities = Array{Float64,1}()

	bias_fn = joinpath(RAW_DIR, string("bias", suffix, ".h5"))
	if !isfile(bias_fn)
		build_AIBS_actual_bias(meta, dir, suffix)
	end
	bias = h5read(bias_fn, "img")
	bias_mean = mean(bias)

	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	function compile_images(i, dir, meta, bias, bias_mean)
		println("Compiling meta $i")
		fn = joinpath(dir, meta[i,1])
		if isfile(fn)
			img = read(fn)
			if !is_AIBS_actual_resin(img)
				img = (img ./ bias) * bias_mean
				return img[:]
			end
		end
		return Array{Float64,1}()
	end

	samples = rand(1:size(meta,1), 30)
	intensities = pmap(compile_images, samples, repeated(dir), repeated(meta), repeated(bias), repeated(bias_mean)) 
	println("Calculating intensity histogram...")
	distr = nquantile(vcat(intensities...), n)
	println(distr)
	writedlm(adj_fn, distr)
end

function build_AIBS_actual_meta(meta_fn, k)
	meta = readdlm(meta_fn)
	suffix = string("_sec", @sprintf("%05d", k))
	meta_new_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))

	# cleanup meta text file
	path = meta[:,1]
	path_no_ext = [split(i, ".")[1] for i in path]
	col = [parse(split(i, "_")[9])+1 for i in path_no_ext]
	row = [parse(split(i, "_")[10])+1 for i in path_no_ext]
	sec = [k+1 for i in path_no_ext]
	# string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
	# path_new_fn = [string("Tile_r", @sprintf("%02d", i), "-c", @sprintf("%02d", j), "_S2-W001_sec", @sprintf("%04d", k), ".h5") for (i,j,k) in zip(row, col, sec)]
	path_new_fn = [string("Tile_r", r, "-c", c, "_S2-W001_sec", s, ".h5") for (r,c,s) in zip(row, col, sec)]
	meta = hcat(meta, [path_new_fn row col sec])
	writedlm(meta_new_fn, meta)
end

function import_AIBS_actual(k; reset=false)
	dir = joinpath("/media/tmacrina/Data/AIBS_actual_trial/Dropbox (MIT)/", "$k")
	files = readdir(dir)
	j = findfirst(i -> contains(i, "trackem"), files)
	if j > 0
		trakem_fn = files[j]
		suffix = string("_sec", @sprintf("%05d", k))
		meta_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))
		if !isfile(meta_fn) || reset
			build_AIBS_actual_meta(joinpath(dir, trakem_fn), k)
		end
		meta = readdlm(meta_fn)

		bias_fn = joinpath(RAW_DIR, string("bias", suffix, ".h5"))
		if !isfile(bias_fn) || reset
			build_AIBS_actual_bias(meta, dir, suffix)
		end
		bias = h5read(bias_fn, "img")
		bias_mean = mean(bias)

		adj_fn = joinpath(RAW_DIR, string("adj", suffix, ".txt"))
		if !isfile(adj_fn) || reset
			build_AIBS_actual_adjustment(meta, dir, suffix)
		end
		adj = readdlm(adj_fn)
		minval = adj[2]
		maxval = adj[length(adj)-1]

		function read(fn)
			return Array{Float64, 2}(Images.load(fn).data')
		end

		function bias_correction(img, bias, bias_mean)
			return (img ./ bias) * bias_mean
		end

		function auto_adjust(img, minval, maxval) #, n=20)
			# distr = nquantile(img[:], n)
			# minval = distr[2]; maxval = distr[n];
			return min(1, max(0, (img-minval) / (maxval-minval)))
		end

		function convert_float64_to_uint8(img)
			return Array{UInt8, 2}(round(UInt8, img * 255))
		end

		function correct_and_save(i, dir, meta, bias, bias_mean, minval, maxval)
			println("# $i / ", size(meta,1))
			old_fn = joinpath(dir, meta[i,1])
			new_fn = joinpath(RAW_DIR, meta[i,5])
			index = (1, meta[i,8], meta[i,6], meta[i,7])
			offset = [meta[i,3], meta[i,2]]

			println(old_fn)
			if isfile(old_fn) # && !isfile(new_fn)
				img = read(old_fn)
				if !is_AIBS_actual_resin(img)
					img = bias_correction(img, bias, bias_mean)
					img = auto_adjust(img, minval, maxval)
					img = convert_float64_to_uint8(img)
					#f = h5open(joinpath(to, import_frame_tif_to_wafer_section_h5(tif, wafer_num, sec_num)), "w")
					f = h5open(new_fn, "w")
					#    	chunksize = div(min(size(img)...), 4);
					chunksize = 1000;
					@time f["img", "chunk", (chunksize,chunksize)] = img
					close(f)
					update_offset(index, offset, [size(img)...])
					return i, false
				end
			end
			return i, true
		end

		is_resin = [true for i in 1:size(meta,1)]
		ignore = pmap(correct_and_save, 1:size(meta,1), repeated(dir), 
							repeated(meta), repeated(bias), repeated(bias_mean), 
							repeated(minval), repeated(maxval))
		indices = [i[1] for i in ignore]
		ignore = [i[2] for i in ignore]
		to_ignore = ignore[sortperm(indices)]
		meta = hcat(meta, to_ignore)
		writedlm(meta_fn, meta)
	end
end

function import_AIBS_actual_tile(k, row, col)
	suffix = string("_sec", @sprintf("%05d", k))
	meta_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))
	meta = readdlm(meta_fn)

end

function reset_offsets(k)
	dir = joinpath("/media/tmacrina/Data/AIBS_actual_trial/Dropbox (MIT)/", "$k")
	files = readdir(dir)
	j = findfirst(i -> contains(i, "trackem"), files)
	dist = []
	if j > 0
		trakem_fn = files[j]
		suffix = string("_sec", @sprintf("%05d", k))
		meta_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))
		if !isfile(meta_fn)
			build_AIBS_actual_meta(joinpath(dir, trakem_fn), k)
		end
		meta = readdlm(meta_fn)

		for i in 1:size(meta, 1)
			println("# $i / ", size(meta,1))
			old_fn = joinpath(dir, meta[i,1])
			new_fn = joinpath(RAW_DIR, meta[i,5])
			index = (1, meta[i,8], meta[i,6], meta[i,7])
			offset = [meta[i,3], meta[i,2]]
			is_resin = meta[i,8]

			println(old_fn)
			if isfile(old_fn) # && !isfile(new_fn)
				if is_resin
					update_offset(index, offset, [3840,3840])
				else
					println("Not importing - this is resin.")
				end
			else
				println("Does not exist!")
			end
		end
	end
end

# function import_AIBS_pilot_v2()
# 	trakem_fn = joinpath(homedir(), "seungmount/research/Julimaps/datasets/AIBS_pilot_v2/234251S6R_01_02_160421_zorder_sectmanifest_columnsappended.txt")
# 	image_dir = joinpath(homedir(), "seungmount/research/Julimaps/datasets/AIBS_pilot_v2/0_raw")

# 	meta = readdlm(trakem_fn)

# 	# cleanup trakem text file
# 	fn_length = map(length, meta[:,1])
# 	fn_start = 3
# 	fn_range = map(range, repeated(fn_start), fn_length-fn_start+1)
# 	meta[:,1] = map(getindex, meta[:,1], fn_range)

# 	# rename files based on ordering
# 	fn_splits = [split(split(i, ".")[1], "_")[end-1:end] for i in meta[:,1]]
# 	row = [parse(i[1]) for i in fn_splits]
# 	col = [parse(i[2]) for i in fn_splits]
# 	meta = hcat(meta, [row col])

# end

function import_trakem_dir_tifs(src_folder, suffix="_prealigned", dst_folder=src_folder)
	for tif in sort_dir(src_folder, "tif")
		println(tif)
		img = Array{Float64, 2}(Images.load(joinpath(src_folder, tif)).data')
		#	  img = img ./ NORMALIZER;
		#	  distr = nquantile(img[:], 16)
		#	  minval = distr[2]; maxval = distr[16];
		#	  img = min(1, max(0, (img - minval) / (maxval-minval)))
		img = Array{UInt8, 2}(round(UInt8, img * 255));
		#f = h5open(joinpath(to, import_frame_tif_to_wafer_section_h5(tif, wafer_num, sec_num)), "w")
		f = h5open(joinpath(dst_folder, import_tif_to_wafer_section_h5(tif, suffix)), "w")
		#    	chunksize = div(min(size(img)...), 4);
		chunksize = 1000;
		@time f["img", "chunk", (chunksize,chunksize)] = img
		close(f)
		update_offset(prealigned(1, int(tif[1:3])), [0, 0], [size(img)...])
	end
end

function import_AIBS_dir_meta(folder, wafer_num, sec_num)
        if !isfile(joinpath(folder, "meta.txt")) return nothing; end
        registry = readdlm(joinpath(folder, "meta.txt"))[:, 1:3];
	registry[:,1] = map(import_frame_tif_to_wafer_section_h5, registry[:, 1], repeated(wafer_num), repeated(sec_num));
	registry[:,2:3] = registry[:, 3:-1:2];
	registry = hcat(registry, fill(3840, size(registry, 1), 2));
	return registry
      end

function import_AIBS_grid(folder, wafer_num, to)
  	sec_folders = readdir(folder);
	fullpaths = map(joinpath, repeated(folder), sec_folders);
	pmap(import_AIBS_dir_tifs, fullpaths, repeated(wafer_num), collect(1:length(sec_folders)), repeated(to))
	regs = map(import_AIBS_dir_meta, fullpaths, repeated(wafer_num), collect(1:length(sec_folders)))
	regs = regs[find(reg -> reg != nothing, regs)]
	return vcat(regs...);
end

function import_frame_tif_to_wafer_section_h5(frame_name, wafer_num, sec_num)
	return string(wafer_num, "-", sec_num, "-", frame_name[7:9], ".h5")
end

function import_tif_to_wafer_section_h5(filename, suffix)
	return string("1", ",", int(filename[1:3]), suffix, ".h5")
end

function clean_up_nicks_data(src_folder, dst_folder)
	for fn in sort_dir(src_folder, ".h5")
		img = h5read(joinpath(src_folder, fn), "main")
		min_px = minimum(img)
		max_px = maximum(img)
		img = Array{UInt8, 2}(round(UInt8, (img-min_px)*(255/(max_px-min_px))));
		f = h5open(joinpath(dst_folder, fn), "w")
		chunksize = 100;
		@time f["img", "chunk", (chunksize,chunksize)] = img
		close(f)
		update_offset(parse_name(fn), [0, 0], [size(img)...])
	end
end

function import_tif_as_h5s_sample()

tuples = Array{Any}(0)
for i in 1:15, j in 1:200, k in 1:10, l in 1:10
       push!(tuples, (i,j,k,l))
       end
@parallel for tuple in tuples
       path = get_path(tuple, ".tif");
       if isfile(path) println(path);
           img = get_image(path);
           f = h5open(joinpath(homedir(), "datasets", "zfish", string(get_name(tuple), ".h5")), "w")
           @time f["img", "chunk", (1000, 1000)] = img
           close(f)
	   #run(`aws s3 cp path s3://seunglab/datasets/zfish_10-12/`)
       end
end
end


