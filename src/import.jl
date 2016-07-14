#NORMALIZER = Array{Float64, 2}(readdlm("calib"))
using StatsBase

function construct_adjustment_image(meta, suffix="")
	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	println("Constructing contrast adjustment image")
	adj_count = 1
	j = 50
	filepath = joinpath(RAW_DIR, meta[j,5])
	img_adj = read(filepath)
	for i in 92:90:size(meta,1)
		filepath = joinpath(RAW_DIR, meta[i,5])
		img_adj += read(filepath)
		adj_count += 1
	end
	img_adj /= adj_count
	fn = joinpath(RAW_DIR, "adjustment_image.txt")
	writedlm(fn, img_adj)
end

function import_AIBS_pilot_v1(construct=false)
	meta_fn = joinpath(RAW_DIR, "2342516R01_ribbonordered.txt")
	meta = readdlm(meta_fn)
	meta_new_fn = joinpath(RAW_DIR, "meta.txt")

	# cleanup meta text file
	path = meta[:,1]
	path_split = map(split, path, repeated("\\"))
	path_grid = [i[1] for i in path_split]
	path_subdir = [i[2] for i in path_split]
	path_old_fn = [i[3] for i in path_split]
	path = [joinpath(a,b,c) for (a,b,c) in zip(path_grid, path_subdir, path_old_fn)]
	path_old_fn_no_ext = [split(i, ".")[1] for i in path_old_fn]
	col = [parse(split(i, "-")[2])+1 for i in path_old_fn_no_ext]
	row = [parse(split(i, "-")[3])+1 for i in path_old_fn_no_ext]
	sec = meta[:,4]+1
	# string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
	# path_new_fn = [string("Tile_r", @sprintf("%02d", i), "-c", @sprintf("%02d", j), "_S2-W001_sec", @sprintf("%04d", k), ".h5") for (i,j,k) in zip(row, col, sec)]
	path_new_fn = [string("Tile_r",i, "-c", j, "_S2-W001_sec", k, ".h5") for (i,j,k) in zip(row, col, sec)]
	meta = hcat(meta, [path path_new_fn row col])
	writedlm(meta_new_fn, meta)

	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	function auto_adjust(img)
		distr = nquantile(img[:], 16)
		minval = distr[2]; maxval = distr[16];
		return min(1, max(0, (img-minval) / (maxval-minval)))
	end

	if construct
		construct_adjustment_image(meta)
	end
	fn = joinpath(RAW_DIR, "adjustment_image.txt")
	img_adj = Array{Float64, 2}(readdlm(fn))

	for i in 1:size(meta, 1)
		if 9 < meta[i,4]+1 < 13
			old_filepath = joinpath(RAW_DIR, meta[i,5])
			new_filepath = joinpath(RAW_DIR, meta[i,6])
			index = (1, meta[i,4]+1, meta[i,7], meta[i,8])
			offset = [meta[i,3], meta[i,2]]

			println(old_filepath)
			img = auto_adjust(read(old_filepath) ./ img_adj)
			img = Array{UInt8, 2}(round(UInt8, img * 255))
			#f = h5open(joinpath(to, import_frame_tif_to_wafer_section_h5(tif, wafer_num, sec_num)), "w")
			f = h5open(new_filepath, "w")
			#    	chunksize = div(min(size(img)...), 4);
			chunksize = 1000;
			@time f["img", "chunk", (chunksize,chunksize)] = img
			close(f)
			update_offset(index, offset, [size(img)...])
		end		
	end
end

function build_AIBS_actual_adjustment(meta, dir, suffix="_sec00000", n=20)
	println("building adjustment file")

	adj_fn = joinpath(RAW_DIR, string("adj", suffix, ".txt"))
	intensities = Array{Float64,1}()

	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	samples = [1:30:size(meta,1)]
	for (k, i) = enumerate(samples)
		println("Compiling meta $i, # $k / ", length(samples))
		old_fn = joinpath(dir, meta[i,1])
		if isfile(old_fn)
			intensities = vcat(intensities, read(old_fn)[:])
		end
	end
	println("Calculating intensity histogram...")
	distr = nquantile(intensities, n)
	println(distr)
	writedlm(adj_fn, distr)
end

function build_AIBS_actual_meta(meta_fn, suffix="_sec00000")
	meta = readdlm(meta_fn)
	meta_new_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))

	# cleanup meta text file
	path = meta[:,1]
	path_no_ext = [split(i, ".")[1] for i in path]
	col = [parse(split(i, "_")[9])+1 for i in path_no_ext]
	row = [parse(split(i, "_")[10])+1 for i in path_no_ext]
	sec = meta[:,4]+1
	# string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
	# path_new_fn = [string("Tile_r", @sprintf("%02d", i), "-c", @sprintf("%02d", j), "_S2-W001_sec", @sprintf("%04d", k), ".h5") for (i,j,k) in zip(row, col, sec)]
	path_new_fn = [string("Tile_r", i, "-c", j, "_S2-W001_sec", k, ".h5") for (i,j,k) in zip(row, col, sec)]
	meta = hcat(meta, [path_new_fn row col])
	writedlm(meta_new_fn, meta)
end

function import_AIBS_actual(dir="/media/tmacrina/Data/AIBS_actual_trial/Dropbox (MIT)/0", suffix="_sec00000", trakem_fn="_trackem_20160630142845_243774_7R_SID_09_redo_0_0_416.txt")
	meta_fn = joinpath(RAW_DIR, string("meta", suffix, ".txt"))
	if !isfile(meta_fn)
		build_AIBS_actual_meta(joinpath(dir, trakem_fn), suffix)
	end
	meta = readdlm(meta_fn)

	adj_fn = joinpath(RAW_DIR, string("adj", suffix, ".txt"))
	if !isfile(adj_fn)
		build_AIBS_actual_adjustment(meta, dir, suffix)
	end
	adj = readdlm(adj_fn)
	minval = adj[2]
	maxval = adj[length(adj)-1]

	function read(fn)
		return Array{Float64, 2}(Images.load(fn).data')
	end

	function auto_adjust(img) #, n=20)
		# distr = nquantile(img[:], n)
		# minval = distr[2]; maxval = distr[n];
		return min(1, max(0, (img-minval) / (maxval-minval)))
	end

	function convert_float64_to_uint8(img)
		return Array{UInt8, 2}(round(UInt8, img * 255))
	end

	for i in 1:size(meta, 1)
		println("# $i / ", size(meta,1))
		old_fn = joinpath(dir, meta[i,1])
		new_fn = joinpath(RAW_DIR, meta[i,5])
		index = (1, meta[i,4]+1, meta[i,6], meta[i,7])
		offset = [meta[i,3], meta[i,2]]

		println(old_fn)
		if isfile(old_fn) # && !isfile(new_fn)
			img = auto_adjust(read(old_fn))
			img = convert_float64_to_uint8(img)
			#f = h5open(joinpath(to, import_frame_tif_to_wafer_section_h5(tif, wafer_num, sec_num)), "w")
			f = h5open(new_fn, "w")
			#    	chunksize = div(min(size(img)...), 4);
			chunksize = 1000;
			@time f["img", "chunk", (chunksize,chunksize)] = img
			close(f)
			update_offset(index, offset, [size(img)...])
		else
			println("Does not exist / already exists!")
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


