#NORMALIZER = Array{Float64, 2}(readdlm("calib"))
using StatsBase

function import_AIBS_pilot_v1()
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

	for i in 1:size(meta, 1)
		old_filepath = joinpath(RAW_DIR, meta[i,5])
		new_filepath = joinpath(RAW_DIR, meta[i,6])
		index = (1, meta[i,4]+1, meta[i,7], meta[i,8])
		offset = [meta[i,3], meta[i,2]]

		println(old_filepath)
		img = Array{Float64, 2}(Images.load(old_filepath).data')
		#	  img = img ./ NORMALIZER;
		#	  distr = nquantile(img[:], 16)
		#	  minval = distr[2]; maxval = distr[16];
		#	  img = min(1, max(0, (img - minval) / (maxval-minval)))
		img = Array{UInt8, 2}(round(UInt8, img * 255));
		#f = h5open(joinpath(to, import_frame_tif_to_wafer_section_h5(tif, wafer_num, sec_num)), "w")
		f = h5open(new_filepath, "w")
		#    	chunksize = div(min(size(img)...), 4);
		chunksize = 1000;
		@time f["img", "chunk", (chunksize,chunksize)] = img
		close(f)
		update_offset(index, offset, [size(img)...])		
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


