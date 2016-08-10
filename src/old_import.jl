function import_trakem_dir_tifs(src_dir)
	for fn in sort_dir(src_dir, "tif")[1:2]
		println(fn)
		img = get_image_disk(joinpath(src_dir, fn), UInt8)
		# img = convert_float64_to_uint8(img)
		index = prealigned(1, fn[2:5])
		dst_fn = get_path(index)
		f = h5open(dst_fn, "w")
		chunksize = 1000;
		@time f["img", "chunk", (chunksize,chunksize)] = img
		close(f)
		update_offset(index, [0,0], [size(img)...])
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