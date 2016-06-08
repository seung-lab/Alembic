#NORMALIZER = Array{Float64, 2}(readdlm("calib"))
using StatsBase

function import_AIBS_dir_tifs(src_folder, suffix="_prealigned", dst_folder=src_folder)
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