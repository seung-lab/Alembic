#NORMALIZER = Array{Float64, 2}(readdlm("calib"))
using StatsBase

function import_AIBS_dir_tifs(folder, wafer_num = 0, sec_num = 0, to = folder)
	for tif in sort_dir(folder, "tif")
	  println(tif)
	  img = Array{Float64, 2}(Images.load(joinpath(folder, tif)).data')
#	  img = img ./ NORMALIZER;
#	  distr = nquantile(img[:], 16)
#	  minval = distr[2]; maxval = distr[16];
#	  img = min(1, max(0, (img - minval) / (maxval-minval)))
	  img = Array{UInt8, 2}(round(UInt8, img * 255));
        #f = h5open(joinpath(to, import_frame_tif_to_wafer_section_h5(tif, wafer_num, sec_num)), "w")
        f = h5open(joinpath(to, import_aligned_tif_to_wafer_section_h5(tif)), "w")
#    	chunksize = div(min(size(img)...), 4);
    	chunksize = 1000;
    	@time f["img", "chunk", (chunksize,chunksize)] = img
        close(f)
	update_offset(aligned(1, int(tif[1:3])), [0, 0], [size(img)...])
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

function import_aligned_tif_to_wafer_section_h5(aligned_name)
	return string("1", ",", int(aligned_name[1:3]), "_", "aligned", ".h5")
end
