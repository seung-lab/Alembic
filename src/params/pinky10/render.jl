addprocs();
using Alembic

# NOTE: Need to change global const IMG_ELTYPE to UInt16/UInt32 for segmentation files
# layers = ["seg_original"] # UInt16
layers = ["seg_invagination_new", "mit_new", "psd"] #UInt32
for l in layers
	load_params(joinpath(homedir(), ".julia/v0.6/Alembic/src/params/pinky10/stitched_vol19-vol34_$l.json"))

	for src_z in 65:164
		PARAMS[:mesh][:z_start] = src_z
		PARAMS[:mesh][:z_stop] = src_z
		ms = compile_meshset();
		mesh = ms.meshes[1]
		src_image = get_image(src_z, :src_image);
		((warped_image, offset), _) = meshwarp_mesh(src_image, mesh);
		@time save_image(:dst_image, warped_image, offset, src_z, mip=get_mip(:dst_image), cdn_cache=false);
		reset_cache();
	end
end
