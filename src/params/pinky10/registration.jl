addprocs();
using Alembic

load_params(joinpath(homedir(), ".julia/v0.6/Alembic/src/params/pinky10/stitched_vol19-vol34_registration.json"))

for src_z in 65:164
	trgt_z = src_z
	mesh = Mesh(src_z);
	src_image = get_image(src_z, :src_image);
	trgt_image = get_image(trgt_z, :match_image);

	match = Match(mesh, mesh, src_z, src_z, src_image, trgt_image);
	println((src_z, count_filtered_correspondences(match), count_rejected_correspondences(match)))

	pre, post = get_correspondences(match, globalized=false, filtered=false)
	src_pt_triangles = find_mesh_triangle(mesh, pre);
	src_pt_weights = get_triangle_weights(mesh, pre, src_pt_triangles);
	matched_inds = [round(Int, maximum(collect(a) .* collect(b))) for (a,b) in zip(src_pt_triangles, src_pt_weights)];
	mesh.dst_nodes[:,matched_inds] = match.dst_points;

	ms = MeshSet();
	add_mesh!(ms, mesh);
	save(ms, :mesh, string(src_z))
	save(match, :match)
	((warped_image, offset), _) = meshwarp_mesh(src_image, mesh);
	@time save_image(:dst_image, warped_image, offset, trgt_z, mip=get_mip(:dst_image), cdn_cache=false);
	reset_cache();
end
