#=
# "approximate" functions operate on a Mesh, and will try to come up with the best fit for the given nodes / nodes_t
# "solve" functions operate on a Match (within a MeshSet), and will try to come up with the best fit for the given matches
=#

"""
Retrieve nodes and nodes_t of mesh, make homogenous, and transpose to Nx3
"""
function get_homogeneous_nodes(mesh::Mesh)
  return homogenize_points(mesh.src_nodes), homogenize_points(mesh.dst_nodes);
end

function get_homogeneous_nodes(meshset::MeshSet, k)
  return get_homogeneous_nodes(meshset.meshes[k]);
end

#filtered ones only
function get_homogeneous_correspondences(match::Match)
  src_points, dst_points = get_filtered_correspondences(match);
  return homogenize_points(src_points), homogenize_points(dst_points);
end

function get_homogeneous_correspondences(meshset::MeshSet, k)
  return get_homogeneous_correspondences(meshset.matches[k])
end

"""
Return right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""
function rigid_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_nodes(M)
  return find_rigid(pts_src, pts_dst)
end

"""
Return the right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""
function affine_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_nodes(M)
  return find_affine(pts_src, pts_dst)
end

"""
Apply some weighted combination of affine and rigid, gauged by lambda
"""
function regularized_approximate(M::Mesh, lambda=0.9)
  affine = affine_approximate(M)
  rigid = rigid_approximate(M)
  return lambda*affine + (1-lambda)*rigid
end

"""
Use in montage?
"""
function affine_approximate(Ms::MeshSet, row, col)
	ind = findfirst(i -> Ms.meshes[i].index[3:4] == (row, col), 1:Ms.N)
 	return affine_approximate(Ms.meshes[ind])
end


"""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""
function affine_solve(Ms::MeshSet, k=1)
  pts_src, pts_dst = get_homogeneous_correspondences(Ms, k)
  return find_affine(pts_src, pts_dst)
end

"""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""
function rigid_solve(Ms::MeshSet, k=1)
  pts_src, pts_dst = get_homogeneous_correspondences(Ms, k)
  return find_rigid(pts_src, pts_dst)
end

"""
Apply some weighted combination of affine and rigid, gauged by lambda
"""
function regularized_solve(Ms::MeshSet; k=1, lambda=0.9)
  affine = affine_solve(Ms, k)
  rigid = rigid_solve(Ms, k)
  return lambda*affine + (1-lambda)*rigid
end

"""
Calculate affine transform of the matches, and set nodes_t of the moving mesh
"""
function affine_solve!(Ms)
	tform = affine_solve(Ms, 1);	
	nodes = Ms.meshes[Ms.matches_pairs[1][1]].nodes
	nodes_t = Points(size(nodes, 1))

  for i in 1:size(nodes, 1)
    h_pt = [nodes[i]; 1]
    nodes_t[i] = (h_pt' * tform)[1:2]
  end
  Ms.meshes[Ms.matches_pairs[1][1]].nodes_t = nodes_t
  stats(Ms)
  return Ms
  offset = mesh.disp
end

function solve!(meshset; method="elastic")
	if method == "elastic" return elastic_solve!(meshset); end
	if method == "translate" return translate_solve!(meshset); end
	if method == "rigid" return rigid_solve!(meshset); end
	if method == "regularized" return regularized_solve!(meshset); end
	if method == "affine" return affine_solve!(meshset); end
end

"""
Elastic solve
"""
function elastic_solve!(meshset)
  params = get_params(meshset)
  match_spring_coeff = params["solve"]["match_spring_coeff"]
  mesh_spring_coeff = params["solve"]["mesh_spring_coeff"]
  ftol_cg = params["solve"]["ftol_cg"]

  println("Solving meshset: $(count_nodes(meshset)) nodes, $(count_edges(meshset)) edges, $(count_filtered_correspondences(meshset)) correspondences");

  nodes = Array{Float64, 2}(2, count_nodes(meshset));
  nodes_fixed = BinaryProperty(count_nodes(meshset))
  edges = spzeros(Float64, count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset));
  edge_lengths = FloatProperty(count_edges(meshset) + count_filtered_correspondences(meshset))
  edge_spring_coeffs = FloatProperty(count_edges(meshset) + count_filtered_correspondences(meshset))

  noderanges = Dict{Any, Any}();
  edgeranges = Dict{Any, Any}();
  meshes = Dict{Any, Any}();
  cum_nodes = 0; cum_edges = 0;

  for mesh in meshset.meshes
  	noderanges[mesh.index] = cum_nodes + (1:count_nodes(mesh))
	edgeranges[mesh.index] = cum_edges + (1:count_edges(mesh))
	meshes[mesh.index] = mesh
	cum_nodes = cum_nodes + count_nodes(mesh);
	cum_edges = cum_edges + count_edges(mesh);
  end

  for match in meshset.matches
	edgeranges[match] = cum_edges + (1:count_filtered_correspondences(match));
	cum_edges = cum_edges + count_filtered_correspondences(match);
  end

  for mesh in meshset.meshes
    nodes[:, noderanges[mesh.index]] = get_globalized_nodes_h(mesh)[1];
    nodes_fixed[noderanges[mesh.index]] = fill(false, count_nodes(mesh));
    edge_lengths[edgeranges[mesh.index]] = get_edge_lengths(mesh);
    edge_spring_coeffs[edgeranges[mesh.index]] = fill(mesh_spring_coeff, count_edges(mesh));
  end

  println("meshes collated: $(count_meshes(meshset)) meshes")

  for mesh in meshset.meshes
    edges[noderanges[mesh.index], edgeranges[mesh.index]] = mesh.edges;
  end

  for match in meshset.matches
    	src_mesh = meshes[match.src_index];
    	dst_mesh = meshes[match.dst_index];
	src_pts, dst_pts = get_filtered_correspondences(match);
	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	src_pt_weights = map(get_triangle_weights, repeated(src_mesh), src_pts, src_pt_triangles);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	dst_pt_weights = map(get_triangle_weights, repeated(dst_mesh), dst_pts, dst_pt_triangles);

	noderange_src = noderanges[match.src_index];
	noderange_dst = noderanges[match.dst_index];
	edgerange = edgeranges[match];

	for ind in 1:count_filtered_correspondences(match)
		edges[noderange_src[src_pt_triangles[ind][1]], edgerange[ind]] = -src_pt_weights[ind][1];
		edges[noderange_src[src_pt_triangles[ind][2]], edgerange[ind]] = -src_pt_weights[ind][2];
		edges[noderange_src[src_pt_triangles[ind][3]], edgerange[ind]] = -src_pt_weights[ind][3];
		edges[noderange_dst[dst_pt_triangles[ind][1]], edgerange[ind]] = dst_pt_weights[ind][1];
		edges[noderange_dst[dst_pt_triangles[ind][2]], edgerange[ind]] = dst_pt_weights[ind][2];
		edges[noderange_dst[dst_pt_triangles[ind][3]], edgerange[ind]] = dst_pt_weights[ind][3];
	end

    	edge_lengths[edgeranges[match]] = fill(0, count_filtered_correspondences(match));
    	edge_spring_coeffs[edgeranges[match]] = fill(match_spring_coeff, count_filtered_correspondences(match));
  end

  println("matches collated: $(count_matches(meshset)) matches")

  SolveMesh2!(nodes, nodes_fixed, edges, edge_spring_coeffs, edge_lengths, ftol_cg)
  dst_nodes = Points(0)
  for i in 1:size(nodes, 2)
          push!(dst_nodes, vec(nodes[:, i]))
        end

  for mesh in meshset.meshes
	mesh.dst_nodes = dst_nodes[noderanges[mesh.index]] - fill(get_offset(mesh), count_nodes(mesh));
  end

  stats(meshset);

end

# RENAME to get_filtered_correspondences
function get_matched_points_t(meshset, ind)
  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
	meshes[mesh.index] = mesh;
  end

  match = meshset.matches[ind];
  
	src_pts, dst_pts = get_filtered_correspondences(match);

	src_mesh = meshes[match.src_index]
	dst_mesh = meshes[match.dst_index]

	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	src_pt_weights = map(get_triangle_weights, repeated(src_mesh), src_pts, src_pt_triangles);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	dst_pt_weights = map(get_triangle_weights, repeated(dst_mesh), dst_pts, dst_pt_triangles);

	src_pts_after = map(get_tripoint_dst, repeated(src_mesh), src_pt_triangles, src_pt_weights);
	dst_pts_after = map(get_tripoint_dst, repeated(dst_mesh), dst_pt_triangles, dst_pt_weights);

	g_src_pts = src_pts + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts = dst_pts + fill(get_offset(match.dst_index), length(dst_pts));

	g_src_pts_after = src_pts_after + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts_after = dst_pts_after + fill(get_offset(match.dst_index), length(dst_pts));

	return g_src_pts_after, g_dst_pts_after;
end

function stats(meshset::MeshSet)

  println("Computing statistics...")

  residuals_pre = Points(0)
  residuals_post = Points(0)
  r_vals = Array{Float64}(0)
  matches_to_review = Array{Match, 1}(0)

  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
	meshes[mesh.index] = mesh;
  end

  	print("index    ")
  	print("src_index       ")
  	print("dst_index   ")
	print("corrs  ")
	print("    ")
	print("rms_pre   ")
	print("avg_pre   ")
	print("sig_pre   ")
	print("max_pre   ")
	print("    ")
	print("rms_post  ")
	print("avg_post  ")
	print("sig_post  ")
	print("max_post  ")
	print("     ")
	print("rms_r     ")
	print("avg_r     ")
	print("sig_r     ")
	print("min_r     ")
	print("flags")
	println();

  for match in meshset.matches
	src_pts, dst_pts = get_filtered_correspondences(match);
	props = get_filtered_correspondence_properties(match);

	src_mesh = meshes[match.src_index]
	dst_mesh = meshes[match.dst_index]

	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	src_pt_weights = map(get_triangle_weights, repeated(src_mesh), src_pts, src_pt_triangles);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	dst_pt_weights = map(get_triangle_weights, repeated(dst_mesh), dst_pts, dst_pt_triangles);

	src_pts_after = Points(map(get_tripoint_dst, repeated(src_mesh), src_pt_triangles, src_pt_weights));
	dst_pts_after = Points(map(get_tripoint_dst, repeated(dst_mesh), dst_pt_triangles, dst_pt_weights));

	g_src_pts = src_pts + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts = dst_pts + fill(get_offset(match.dst_index), length(dst_pts));

	g_src_pts_after = src_pts_after + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts_after = dst_pts_after + fill(get_offset(match.dst_index), length(dst_pts));

	residuals_match_pre = g_dst_pts - g_src_pts;
	residuals_match_post = g_dst_pts_after - g_src_pts_after;
	r_vals_match = Array{Float64}(map(getindex, props, repeated("r_val")));

   	res_norm = map(norm, residuals_match_pre)
   	rms_pre = sqrt(mean(res_norm.^2))
   	avg_pre = mean(res_norm)
   	sig_pre = std(res_norm)
   	max_pre = maximum(res_norm)

   	res_norm_post = map(norm, residuals_match_post)
   	rms_post = sqrt(mean(res_norm_post.^2))
   	avg_post = mean(res_norm_post)
   	sig_post = std(res_norm_post)
   	max_post = maximum(res_norm_post)

   	rms_r = sqrt(mean(r_vals_match.^2))
   	avg_r = mean(r_vals_match)
   	sig_r = std(r_vals_match)
   	min_r = minimum(r_vals_match)

   	rms_pre_s = @sprintf("%10.2f", rms_pre)
   	avg_pre_s = @sprintf("%10.2f", avg_pre) 
   	sig_pre_s = @sprintf("%10.2f", sig_pre)
   	max_pre_s = @sprintf("%10.2f", max_pre)

   	rms_post_s = @sprintf("%10.2f", rms_post)
   	avg_post_s = @sprintf("%10.2f", avg_post) 
   	sig_post_s = @sprintf("%10.2f", sig_post)
   	max_post_s = @sprintf("%10.2f", max_post)

   	rms_r_s = @sprintf("%10.3f", rms_r) 
   	avg_r_s = @sprintf("%10.3f", avg_r)
   	sig_r_s = @sprintf("%10.3f", sig_r)
   	min_r_s = @sprintf("%10.3f", min_r)


	print(@sprintf("%4i", findfirst(this -> meshset.matches[this] == match, 1:count_matches(meshset))));
	print(@sprintf("%14s", string(src_mesh.index)))
	print("->")
	print(@sprintf("%14s", string(dst_mesh.index)))
	print(@sprintf("%6i", count_filtered_correspondences(match)))
	print("    ")
	print(rms_pre_s)
	print(avg_pre_s)
	print(sig_pre_s)
	print(max_pre_s)
	print("    ")
	print(rms_post_s)
	print(avg_post_s)
	print(sig_post_s)
	print(max_post_s)
	print("    ")
	print(rms_r_s)
	print(avg_r_s)
	print(sig_r_s)
	print(min_r_s)
	print("       ")

	# FLAG PARAMETERS
	if rms_post > 10 || avg_post > 5 || max_post > 25 || min_r < .5
		print("*")
		push!(matches_to_review, match)
	end
	println()
	      
	append!(residuals_pre, residuals_match_pre)
	append!(residuals_post, residuals_match_post)
	append!(r_vals, r_vals_match)
  end


   	res_norm = map(norm, residuals_pre)
   	rms_pre = sqrt(mean(res_norm.^2))
   	avg_pre = mean(res_norm)
   	sig_pre = std(res_norm)
   	max_pre = maximum(res_norm)

   	res_norm_post = map(norm, residuals_post)
   	rms_post = sqrt(mean(res_norm_post.^2))
   	avg_post = mean(res_norm_post)
   	sig_post = std(res_norm_post)
   	max_post = maximum(res_norm_post)

   	rms_r = sqrt(mean(r_vals.^2))
   	avg_r = mean(r_vals)
   	sig_r = std(r_vals)
   	min_r = minimum(r_vals)

   	rms_pre_s = @sprintf("%10.2f", rms_pre)
   	avg_pre_s = @sprintf("%10.2f", avg_pre) 
   	sig_pre_s = @sprintf("%10.2f", sig_pre)
   	max_pre_s = @sprintf("%10.2f", max_pre)

   	rms_post_s = @sprintf("%10.2f", rms_post)
   	avg_post_s = @sprintf("%10.2f", avg_post) 
   	sig_post_s = @sprintf("%10.2f", sig_post)
   	max_post_s = @sprintf("%10.2f", max_post)

   	rms_r_s = @sprintf("%10.3f", rms_r) 
   	avg_r_s = @sprintf("%10.3f", avg_r)
   	sig_r_s = @sprintf("%10.3f", sig_r)
   	min_r_s = @sprintf("%10.3f", min_r)

   println("==============")
   println("Statistics across all matches")
   println("Residuals before solving: rms: $rms_pre_s,  mean: $avg_pre_s,  sigma = $sig_pre_s,  max = $max_pre_s")
   println("Residuals after solving:  rms: $rms_post_s,  mean: $avg_post_s,  sigma = $sig_post_s,  max = $max_post_s")
   println("r-values:                 rms: $rms_r_s,  mean: $avg_r_s,  sigma = $sig_r_s,  min = $min_r_s")

   println("==============")
   println("$(length(matches_to_review)) matches flagged for review")
   println() 
end

function decomp_affine(tform::Array{Float64, 2})
  # http://math.stackexchange.com/questions/78137/decomposition-of-a-nonsquare-affine-matrix
  a = tform[1, 1];
  b = tform[2, 1];
  c = tform[1, 2];
  d = tform[2, 2];

  p = norm([a, b])
  r = det(tform[1:2, 1:2]) / p
  q = (a * c + b * d) / det(tform)
  theta = rad2deg(atan(b / a))

  t_i = tform[3, 1]
  t_j = tform[3, 2]

  println("Affine decomposition: in right-to-left order with the transformation matrix being v*T,")
  println("Translation: i-translation: $t_i, j-translation: $t_j")
  println("Scaling:  i-scaling: $p, j-scaling: $r")
  println("Shearing: j-shear: $q")
  println("Rotation: $theta deg.")

  return t_i, t_j, p, r, q, theta
end
