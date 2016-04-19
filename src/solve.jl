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
function get_homogeneous_correspondences(match::Match; globalized = false)
  src_points, dst_points = get_filtered_correspondences(match; globalized=globalized);
  return homogenize_points(src_points), homogenize_points(dst_points);
end

function get_homogeneous_correspondences(meshset::MeshSet, k; globalized = false)
  return get_homogeneous_correspondences(meshset.matches[k]; globalized=globalized)
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
function affine_solve(Ms::MeshSet; k=1, globalized=false)
  pts_src, pts_dst = get_homogeneous_correspondences(Ms, k; globalized=globalized)
	for ind in 1:size(pts_dst, 1)
  	pts_dst[ind, 1:2] = pts_dst[ind, 1:2] - get_offset(Ms.matches[k].src_index)'
	end
  return find_affine(pts_src, pts_dst)
end

"""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""
function rigid_solve(Ms::MeshSet; k=1, globalized=false)
  pts_src, pts_dst = get_homogeneous_correspondences(Ms, k; globalized=globalized)
	for ind in 1:size(pts_dst, 1)
  	pts_dst[ind, 1:2] = pts_dst[ind, 1:2] - get_offset(Ms.matches[k].src_index)'
	end
  return find_rigid(pts_src, pts_dst)
end

"""
Apply some weighted combination of affine and rigid, gauged by lambda
"""
function regularized_solve(Ms::MeshSet; k=1, lambda=0.9, globalized=false)
  affine = affine_solve(Ms; k=k, globalized=globalized)
  rigid = rigid_solve(Ms; k=k, globalized=globalized)
  return lambda*affine + (1-lambda)*rigid
end

"""
ONLY WORKS ON CASES WHERE MATCHES[1] = MESHES[1] -> MESHES[2]
"""
function regularized_solve!(Ms::MeshSet; k=1, lambda=0.9)
	tform = regularized_solve(Ms, k=k, lambda=lambda);
	for ind in 1:count_nodes(Ms.meshes[k])
		Ms.meshes[k].dst_nodes[ind] = ([Ms.meshes[k].src_nodes[ind]; 1]' * tform)[1:2]
	end
	stats(Ms);
end

function solve!(meshset)
  method=meshset.properties["params"]["solve"]["method"]
  solve!(meshset; method=method)
  mark_solved(meshset)
end

function solve!(meshset; method="elastic")
	sanitize!(meshset);
  assert(count_matches(meshset) != 0)
  assert(count_filtered_correspondences(meshset) != 0)

	if method == "elastic" return elastic_solve!(meshset); end
	if method == "translate" return translate_solve!(meshset); end
	if method == "rigid" return rigid_solve!(meshset); end
	if method == "regularized" return regularized_solve!(meshset; lambda = meshset.properties["params"]["solve"]["lambda"]); end
	if method == "affine" return affine_solve!(meshset); end
end

"""
Elastic solve
"""
function elastic_solve!(meshset; from_current=false)
  params = get_params(meshset)
  #fixed = get_fixed(meshset)
  match_spring_coeff = params["solve"]["match_spring_coeff"]
  mesh_spring_coeff = params["solve"]["mesh_spring_coeff"]
  max_iters = params["solve"]["max_iters"]
  ftol_cg = params["solve"]["ftol_cg"]
  eta_gd = params["solve"]["eta_gd"] 
  ftol_gd = params["solve"]["ftol_gd"]
  eta_newton = params["solve"]["eta_newton"]
  ftol_newton = params["solve"]["ftol_newton"]

  println("Solving meshset: $(count_nodes(meshset)) nodes, $(count_edges(meshset)) edges, $(count_filtered_correspondences(meshset)) correspondences");

  nodes = Array{Float64, 2}(2, count_nodes(meshset));
  nodes_fixed = BinaryProperty(falses(count_nodes(meshset)));

  edge_lengths = FloatProperty(count_edges(meshset) + count_filtered_correspondences(meshset))
  edge_spring_coeffs = FloatProperty(count_edges(meshset) + count_filtered_correspondences(meshset))

  noderanges = Dict{Any, Any}();
  edgeranges = Dict{Any, Any}();
  meshes = Dict{Any, Any}();
  meshes_order = Dict{Any, Int64}();
  cum_nodes = 0; cum_edges = 0;

  meshes_ref = Array{RemoteRef, 1}()
  matches_ref = Array{RemoteRef, 1}()
  src_indices = Array{Any, 1}();
  dst_indices = Array{Any, 1}();

  for (index, mesh) in enumerate(meshset.meshes)
  	noderanges[mesh.index] = cum_nodes + (1:count_nodes(mesh))
	edgeranges[mesh.index] = cum_edges + (1:count_edges(mesh))
	meshes[mesh.index] = mesh
	meshes_order[mesh.index] = index;
	cum_nodes = cum_nodes + count_nodes(mesh);
	cum_edges = cum_edges + count_edges(mesh);
	mesh_ref = RemoteRef(); 
	put!(mesh_ref, mesh); push!(meshes_ref, mesh_ref);
  end

  for match in meshset.matches
	edgeranges[match] = cum_edges + (1:count_filtered_correspondences(match));
	cum_edges = cum_edges + count_filtered_correspondences(match);
	match_ref = RemoteRef(); 
	put!(match_ref, match); push!(matches_ref, match_ref);
	push!(src_indices, get_src_index(match))
	push!(dst_indices, get_dst_index(match))
  end

  for mesh in meshset.meshes
    nodes[:, noderanges[mesh.index]] = get_globalized_nodes_h(mesh)[2];
    if is_fixed(mesh)
    nodes_fixed[noderanges[mesh.index]] = fill(true, count_nodes(mesh));
   # else
    #nodes_fixed[noderanges[mesh.index]] = fill(false, count_nodes(mesh));
    end
    edge_lengths[edgeranges[mesh.index]] = get_edge_lengths(mesh);
    edge_spring_coeffs[edgeranges[mesh.index]] = fill(mesh_spring_coeff, count_edges(mesh));
  end

  println("meshes collated: $(count_meshes(meshset)) meshes")

  for match in meshset.matches
    	edge_lengths[edgeranges[match]] = fill(0, count_filtered_correspondences(match));
    	edge_spring_coeffs[edgeranges[match]] = fill(match_spring_coeff, count_filtered_correspondences(match));
  end

  function compute_sparse_entries(match_ref, src_mesh_ref, dst_mesh_ref, noderange_src, noderange_dst, edgerange)

  	match = fetch(match_ref)
  	println("match $(match.src_index)->$(match.dst_index) being collated...")
  	src_mesh = fetch(src_mesh_ref)
  	dst_mesh = fetch(dst_mesh_ref)

	src_pts, dst_pts = get_filtered_correspondences(match);
	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	src_pt_weights = map(get_triangle_weights, repeated(src_mesh), src_pts, src_pt_triangles);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	dst_pt_weights = map(get_triangle_weights, repeated(dst_mesh), dst_pts, dst_pt_triangles);


	edges_to_add = Array{Tuple{Int64, Int64, Float64}, 1}();
	for ind in 1:count_filtered_correspondences(match)
		if src_pt_triangles[ind] == NO_TRIANGLE || dst_pt_triangles[ind] == NO_TRIANGLE continue; end
	        for i in 1:3
			push!(edges_to_add, (noderange_src[src_pt_triangles[ind][i]], edgerange[ind], -src_pt_weights[ind][i]))
			push!(edges_to_add, (noderange_dst[dst_pt_triangles[ind][i]], edgerange[ind], dst_pt_weights[ind][i]))
		end
	end
	return edges_to_add;
  end

  edges_to_add = pmap(compute_sparse_entries, matches_ref, meshes_ref[map(getindex, repeated(meshes_order), src_indices)], meshes_ref[map(getindex, repeated(meshes_order), dst_indices)], map(getindex, repeated(noderanges), src_indices), map(getindex, repeated(noderanges), dst_indices), map(getindex, repeated(edgeranges), meshset.matches));
  
  edges_to_add = vcat(edges_to_add...)
  num_entries = length(edges_to_add);
  println("populating sparse matrix: $(num_entries) entries")

  @time begin
  node_inds = Array{Int64, 1}(num_entries);
  edge_inds = Array{Int64, 1}(num_entries);
  tri_weights = Array{Float64, 1}(num_entries);

  for (index, edge) in enumerate(edges_to_add)
    node_inds[index], edge_inds[index], tri_weights[index] = edge
  end
    edges = sparse(node_inds, edge_inds, tri_weights, count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset));
  #edges = spzeros(Float64, count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset));

  for mesh in meshset.meshes
    edges[noderanges[mesh.index], edgeranges[mesh.index]] = mesh.edges;
  end

  println("matches collated: $(count_matches(meshset)) matches")

	end #time
#return nodes, nodes_fixed, edges, edge_spring_coeffs, edge_lengths, max_iters, ftol_cg;

 #println(find(this -> this == true, nodes_fixed));
  if params["solve"]["use_conjugate_gradient"]
    @time SolveMeshConjugateGradient!(nodes, nodes_fixed, edges, edge_spring_coeffs, edge_lengths, max_iters, ftol_cg)
  else
    @time SolveMeshGDNewton!(nodes, nodes_fixed, edges, edge_spring_coeffs, edge_lengths, eta_gd, ftol_gd, eta_newton, ftol_newton)
  end
  dst_nodes = Points(0)
  for i in 1:size(nodes, 2)
          push!(dst_nodes, vec(nodes[:, i]))
        end

  for mesh in meshset.meshes
	mesh.dst_nodes = dst_nodes[noderanges[mesh.index]] - fill(get_offset(mesh), count_nodes(mesh));
  end

  stats(meshset);

end

# may be invalid as well
function get_globalized_correspondences(meshset, ind)
  	match = meshset.matches[ind];
  
	src_pts, dst_pts = get_correspondences(match);
	filtered_inds = get_filtered_indices(match);

	if is_montaged(meshset.meshes[1].index) || (haskey(meshset.properties["params"], "registry") && !meshset.properties["params"]["registry"]["global_offsets"])
	g_src_pts = src_pts + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts = dst_pts;
	else
	g_src_pts = src_pts + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts = dst_pts + fill(get_offset(match.dst_index), length(dst_pts));
	end

	return g_src_pts, g_dst_pts, filtered_inds;
end

# invalids set to NO_POINT
function get_globalized_correspondences_post(meshset, ind)
  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
	meshes[mesh.index] = mesh;
  end

  match = meshset.matches[ind];
  
	src_pts, dst_pts = get_correspondences(match);
	filtered_inds = get_filtered_indices(match);

	src_mesh = meshes[match.src_index]
	dst_mesh = meshes[match.dst_index]

	src_pt_triangles = find_mesh_triangle(src_mesh, src_pts);
	src_pt_weights = get_triangle_weights(src_mesh, src_pts, src_pt_triangles);
	dst_pt_triangles = find_mesh_triangle(dst_mesh, dst_pts);
	dst_pt_weights = get_triangle_weights(dst_mesh, dst_pts, dst_pt_triangles);

	src_pts_after = map(get_tripoint_dst, repeated(src_mesh), src_pt_triangles, src_pt_weights);
	dst_pts_after = map(get_tripoint_dst, repeated(dst_mesh), dst_pt_triangles, dst_pt_weights);

	if is_montaged(meshset.meshes[1].index) || (haskey(meshset.properties["params"], "registry") && !meshset.properties["params"]["registry"]["global_offsets"])
	g_src_pts_after = src_pts_after + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts_after = dst_pts_after;
	else
	g_src_pts_after = src_pts_after + fill(get_offset(match.src_index), length(src_pts));
	g_dst_pts_after = dst_pts_after + fill(get_offset(match.dst_index), length(dst_pts));
	end

	return g_src_pts_after, g_dst_pts_after, filtered_inds;
end

function stats(meshset::MeshSet)

  println("Computing statistics...")

  params = get_params(meshset)

  residuals_pre = Points(0)
  residuals_post = Points(0)
  r_maxs = Array{Float64}(0)
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
	print("std_pre   ")
	print("max_pre   ")
	print("    ")
	print("rms_post  ")
	print("avg_post  ")
	print("std_post  ")
	print("max_post  ")
	print("     ")
	print("rms_r     ")
	print("avg_r     ")
	print("std_r     ")
	print("min_r     ")
	print("flags")
	println();

  for match in meshset.matches
  	src_mesh = meshes[match.src_index]
  	dst_mesh = meshes[match.dst_index]

  	# handle for empty match
  	if count_filtered_correspondences(match) == 0
  		print(@sprintf("%4i", findfirst(this -> meshset.matches[this] == match, 1:count_matches(meshset))));
  		print(@sprintf("%14s", string(src_mesh.index)))
  		print("->")
  		print(@sprintf("%14s", string(dst_mesh.index)))
  		print(@sprintf("%6i", count_filtered_correspondences(match)))
  		println()
  		continue;
  	end

  	src_pts, dst_pts = get_filtered_correspondences(match);
  	props = get_filtered_correspondence_properties(match);

  	src_pt_triangles = find_mesh_triangle(src_mesh, src_pts);
  	src_pt_weights = get_triangle_weights(src_mesh, src_pts, src_pt_triangles);
  	dst_pt_triangles = find_mesh_triangle(dst_mesh, dst_pts);
  	dst_pt_weights = get_triangle_weights(dst_mesh, dst_pts, dst_pt_triangles);

  	src_pts_after = Points(map(get_tripoint_dst, repeated(src_mesh), src_pt_triangles, src_pt_weights));
  	dst_pts_after = Points(map(get_tripoint_dst, repeated(dst_mesh), dst_pt_triangles, dst_pt_weights));

  	g_src_pts = src_pts + fill(get_offset(match.src_index), length(src_pts));
  	g_src_pts_after = src_pts_after + fill(get_offset(match.src_index), length(src_pts));

  	if !haskey(params, "registry") || params["registry"]["global_offsets"]
  		g_dst_pts = dst_pts + fill(get_offset(match.dst_index), length(dst_pts));
  		g_dst_pts_after = dst_pts_after + fill(get_offset(match.dst_index), length(dst_pts));
  	else
  		g_dst_pts = dst_pts;
  		g_dst_pts_after = dst_pts_after;

  	end

  	residuals_match_pre = g_dst_pts - g_src_pts;
  	residuals_match_post = g_dst_pts_after - g_src_pts_after;
  	r_maxs_match = Array{Float64}(map(get_dfs, props, repeated("r_max")));

   	res_norm = map(norm, residuals_match_pre)
   	rms_pre = sqrt(mean(res_norm.^2))
   	avg_pre = mean(res_norm)
   	std_pre = std(res_norm)
   	max_pre = maximum(res_norm)

   	res_norm_post = map(norm, residuals_match_post)
   	rms_post = sqrt(mean(res_norm_post.^2))
   	avg_post = mean(res_norm_post)
   	std_post = std(res_norm_post)
   	max_post = maximum(res_norm_post)

   	rms_r = sqrt(mean(r_maxs_match.^2))
   	avg_r = mean(r_maxs_match)
   	std_r = std(r_maxs_match)
   	min_r = minimum(r_maxs_match)

   	rms_pre_s = @sprintf("%10.2f", rms_pre)
   	avg_pre_s = @sprintf("%10.2f", avg_pre) 
   	std_pre_s = @sprintf("%10.2f", std_pre)
   	max_pre_s = @sprintf("%10.2f", max_pre)

   	rms_post_s = @sprintf("%10.2f", rms_post)
   	avg_post_s = @sprintf("%10.2f", avg_post) 
   	std_post_s = @sprintf("%10.2f", std_post)
   	max_post_s = @sprintf("%10.2f", max_post)

   	rms_r_s = @sprintf("%10.3f", rms_r) 
   	avg_r_s = @sprintf("%10.3f", avg_r)
   	std_r_s = @sprintf("%10.3f", std_r)
   	min_r_s = @sprintf("%10.3f", min_r)

  	print(@sprintf("%4i", findfirst(this -> meshset.matches[this] == match, 1:count_matches(meshset))));
  	print(@sprintf("%14s", string(src_mesh.index)))
  	print("->")
  	print(@sprintf("%14s", string(dst_mesh.index)))
  	print(@sprintf("%6i", count_filtered_correspondences(match)))
  	print("    ")
  	print(rms_pre_s)
  	print(avg_pre_s)
  	print(std_pre_s)
  	print(max_pre_s)
  	print("    ")
  	print(rms_post_s)
  	print(avg_post_s)
  	print(std_post_s)
  	print(max_post_s)
  	print("    ")
  	print(rms_r_s)
  	print(avg_r_s)
  	print(std_r_s)
  	print(min_r_s)
  	print("       ")

  	# FLAG PARAMETERS
  	if is_flagged(match)
  		print("*")
  		push!(matches_to_review, match)
  	end
  	println()
  	      
  	append!(residuals_pre, residuals_match_pre)
  	append!(residuals_post, residuals_match_post)
  	append!(r_maxs, r_maxs_match)
  end

 	res_norm = map(norm, residuals_pre)
 	rms_pre = sqrt(mean(res_norm.^2))
 	avg_pre = mean(res_norm)
 	std_pre = std(res_norm)
 	max_pre = maximum(res_norm)

 	res_norm_post = map(norm, residuals_post)
 	rms_post = sqrt(mean(res_norm_post.^2))
 	avg_post = mean(res_norm_post)
 	std_post = std(res_norm_post)
 	max_post = maximum(res_norm_post)

 	rms_r = sqrt(mean(r_maxs.^2))
 	avg_r = mean(r_maxs)
 	std_r = std(r_maxs)
 	min_r = minimum(r_maxs)

 	rms_pre_s = @sprintf("%10.2f", rms_pre)
 	avg_pre_s = @sprintf("%10.2f", avg_pre) 
 	std_pre_s = @sprintf("%10.2f", std_pre)
 	max_pre_s = @sprintf("%10.2f", max_pre)

 	rms_post_s = @sprintf("%10.2f", rms_post)
 	avg_post_s = @sprintf("%10.2f", avg_post) 
 	std_post_s = @sprintf("%10.2f", std_post)
 	max_post_s = @sprintf("%10.2f", max_post)

 	rms_r_s = @sprintf("%10.3f", rms_r) 
 	avg_r_s = @sprintf("%10.3f", avg_r)
 	std_r_s = @sprintf("%10.3f", std_r)
 	min_r_s = @sprintf("%10.3f", min_r)

  println("==============")
  println("Statistics across all matches")
  println("Residuals before solving: rms: $rms_pre_s,  mean: $avg_pre_s,  std: $std_pre_s,  max: $max_pre_s")
  println("Residuals after solving:  rms: $rms_post_s,  mean: $avg_post_s,  std: $std_post_s,  max: $max_post_s")
  println("r-values:                 rms: $rms_r_s,  mean: $avg_r_s,  std: $std_r_s,  min: $min_r_s")

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
