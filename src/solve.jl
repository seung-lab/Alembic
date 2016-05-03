#=
# "approximate" functions operate on a Mesh, and will try to come up with the best fit for the given nodes / nodes_t
# "solve" functions operate on a Match (within a MeshSet), and will try to come up with the best fit for the given matches
=#

#="""
Retrieve nodes and nodes_t of mesh, make homogenous, and transpose to Nx3
"""=#
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

#="""
Return right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""=#
function rigid_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_nodes(M)
  return find_rigid(pts_src, pts_dst)
end

#="""
Return the right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""=#
function affine_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_nodes(M)
  return find_affine(pts_src, pts_dst)
end

#="""
Apply some weighted combination of affine and rigid, gauged by lambda
"""=#
function regularized_approximate(M::Mesh, lambda=0.9)
  affine = affine_approximate(M)
  rigid = rigid_approximate(M)
  return lambda*affine + (1-lambda)*rigid
end

#="""
Use in montage?
"""=#
function affine_approximate(Ms::MeshSet, row, col)
	ind = findfirst(i -> Ms.meshes[i].index[3:4] == (row, col), 1:Ms.N)
 	return affine_approximate(Ms.meshes[ind])
end


#="""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""=#
function affine_solve(Ms::MeshSet; k=1, globalized=false)
  pts_src, pts_dst = get_homogeneous_correspondences(Ms, k; globalized=globalized)
	for ind in 1:size(pts_dst, 1)
  	pts_dst[ind, 1:2] = pts_dst[ind, 1:2] - get_offset(Ms.matches[k].src_index)'
	end
  return find_affine(pts_src, pts_dst)
end

#="""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""=#
function rigid_solve(Ms::MeshSet; k=1, globalized=false)
  pts_src, pts_dst = get_homogeneous_correspondences(Ms, k; globalized=globalized)
	for ind in 1:size(pts_dst, 1)
  	pts_dst[ind, 1:2] = pts_dst[ind, 1:2] - get_offset(Ms.matches[k].src_index)'
	end
  return find_rigid(pts_src, pts_dst)
end

#="""
Apply some weighted combination of affine and rigid, gauged by lambda
"""=#
function regularized_solve(Ms::MeshSet; k=1, lambda=0.9, globalized=false)
  affine = affine_solve(Ms; k=k, globalized=globalized)
  rigid = rigid_solve(Ms; k=k, globalized=globalized)
  return lambda*affine + (1-lambda)*rigid
end

#="""
ONLY WORKS ON CASES WHERE MATCHES[1] = MESHES[1] -> MESHES[2]
"""=#
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

#="""
Elastic solve
"""=#
function elastic_solve!(meshset; from_current =false)
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
  #edges = spzeros(count_nodes(meshset), 0)

  noderanges = Dict{Any, Any}();
  edgeranges = Dict{Any, Any}();
  meshes = Dict{Any, Any}();
  meshes_order = Dict{Any, Int64}();
  cum_nodes = 0; cum_edges = 0;

  meshes_ref = Array{RemoteRef, 1}()
  matches_ref = Array{RemoteRef, 1}()
  src_indices = Array{Any, 1}();
  dst_indices = Array{Any, 1}();

  @fastmath @inbounds begin

  for (index, mesh) in enumerate(meshset.meshes)
  	noderanges[get_index(mesh)] = cum_nodes + (1:count_nodes(mesh))
	edgeranges[get_index(mesh)] = cum_edges + (1:count_edges(mesh))
	meshes[get_index(mesh)] = mesh
	meshes_order[get_index(mesh)] = index;
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
    if from_current
    nodes[:, noderanges[get_index(mesh)]] = get_globalized_nodes_h(mesh)[2];
  else
    nodes[:, noderanges[get_index(mesh)]] = get_globalized_nodes_h(mesh)[1];
  end
    if is_fixed(mesh)
    nodes_fixed[noderanges[get_index(mesh)]] = fill(true, count_nodes(mesh));
   # else
    #nodes_fixed[noderanges[get_index(mesh)]] = fill(false, count_nodes(mesh));
    end
    edge_lengths[edgeranges[get_index(mesh)]] = get_edge_lengths(mesh);
    edge_spring_coeffs[edgeranges[get_index(mesh)]] = fill(mesh_spring_coeff, count_edges(mesh));
  end

  end # @fm @ib 

  noderange_list = Array{UnitRange, 1}([getindex(noderanges, get_index(mesh)) for mesh in meshset.meshes]);
  edgerange_list = Array{UnitRange, 1}([getindex(edgeranges, get_index(mesh)) for mesh in meshset.meshes]);

  @inbounds @fastmath function make_local_sparse(num_nodes, num_edges)
	global LOCAL_SPM = spzeros(num_nodes, num_edges)
  end

@sync begin
   @async for proc in procs() remotecall_fetch(proc, make_local_sparse, count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset)); end 
 end

  function copy_sparse_matrix(mesh_ref, noderange, edgerange)
    mesh = fetch(mesh_ref)
    @inbounds LOCAL_SPM[noderange, edgerange] = mesh.edges;
  end

  pmap(copy_sparse_matrix, meshes_ref, noderange_list, edgerange_list);

#  edges_subarrays_meshes = Array{SparseMatrixCSC{Float64, Int64}, 1}(pmap(pad_sparse_matrix, meshes_ref, repeated(count_nodes(meshset)), noderange_list))

#  @time @inbounds @fastmath edges_subarrays_meshes = [hcat(edges_subarrays_meshes..., spzeros(count_nodes(meshset), count_filtered_correspondences(meshset)))]


  println("meshes collated: $(count_meshes(meshset)) meshes")

  for match in meshset.matches
    	edge_lengths[edgeranges[match]] = fill(0, count_filtered_correspondences(match));
    	edge_spring_coeffs[edgeranges[match]] = fill(match_spring_coeff, count_filtered_correspondences(match));
  end

  function compute_sparse_matrix(match_ref, src_mesh_ref, dst_mesh_ref, noderange_src, noderange_dst, edgerange)
	@inbounds begin
  	match = fetch(match_ref)
  	println("match $(get_src_index(match))->$(get_dst_index(match)) being collated...")
  	src_mesh = fetch(src_mesh_ref)
  	dst_mesh = fetch(dst_mesh_ref)

	src_pts, dst_pts = get_filtered_correspondences(match);
	src_pt_triangles = find_mesh_triangle(src_mesh, src_pts);
	dst_pt_triangles = find_mesh_triangle(dst_mesh, dst_pts);
	src_pt_weights = get_triangle_weights(src_mesh, src_pts, src_pt_triangles);
	dst_pt_weights = get_triangle_weights(dst_mesh, dst_pts, dst_pt_triangles);

	for ind in 1:count_filtered_correspondences(match)
		if src_pt_triangles[ind] == NO_TRIANGLE || dst_pt_triangles[ind] == NO_TRIANGLE continue; end
	        for i in 1:3
		  	if src_pt_weights[ind][i] > eps
			LOCAL_SPM[noderange_src[src_pt_triangles[ind][i]], edgerange[ind]] = -src_pt_weights[ind][i]
		        end
		  	if dst_pt_weights[ind][i] > eps
			LOCAL_SPM[noderange_dst[dst_pt_triangles[ind][i]], edgerange[ind]] = dst_pt_weights[ind][i]
		      end
		end
	end
      end #inbounds
  end


  noderange_src_list = Array{UnitRange, 1}(map(getindex, repeated(noderanges), src_indices))
  noderange_dst_list = Array{UnitRange, 1}(map(getindex, repeated(noderanges), dst_indices))
  
  edgerange_list = Array{UnitRange, 1}(map(getindex, repeated(edgeranges), meshset.matches))


  pmap(compute_sparse_matrix, matches_ref, meshes_ref[map(getindex, repeated(meshes_order), src_indices)], meshes_ref[map(getindex, repeated(meshes_order), dst_indices)], noderange_src_list, noderange_dst_list, edgerange_list);

  println("matches collated: $(count_matches(meshset)) matches. populating sparse matrix....")

  function get_local_sparse()
	return LOCAL_SPM;
  end
  
  edges_subarrays = Array{SparseMatrixCSC{Float64, Int64}, 1}(length(procs()))

@sync begin
   @async for proc in procs() edges_subarrays[proc] = remotecall_fetch(proc, get_local_sparse); end 
 end

  function add_local_sparse(sp_a, sp_b)
    global LOCAL_SPM = 0;
    global LOCAL_SPM = 0;
    gc();
    @fastmath global LOCAL_SPM = sp_a + sp_b
	return LOCAL_SPM;
  end

  @time @inbounds @fastmath while length(edges_subarrays) != 1
    println(length(edges_subarrays));
    if isodd(length(edges_subarrays)) push!(edges_subarrays, spzeros(count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset))) end
    @time edges_subarrays = Array{SparseMatrixCSC{Float64, Int64}, 1}(pmap(add_local_sparse, edges_subarrays[1:div(length(edges_subarrays), 2)], edges_subarrays[div(length(edges_subarrays),2)+1:end]))
  end

  edges = edges_subarrays[1];

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
	mesh.dst_nodes = dst_nodes[noderanges[get_index(mesh)]] - fill(get_offset(mesh), count_nodes(mesh));
  end

 # stats(meshset; summary = true);

end

function get_globalized_correspondences(meshset::MeshSet, ind::Int64)
  match = meshset.matches[ind];
  src_mesh = meshset.meshes[find_mesh_index(meshset, get_src_index(match))];
  dst_mesh = meshset.meshes[find_mesh_index(meshset, get_dst_index(match))];
	return get_globalized_correspondences(match, meshset.properties["params"]["registry"]["global_offsets"])
end

# may be invalid as well
function get_globalized_correspondences(match::Match, global_offsets::Bool = true)
  
	src_pts, dst_pts = get_correspondences(match);
	filtered_inds = get_filtered_indices(match);

	if !global_offsets
	@fastmath g_src_pts = src_pts + fill(get_offset(get_src_index(match)), length(src_pts));
	g_dst_pts = dst_pts;
	else
	@fastmath g_src_pts = src_pts + fill(get_offset(get_src_index(match)), length(src_pts));
	@fastmath g_dst_pts = dst_pts + fill(get_offset(get_dst_index(match)), length(dst_pts));
	end

	return g_src_pts, g_dst_pts, filtered_inds;
end

function get_globalized_correspondences_post(match::Match, src_mesh::Mesh, dst_mesh::Mesh, global_offsets::Bool = true)
	src_pts, dst_pts = get_correspondences(match);
	filtered_inds = get_filtered_indices(match);

	src_pt_triangles = find_mesh_triangle(src_mesh, src_pts);
	dst_pt_triangles = find_mesh_triangle(dst_mesh, dst_pts);

	src_pt_weights = get_triangle_weights(src_mesh, src_pts, src_pt_triangles);
	dst_pt_weights = get_triangle_weights(dst_mesh, dst_pts, dst_pt_triangles);

	src_pts_after = get_tripoint_dst(src_mesh, src_pt_triangles, src_pt_weights);
	dst_pts_after = get_tripoint_dst(dst_mesh, dst_pt_triangles, dst_pt_weights);

	if !global_offsets
	@fastmath g_src_pts_after = src_pts_after + fill(get_offset(get_src_index(match)), length(src_pts));
	g_dst_pts_after = dst_pts_after;
	else
	@fastmath g_src_pts_after = src_pts_after + fill(get_offset(get_src_index(match)), length(src_pts));
	@fastmath g_dst_pts_after = dst_pts_after + fill(get_offset(get_dst_index(match)), length(dst_pts));
	end

	return g_src_pts_after, g_dst_pts_after, filtered_inds;
end

# invalids set to NO_POINT
function get_globalized_correspondences_post(meshset::MeshSet, ind)
  match = meshset.matches[ind];
  src_mesh = meshset.meshes[find_mesh_index(meshset, get_src_index(match))];
  dst_mesh = meshset.meshes[find_mesh_index(meshset, get_dst_index(match))];
	return get_globalized_correspondences_post(match, src_mesh, dst_mesh, meshset.properties["params"]["registry"]["global_offsets"])
end

function rectify_drift(meshset::MeshSet)
  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
	meshes[mesh.index] = mesh;
  end

  drifts = Points(fill(Point([0,0]), count_meshes(meshset)));

  for match in meshset.matches
  	src_mesh = meshes[get_src_index(match)]
  	dst_mesh = meshes[get_dst_index(match)]

  	# handle for empty match
  	if count_filtered_correspondences(match) == 0
  		continue;
  	end

	g_src_pts_after, g_dst_pts_after, filtered_after = get_globalized_correspondences_post(match, src_mesh, dst_mesh, meshset.properties["params"]["registry"]["global_offsets"])

  	residuals_match_post[filtered_after] = g_dst_pts_after - g_src_pts_after;

	avg_drift = mean(residuals_match_post);

	if is_preceding(get_index(src_mesh), get_index(dst_mesh), 5)
	  drifts[find_mesh_index(meshset, get_dst_index(match))] -= avg_drift/2;
	else
	  drifts[find_mesh_index(meshset, get_src_index(match))] += avg_drift/2;
	end

  end
  cum_drift = Point([0,0]);
  for (ind, drift) in enumerate(drifts)
    cum_drift += drift;
    update_offset(get_index(meshset.meshes[ind]), get_offset(get_index(meshset.meshes[ind])) + cum_drift)
  end


end

@inbounds function stats(meshset::MeshSet, first_ind = 1, last_ind = count_matches(meshset); flagged_only::Bool = false, summary::Bool = false)

  println("Computing statistics... sigma is computed at beta = 0.8")

  params = get_params(meshset)

  # these arrays have been reversed as necessary to get the right direction
  residuals_pre = Points(0)
  residuals_post = Points(0)
  avg_drifts = Points(0)
  total_corresps = 0;

  r_maxs = Array{Float64}(0)
  sigmas = Array{Float64}(0)
  matches_to_review = Array{Int64, 1}(0)

  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
	meshes[get_index(mesh)] = mesh;
  end

	print("index   ")
	print("src_index      ")
	print("dst_index      ")
	print("corrs")
	print(" | ")
	print("rms_pre  ")
	print("avg_pre  ")
	print("std_pre  ")
	print("max_pre")
	print(" | ")
	print("rms_post ")
	print("avg_post ")
	print("std_post ")
	print("max_post ")
	print("drift_i ")
	print("drift_j")
	print(" | ")
	print("avg_r  ")
	print("std_r  ")
	print("min_r   ")

	print("avg_σ  ")
	print("std_σ  ")
	print("max_σ ")
	print("  |  ")
	print("flags")
	println();
println(join(fill('-', 190)))

  for match in meshset.matches[first_ind:last_ind]
        if flagged_only && !is_flagged(match) continue; end
  	src_mesh = meshes[get_src_index(match)]
  	dst_mesh = meshes[get_dst_index(match)]

  	# handle for empty match
  	if count_filtered_correspondences(match) < 2
  		print(@sprintf("%4i", findfirst(this -> meshset.matches[this] == match, 1:count_matches(meshset))));
  		print(@sprintf("%14s", string(get_index(src_mesh))))
  		print("->")
  		print(@sprintf("%14s", string(get_index(dst_mesh))))
  		print(@sprintf("%6i", count_filtered_correspondences(match)))
  		println()
  		continue;
  	end

	g_src_pts, g_dst_pts, filtered = get_globalized_correspondences(match, params["registry"]["global_offsets"])
	g_src_pts_after, g_dst_pts_after, filtered_after = get_globalized_correspondences_post(match, src_mesh, dst_mesh, params["registry"]["global_offsets"])

  	props = get_filtered_correspondence_properties(match);
  	r_maxs_match = Array{Float64}(map(get_dfs, props, repeated("r_max")));
  	sigs_match = Array{Float64}(map(get_dfs, props, repeated(0.8)));

	g_src_pts = g_src_pts[filtered];
	g_dst_pts = g_dst_pts[filtered];

	g_src_pts_after = g_src_pts_after[filtered_after];
	g_dst_pts_after = g_dst_pts_after[filtered_after];

	@fastmath @inbounds begin

  	residuals_match_pre = g_dst_pts - g_src_pts;
  	residuals_match_post = g_dst_pts_after - g_src_pts_after;

	avg_drift = mean(residuals_match_post);

	if !summary
   	res_norm = Array{Float64}(map(norm, residuals_match_pre))
   	rms_pre = sqrt(mean(res_norm.^2))
   	avg_pre = mean(res_norm)
   	std_pre = std(res_norm)
   	max_pre = maximum(res_norm)

   	res_norm_post = Array{Float64}(map(norm, residuals_match_post))
   	rms_post = sqrt(mean(res_norm_post.^2))
   	avg_post = mean(res_norm_post)
   	std_post = std(res_norm_post)
   	max_post = maximum(res_norm_post)

   	avg_r = mean(r_maxs_match)
   	std_r = std(r_maxs_match)
   	min_r = minimum(r_maxs_match)

   	avg_sig = mean(sigs_match)
   	std_sig = std(sigs_match)
   	max_sig = maximum(sigs_match)

   	rms_pre_s = @sprintf("%9.2f", rms_pre)
   	avg_pre_s = @sprintf("%9.2f", avg_pre) 
   	std_pre_s = @sprintf("%9.2f", std_pre)
   	max_pre_s = @sprintf("%9.2f", max_pre)

   	rms_post_s = @sprintf("%9.2f", rms_post)
   	avg_post_s = @sprintf("%9.2f", avg_post) 
   	std_post_s = @sprintf("%9.2f", std_post)
   	max_post_s = @sprintf("%9.2f", max_post)

   	avg_drift_di_s = @sprintf("%8.2f", avg_drift[1])
   	avg_drift_dj_s = @sprintf("%8.2f", avg_drift[2])

   	avg_r_s = @sprintf("%7.3f", avg_r)
   	std_r_s = @sprintf("%7.3f", std_r)
   	min_r_s = @sprintf("%7.3f", min_r)

   	avg_sig_s = @sprintf("%7.1f", avg_sig)
   	std_sig_s = @sprintf("%7.1f", std_sig)
   	max_sig_s = @sprintf("%7.1f", max_sig)

  	print(@sprintf("%4i", find_match_index(meshset, match)));
	print(" ")
  	print(@sprintf("%13s", string(get_index(src_mesh))))
  	print("->")
  	print(@sprintf("%13s", string(get_index(dst_mesh))))
  	print(@sprintf("%10i", count_filtered_correspondences(match)))
  	#print(" ")
  	print(rms_pre_s)
  	print(avg_pre_s)
  	print(std_pre_s)
  	print(max_pre_s)
  	print(" ")
  	print(rms_post_s)
  	print(avg_post_s)
  	print(std_post_s)
  	print(max_post_s)
  	print(" ")
  	print(avg_drift_di_s)
  	print(avg_drift_dj_s)
  	print("  ")
  	print(avg_r_s)
  	print(std_r_s)
  	print(min_r_s)
  	print(" ")
  	print(avg_sig_s)
  	print(std_sig_s)
  	print(max_sig_s)
  	print("        ")

  	# FLAG PARAMETERS
  	if is_flagged(match)
  		print("*")
  		push!(matches_to_review, find_match_index(meshset, match))
  	end
  	println()

      end

	# turning residuals around for drift calculation
	if !is_preceding(get_index(src_mesh), get_index(dst_mesh), 5)
		residuals_pre = -1 * residuals_pre;
		residuals_post = -1 * residuals_post;
		avg_drift = -1 * avg_drift;
	end
  	      
  	append!(residuals_pre, residuals_match_pre)
  	append!(residuals_post, residuals_match_post)
  	push!(avg_drifts, avg_drift)
  	append!(r_maxs, r_maxs_match)
  	append!(sigmas, sigs_match)
	total_corresps += count_filtered_correspondences(match);

      end #fmib
  end
  #println("---SUMMARY ACROSS $(length(first_ind:last_ind)) MATCHES-------------------------------------------------------------------------------")
println(join(fill('=', 190)))

   	res_norm = Array{Float64}(map(norm, residuals_pre))
   	rms_pre = sqrt(mean(res_norm.^2))
   	avg_pre = mean(res_norm)
   	std_pre = std(res_norm)
   	max_pre = maximum(res_norm)

   	res_norm_post = Array{Float64}(map(norm, residuals_post))
   	rms_post = sqrt(mean(res_norm_post.^2))
   	avg_post = mean(res_norm_post)
   	std_post = std(res_norm_post)
   	max_post = maximum(res_norm_post)

  	avg_drifts = mean(avg_drifts)

   	avg_r = mean(r_maxs)
   	std_r = std(r_maxs)
   	min_r = minimum(r_maxs)

   	avg_sig = mean(sigmas)
   	std_sig = std(sigmas)
   	max_sig = maximum(sigmas)

   	rms_pre_s = @sprintf("%9.2f", rms_pre)
   	avg_pre_s = @sprintf("%9.2f", avg_pre) 
   	std_pre_s = @sprintf("%9.2f", std_pre)
   	max_pre_s = @sprintf("%9.2f", max_pre)

   	rms_post_s = @sprintf("%9.2f", rms_post)
   	avg_post_s = @sprintf("%9.2f", avg_post) 
   	std_post_s = @sprintf("%9.2f", std_post)
   	max_post_s = @sprintf("%9.2f", max_post)

   	avg_drift_di_s = @sprintf("%8.2f", avg_drifts[1])
   	avg_drift_dj_s = @sprintf("%8.2f", avg_drifts[2])

   	avg_r_s = @sprintf("%7.3f", avg_r)
   	std_r_s = @sprintf("%7.3f", std_r)
   	min_r_s = @sprintf("%7.3f", min_r)

   	avg_sig_s = @sprintf("%7.1f", avg_sig)
   	std_sig_s = @sprintf("%7.1f", std_sig)
   	max_sig_s = @sprintf("%7.1f", max_sig)

	print(" ALL")
	print(" ")
  	print(@sprintf("%13s", string(get_index(meshset.meshes[1]))))
  	print("->")
  	print(@sprintf("%13s", string(get_index(last(meshset.meshes)))))
  	print(@sprintf("%10i", total_corresps))
  	#print(" ")
  	print(rms_pre_s)
  	print(avg_pre_s)
  	print(std_pre_s)
  	print(max_pre_s)
  	print(" ")
  	print(rms_post_s)
  	print(avg_post_s)
  	print(std_post_s)
  	print(max_post_s)
  	print(" ")
  	print(avg_drift_di_s)
  	print(avg_drift_dj_s)
  	print("  ")
  	print(avg_r_s)
  	print(std_r_s)
  	print(min_r_s)
  	print(" ")
  	print(avg_sig_s)
  	print(std_sig_s)
  	print(max_sig_s)
  	print("        ")

  	# FLAG PARAMETERS
  	if is_flagged(meshset)
  		print("*")
  	end
  	println()

println(join(fill('-', 190)))

  print("Statistics on $(length(first_ind:last_ind)) matches from $first_ind -> $last_ind")
  if flagged_only print(" --- only the flagged matches are included in the summary") end
  println();
  println("$(length(matches_to_review)) matches flagged for review: $matches_to_review")
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
