#=
# "approximate" functions operate on a Mesh, and will try to come up with the best fit for the given nodes / nodes_t
# :solve functions operate on a Match (within a MeshSet), and will try to come up with the best fit for the given matches
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
  src_points, dst_points = get_correspondences(match; globalized=globalized, filtered=true);
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
  return calculate_rigid(pts_src, pts_dst)
end

#="""
Return the right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""=#
function affine_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_nodes(M)
  return calculate_affine(pts_src, pts_dst)
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
function affine_approximate(ms::MeshSet, row, col)
	ind = findfirst(i -> ms.meshes[i].index[3:4] == (row, col), 1:ms.N)
 	return affine_approximate(ms.meshes[ind])
end


#="""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""=#
function affine_solve(ms::MeshSet; k=1, globalized=false)
  pts_src, pts_dst = get_homogeneous_correspondences(ms, k; globalized=globalized)
	for ind in 1:size(pts_dst, 1)
  	pts_dst[ind, 1:2] = pts_dst[ind, 1:2] - get_offset(ms.matches[k].src_index)'
	end
  return calculate_affine(pts_src, pts_dst)
end

function affine_solve!(ms::MeshSet; k=1, globalized=false)
  tform = affine_solve(ms; k=k, globalized=globalized)
  apply_transform!(ms, tform, k)
end

#="""
Return right-hand matrix for the matches
`pts_src*T ~= pts_dst`
"""=#
function rigid_solve(ms::MeshSet; k=1, globalized=false)
  pts_src, pts_dst = get_homogeneous_correspondences(ms, k; globalized=globalized)
	for ind in 1:size(pts_dst, 1)
  	pts_dst[ind, 1:2] = pts_dst[ind, 1:2] - get_offset(ms.matches[k].src_index)'
	end
  return calculate_rigid(pts_src, pts_dst)
end

function rigid_solve!(ms::MeshSet; k=1, globalized=false)
  tform = rigid_solve(ms; k=k, globalized=globalized)
  apply_transform!(ms, tform, k)
end

#="""
Apply some weighted combination of affine and rigid, gauged by lambda
"""=#
function regularized_solve(ms::MeshSet; k=1, lambda=ms.properties[:params][:solve][:lambda], globalized=false)
  affine = affine_solve(ms; k=k, globalized=globalized)
  rigid = rigid_solve(ms; k=k, globalized=globalized)
  return lambda*affine + (1-lambda)*rigid
end

#="""
ONLY WORKS ON CASES WHERE MATCHES[1] = MESHES[1] -> MESHES[2]
"""=#
function regularized_solve!(ms::MeshSet; k=1, lambda=0.9)
	tform = regularized_solve(ms, k=k, lambda=lambda);
	apply_transform!(ms, tform, k)
end

function apply_transform!(ms::MeshSet, tform, k::Int64)
  for ind in 1:count_nodes(ms.meshes[k])
    ms.meshes[k].dst_nodes[ind] = ([ms.meshes[k].src_nodes[ind]; 1]' * tform)[1:2]
  end
end

function translation_solve(ms::MeshSet; globalized=false)
  indices = map(get_index, ms.meshes)
  indices = spiral_sort(indices)
  # incomplete
end

function solve!(meshset; method=meshset.properties[:params][:solve][:method])
	sanitize!(meshset);
  assert(count_matches(meshset) != 0)
  assert(count_filtered_correspondences(meshset) >= 3)

	if method == :elastic elastic_solve!(meshset); end
	if method == :translation translation_solve!(meshset); end
	if method == :rigid rigid_solve!(meshset); end
	if method == :regularized regularized_solve!(meshset; lambda = meshset.properties[:params][:solve][:lambda]); end
	if method == :affine affine_solve!(meshset); end
  mark_solved!(meshset)
end

function elastic_solve_piecewise!(meshset::MeshSet; from_current = true)
	sanitize!(meshset);
	meshsets = make_submeshsets(meshset);
	for (index, cur_meshset) in enumerate(meshsets)
	  println("SOLVING SUBMESHSET $index OF $(length(meshsets))");
	  elastic_solve!(cur_meshset; from_current = from_current);
	end
	return meshset;
end

function elastic_collate(meshset; from_current=true, manual_only=false, batch_size=1)
	sanitize!(meshset);
  match_spring_coeff = PARAMS[:solve][:match_spring_coeff]
  mesh_spring_coeff = PARAMS[:solve][:mesh_spring_coeff]
  max_iters = PARAMS[:solve][:max_iters]
  ftol_cg = PARAMS[:solve][:ftol_cg]
  eta_gd = PARAMS[:solve][:eta_gd] 
  ftol_gd = PARAMS[:solve][:ftol_gd]
  eta_newton = PARAMS[:solve][:eta_newton]
  ftol_newton = PARAMS[:solve][:ftol_newton]

  println("Solving meshset: $(count_nodes(meshset)) nodes, $(count_edges(meshset)) edges, $(count_filtered_correspondences(meshset, manual_only=manual_only)) correspondences");

  nodes = Array{Float64, 2}(2, count_nodes(meshset));
  nodes_fixed = BinaryProperty(falses(count_nodes(meshset)));

  edge_lengths = FloatProperty(count_edges(meshset) + count_filtered_correspondences(meshset, manual_only=manual_only))
  edge_spring_coeffs = FloatProperty(count_edges(meshset) + count_filtered_correspondences(meshset, manual_only=manual_only))
  #edges = spzeros(count_nodes(meshset), 0)

  noderanges = Dict{Any, Any}();
  edgeranges = Dict{Any, Any}();
  meshes = Dict{Any, Any}();
  meshes_order = Dict{Any, Int64}();
  cum_nodes = 0; cum_edges = 0;

  meshes_ref = Array{RemoteChannel, 1}()
  matches_ref = Array{RemoteChannel, 1}()
  src_indices = Array{Any, 1}();
  dst_indices = Array{Any, 1}();

  @fastmath @inbounds begin
    @fastmath @inbounds for (k, mesh) in enumerate(meshset.meshes)
      noderanges[get_index(mesh)] = cum_nodes + (1:count_nodes(mesh))
      edgeranges[get_index(mesh)] = cum_edges + (1:count_edges(mesh))
      meshes[get_index(mesh)] = mesh
      meshes_order[get_index(mesh)] = k;
      cum_nodes = cum_nodes + count_nodes(mesh);
      cum_edges = cum_edges + count_edges(mesh);
      mesh_ref = RemoteChannel(); 
      put!(mesh_ref, mesh); push!(meshes_ref, mesh_ref);
    end

    println("meshes referenced")

    @fastmath @inbounds for match in meshset.matches
    	edgeranges[get_src_and_dst_indices(match)] = cum_edges + (1:count_filtered_correspondences(match, manual_only=manual_only));
    	cum_edges = cum_edges + count_filtered_correspondences(match, manual_only=manual_only);
    	match_ref = RemoteChannel(); 
    	put!(match_ref, match); push!(matches_ref, match_ref);
    	push!(src_indices, get_src_index(match))
    	push!(dst_indices, get_dst_index(match))
    end

    println("matches referenced")

    for mesh in meshset.meshes
      if from_current
        @fastmath @inbounds nodes[:, noderanges[get_index(mesh)]] = get_nodes(mesh; globalized = true, use_post = true)[:];
      else
        @fastmath @inbounds nodes[:, noderanges[get_index(mesh)]] = get_nodes(mesh; globalized = true, use_post = false)[:];
      end
      if is_fixed(mesh)
        @fastmath @inbounds nodes_fixed[noderanges[get_index(mesh)]] = fill(true, count_nodes(mesh));
      end
      #@inbounds edge_lengths[edgeranges[get_index(mesh)]] = get_edge_lengths(mesh);
      #@inbounds edge_spring_coeffs[edgeranges[get_index(mesh)]] = fill(mesh_spring_coeff, count_edges(mesh));

      edge_lengths[edgeranges[get_index(mesh)]] = get_edge_lengths(mesh);
      edge_spring_coeffs[edgeranges[get_index(mesh)]] = mesh_spring_coeff
      # removed_edges = get_removed_edge_indices(mesh)
      # edge_spring_coeffs[edgeranges[get_index(mesh)][removed_edges]] = 0
      # fixed_edges = get_fixed_edge_indices(mesh)
      # edge_spring_coeffs[edgeranges[get_index(mesh)][fixed_edges]] = 10000
    end

  end # @fm @ib 

  println("mesh nodes collated: $(count_meshes(meshset)) meshes.")

  @fastmath noderange_mesh_list = Array{UnitRange, 1}([getindex(noderanges, get_index(mesh)) for mesh in meshset.meshes]);
  @fastmath edgerange_mesh_list = Array{UnitRange, 1}([getindex(edgeranges, get_index(mesh)) for mesh in meshset.meshes]);


  function compute_sparse_matrix_meshes(mesh_ref, noderange, edgerange)
	@inbounds begin
	  mesh = fetch(mesh_ref)
	  i_inds = copy(mesh.edges_I)
	  j_inds = copy(mesh.edges_J)
	  vals = copy(mesh.edges_V)

	noderange_start = noderange[1] - 1
	edgerange_start = edgerange[1] - 1
				
	for i in 1:length(i_inds)
		@fastmath @inbounds i_inds[i] += noderange_start
	end
	for i in 1:length(j_inds)
		@fastmath @inbounds j_inds[i] += edgerange_start
	end
	end #ib
	return i_inds, j_inds, vals
  end

  res_meshes = pmap(compute_sparse_matrix_meshes, meshes_ref, noderange_mesh_list, edgerange_mesh_list, batch_size=batch_size)

  i_inds_mesh = vcat([r[1] for r in res_meshes]...)
  j_inds_mesh = vcat([r[2] for r in res_meshes]...)
  vals_mesh = vcat([r[3] for r in res_meshes]...)

  println("mesh edges collated: $(count_meshes(meshset)) meshes.")

  # initialize the arrays for edge lengths and the spring coeffs
  for match in meshset.matches
    	@inbounds edge_lengths[edgeranges[get_src_and_dst_indices(match)]] = fill(0, count_filtered_correspondences(match, manual_only=manual_only));
      z_diff = abs(get(Z_MAP, match.dst_index, match.dst_index) - get(Z_MAP, match.src_index, match.src_index))
    	@inbounds edge_spring_coeffs[edgeranges[get_src_and_dst_indices(match)]] = fill(match_spring_coeff / z_diff, count_filtered_correspondences(match, manual_only=manual_only));
  end

  function compute_sparse_matrix_matches(match_ref, src_mesh_ref, dst_mesh_ref, noderange_src, noderange_dst, edgerange)
	@inbounds begin
	@time begin
  	match = fetch(match_ref)
  	print("match $(get_src_index(match))->$(get_dst_index(match)) being collated... fetch: ")
  	src_mesh = fetch(src_mesh_ref)
  	dst_mesh = fetch(dst_mesh_ref)
      	end;
	
	print("triangulation: "); @time begin
	src_pts, dst_pts = get_correspondences(match; filtered = true, manual_only=manual_only);
	src_pt_triangles = find_mesh_triangle(src_mesh, src_pts);
	dst_pt_triangles = find_mesh_triangle(dst_mesh, dst_pts);
	src_pt_weights = get_triangle_weights(src_mesh, src_pts, src_pt_triangles);
	dst_pt_weights = get_triangle_weights(dst_mesh, dst_pts, dst_pt_triangles);
      end  


	print("allocation: "); @time begin
	i_inds_src = Int64[]
	j_inds_src = Int64[]
	vals_src = Float64[]
	i_inds_dst = Int64[]
	j_inds_dst = Int64[]
	vals_dst = Float64[]
	for ind in 1:count_filtered_correspondences(match, manual_only=manual_only)
	  src_t1, src_t2, src_t3 = src_pt_triangles[ind];
	  dst_t1, dst_t2, dst_t3 = dst_pt_triangles[ind];
	  if src_t1 == 0 || dst_t1 == 0 continue; end
#	  if src_tri == NO_TRIANGLE || dst_tri == NO_TRIANGLE continue; end
	  src_pt_w1, src_pt_w2, src_pt_w3 = src_pt_weights[ind]
	  dst_pt_w1, dst_pt_w2, dst_pt_w3 = dst_pt_weights[ind]
	  if src_pt_w1 > eps
	    	push!(i_inds_src, src_t1);
	    	push!(j_inds_src, ind);
	    	push!(vals_src, -src_pt_w1);
	  end
	  if src_pt_w2 > eps
	    	push!(i_inds_src, src_t2);
	    	push!(j_inds_src, ind);
	    	push!(vals_src, -src_pt_w2);
	  end
	  if src_pt_w3 > eps
	    	push!(i_inds_src, src_t3);
	    	push!(j_inds_src, ind);
	    	push!(vals_src, -src_pt_w3);
	  end
	  if dst_pt_w1 > eps
	    	push!(i_inds_dst, dst_t1);
	    	push!(j_inds_dst, ind);
	    	push!(vals_dst, dst_pt_w1);
	  end
	  if dst_pt_w2 > eps
	    	push!(i_inds_dst, dst_t2);
	    	push!(j_inds_dst, ind);
	    	push!(vals_dst, dst_pt_w2);
	  end
	  if dst_pt_w3 > eps
	    	push!(i_inds_dst, dst_t3);
	    	push!(j_inds_dst, ind);
	    	push!(vals_dst, dst_pt_w3);
	  end
	end

	noderange_src_start = noderange_src[1] - 1
	noderange_dst_start = noderange_dst[1] - 1
	edgerange_start = edgerange[1] - 1
				
	for i in 1:length(i_inds_src)
		@fastmath @inbounds i_inds_src[i] += noderange_src_start
	end
	for i in 1:length(i_inds_dst)
		@fastmath @inbounds i_inds_dst[i] += noderange_dst_start
	end	
	for i in 1:length(j_inds_src)
		@fastmath @inbounds j_inds_src[i] += edgerange_start
	end
	for i in 1:length(j_inds_dst)
		@fastmath @inbounds j_inds_dst[i] += edgerange_start
	end
				
	i_inds = vcat(i_inds_src, i_inds_dst)
	j_inds = vcat(j_inds_src, j_inds_dst)
	vals = vcat(vals_src, vals_dst)
	
  end 
      end #inbounds

	return i_inds, j_inds, vals
  end


  noderange_src_list = Array{UnitRange, 1}(map(getindex, repeated(noderanges), src_indices))
  noderange_dst_list = Array{UnitRange, 1}(map(getindex, repeated(noderanges), dst_indices))
  edgerange_list = Array{UnitRange, 1}(map(getindex, repeated(edgeranges), map(get_src_and_dst_indices,meshset.matches)))

  res = pmap(compute_sparse_matrix_matches, matches_ref, 
  	meshes_ref[map(getindex, repeated(meshes_order), src_indices)], meshes_ref[map(getindex, repeated(meshes_order), dst_indices)], 
	noderange_src_list, noderange_dst_list, edgerange_list, batch_size=batch_size);

	i_inds_match = vcat([r[1] for r in res]...)
	j_inds_match = vcat([r[2] for r in res]...)
	vals_match = vcat([r[3] for r in res]...)

  println("matches collated: $(count_matches(meshset)) matches.")

	i_inds = vcat(i_inds_mesh, i_inds_match)
	j_inds = vcat(j_inds_mesh, j_inds_match)
	vals = vcat(vals_mesh, vals_match)


	print("populating sparse matrix:"); 
  @time edges = sparse(i_inds, j_inds, vals, count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset, manual_only=manual_only))

  collation = nodes, nodes_fixed, edges, edge_spring_coeffs, edge_lengths, max_iters, ftol_cg
  collation_with_ranges = collation, noderanges, edgeranges

  return collation_with_ranges;
end
"""
Elastic solve
"""
function elastic_solve!(ms; from_current=true, manual_only=false, batch_size=1)
	sanitize!(ms);
	collation, noderanges, edgeranges = elastic_collate(ms; from_current=from_current, manual_only=manual_only, batch_size=batch_size)
	BLAS.set_num_threads(Sys.CPU_CORES);
  @time SolveMeshConjugateGradient!(collation...)
  setblas();
  for m in ms.meshes
    dst_nodes = collation[1][:, noderanges[get_index(m)]]
  	broadcast!(-, m.dst_nodes, dst_nodes, get_offset(m));
  end
end

# invalids set to NO_POINT
function get_correspondences(meshset::MeshSet, ind::Int64; filtered=false, globalized::Bool=false, use_post=false, manual_only=false)
	if use_post
	  src_mesh = meshset.meshes[find_mesh_index(meshset,get_src_index(meshset.matches[ind]))]
	  dst_mesh = meshset.meshes[find_mesh_index(meshset,get_dst_index(meshset.matches[ind]))]
	  return get_correspondences(meshset.matches[ind]; filtered=filtered, globalized=globalized, use_post = use_post, src_mesh=src_mesh, dst_mesh=dst_mesh, manual_only=manual_only)
	end
	return get_correspondences(meshset.matches[ind]; filtered=filtered, globalized=globalized, use_post = use_post, manual_only=manual_only)
end

function get_displacements_post(meshset::MeshSet, ind)
  src_nodes, dst_nodes, filtered_inds = get_correspondences(meshset, ind; globalized = true, use_post = true)
  return src_nodes - dst_nodes, filtered_inds
end

function get_norms_post(meshset::MeshSet, ind)
  displacements, filtered_inds = get_displacements_post(meshset, ind)
  return map(norm, displacements), filtered_inds
end

function calculate_post_statistics!(meshset::MeshSet, match_ind)
  dv, _ = get_displacements_post(meshset, match_ind)
  norms, _ = get_norms_post(meshset, match_ind)
  for (cp, dv_i, norm_i) in zip(meshset.matches[match_ind].correspondence_properties, dv, norms)
    cp["posts"] = Dict("dv_post" => dv_i, "norm_post" => norm_i)
  end
end

function rectify_drift(meshset::MeshSet, start_ind = 1, final_ind = count_meshes(meshset), bias = [0.0, 0.0]; use_post = true, rectify_aligned = false)
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
	if use_post
	g_src_pts_after, g_dst_pts_after, filtered_after = get_correspondences(match; src_mesh = src_mesh, dst_mesh = dst_mesh, globalized = true, use_post = true)
  	residuals_match_post = g_dst_pts_after[filtered_after] - g_src_pts_after[filtered_after]; 
	avg_drift = mean(residuals_match_post);
      else
	g_src_pts, g_dst_pts, filtered = get_correspondences(match; globalized = true)
  	residuals_match = g_dst_pts[filtered] - g_src_pts[filtered]; 
	avg_drift = mean(residuals_match);
      end

      if PARAMS[:match][:symmetric]
	if is_preceding(get_index(src_mesh), get_index(dst_mesh), 5)
	  drifts[find_mesh_index(meshset, get_dst_index(match))] -= avg_drift/2;
	else
	  drifts[find_mesh_index(meshset, get_src_index(match))] += avg_drift/2;
	end
	else
	if is_preceding(get_index(src_mesh), get_index(dst_mesh), 5)
	  drifts[find_mesh_index(meshset, get_dst_index(match))] -= avg_drift;
	else
	  drifts[find_mesh_index(meshset, get_src_index(match))] += avg_drift;
	end
      end

  end
  cum_drift = Point([0,0]);
  for i in setdiff(1:length(drifts), start_ind:final_ind) drifts[i] = Point([0,0]); end
  for (ind, drift) in enumerate(drifts)
    cum_drift += (drift + bias);
    println("$(meshset.meshes[ind].index) is rectified against local drift $drift, cumulative drift $cum_drift...")
    for i in 1:count_nodes(meshset.meshes[ind])
	@fastmath @inbounds meshset.meshes[ind].dst_nodes[i] += cum_drift;
    end
    #update_offset(get_index(meshset.meshes[ind]), round(Int64, get_offset(get_index(meshset.meshes[ind])) + cum_drift))
    #=if rectify_aligned
    update_offset(aligned(get_index(meshset.meshes[ind])), round(Int64, get_offset(aligned(get_index(meshset.meshes[ind]))) + cum_drift))
  end=#
  end


end

function load_stats(ms::MeshSet)
  index = nextstage(get_index(ms.meshes[1]))
  path = get_path("stats", index)
  if !isfile(path)
    return calculate_stats(ms)
  else
    return load("stats", index)
  end
end  

function compile_stats(index_range)
  stats = Dict()
  for index in index_range
    path = get_path("stats", index)
    if !isfile(path)
      ms = load(MeshSet, index)
      stats[index] = calculate_stats(ms)
    else
      stats[index] = load("stats", index)
    end
  end
  return stats
end

"""
Filters is array of strings
"""
function filter_stats(stats, filters)
  f = []
  if haskey(stats, "matches")
    for (k, v) in stats["matches"]
      m = []
      for filter in filters
        push!(m, Base.get(v, filter, 0))
      end
      push!(f, [k, m...])
    end
  else
    for (k, v) in stats
      m = []
      for filter in filters
        push!(m, Base.get(v["summary"], filter, 0))
      end
      push!(f, [(k...), m...])
    end
  end
  f = hcat(f...)'
  return f[sortperm(f[:,1]),:]
end

function calculate_stats(index::FourTupleIndex)
  ms = load(MeshSet, index)
  return calculate_stats(ms)
end

@inbounds function calculate_stats(ms::MeshSet)

  stats = Dict()
  stats["matches"] = Dict()

  function make_stats_dict(src_index, dst_index, filtered_count,
                            rms_pre=0, avg_pre=0, std_pre=0, max_pre=0,
                            rms_post=0, avg_post=0, std_post=0, max_post=0,
                            drift_i=0, drift_j=0,
                            avg_r=0, std_r=0, min_r=0,
                            avg_95sig=0, std_95sig=0, max_95sig=0,
                            avg_75sig=0, std_75sig=0, max_75sig=0,
                            avg_50sig=0, std_50sig=0, max_50sig=0,
                            flagged=false)
    m = Dict()
    m["src_index"] = src_index
    m["dst_index"] = dst_index
    m["filtered_count"] = filtered_count

    m["rms_pre"] = rms_pre
    m["avg_pre"] = avg_pre
    m["std_pre"] = std_pre
    m["max_pre"] = max_pre
    m["rms_post"] = rms_post
    m["avg_post"] = avg_post
    m["std_post"] = std_post
    m["max_post"] = max_post
    m["drift_i"] = drift_i
    m["drift_j"] = drift_j
    m["avg_r"] = avg_r
    m["std_r"] = std_r
    m["min_r"] = min_r
    m["avg_95sig"] = avg_95sig
    m["std_95sig"] = std_95sig
    m["max_95sig"] = max_95sig
    m["avg_75sig"] = avg_75sig
    m["std_75sig"] = std_75sig
    m["max_75sig"] = max_75sig
    m["avg_50sig"] = avg_50sig
    m["std_50sig"] = std_50sig
    m["max_50sig"] = max_50sig
    m[:flagged] = flagged
    return m
  end

  println("Computing statistics...")

  params = get_params(ms)

  # these arrays have been reversed as necessary to get the right direction
  residuals_pre = Points(0)
  residuals_post = Points(0)
  avg_drifts = Points(0)
  total_corresps = 0;

  r_maxs = Array{Float64}(0)
  sigmas95 = Array{Float64}(0)
  sigmas75 = Array{Float64}(0)
  sigmas50 = Array{Float64}(0)
  matches_to_review = Array{Int64, 1}(0)

  meshes = Dict{Any, Any}();
  for mesh in ms.meshes
    meshes[get_index(mesh)] = mesh;
  end

  for match in ms.matches
    src_mesh = meshes[get_src_index(match)]
    dst_mesh = meshes[get_dst_index(match)]

    m = make_stats_dict(get_index(src_mesh), get_index(dst_mesh), count_filtered_correspondences(match))
    
    if count_filtered_correspondences(match) > 1

      g_src_pts, g_dst_pts, filtered = get_correspondences(match; globalized = true)
      g_src_pts_after, g_dst_pts_after, filtered_after = get_correspondences(match; src_mesh = src_mesh, dst_mesh = dst_mesh, globalized = true, use_post = true)

      props = get_filtered_correspondence_properties(match);
      r_maxs_match = Array{Float64}(map(get_dfs, props, repeated("r_max")));
      sigs95_match = Array{Float64}(map(get_dfs, props, repeated(Symbol("xcorr_sigma_0.95"))));
      sigs75_match = Array{Float64}(map(get_dfs, props, repeated(Symbol("xcorr_sigma_0.75"))));
      sigs50_match = Array{Float64}(map(get_dfs, props, repeated(Symbol("xcorr_sigma_0.5"))));

      g_src_pts = g_src_pts[filtered];
      g_dst_pts = g_dst_pts[filtered];

      g_src_pts_after = g_src_pts_after[filtered_after];
      g_dst_pts_after = g_dst_pts_after[filtered_after];

      @fastmath @inbounds begin

        residuals_match_pre = g_dst_pts - g_src_pts;
        residuals_match_post = g_dst_pts_after - g_src_pts_after;

        avg_drift = mean(residuals_match_post);
        res_norm = Array{Float64}(map(norm, residuals_match_pre))
        res_norm_post = Array{Float64}(map(norm, residuals_match_post))
        
        m = make_stats_dict(get_index(src_mesh), get_index(dst_mesh), count_filtered_correspondences(match),
                            sqrt(mean(res_norm.^2)), mean(res_norm), std(res_norm), maximum(res_norm),
                            sqrt(mean(res_norm_post.^2)), mean(res_norm_post), std(res_norm_post), maximum(res_norm_post),
                            avg_drift[1], avg_drift[2],
                            mean(r_maxs_match), std(r_maxs_match), minimum(r_maxs_match),
                            mean(sigs95_match), std(sigs95_match), maximum(sigs95_match),
                            mean(sigs75_match), std(sigs75_match), maximum(sigs75_match),
                            mean(sigs50_match), std(sigs50_match), maximum(sigs50_match),
                            is_flagged(match))


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
        append!(sigmas95, sigs95_match)
        append!(sigmas75, sigs75_match)
        append!(sigmas50, sigs50_match)
        total_corresps += count_filtered_correspondences(match);

      end #fmib
    end
    stats["matches"][find_match_index(ms, match)] = m
  end

  if total_corresps == 0
    m = make_stats_dict(get_index(ms.meshes[1]), get_index(ms.meshes[2]), total_corresps)
  else
    res_norm = Array{Float64}(map(norm, residuals_pre))
    res_norm_post = Array{Float64}(map(norm, residuals_post))
    avg_drifts = mean(avg_drifts)

    m = make_stats_dict(get_index(ms.meshes[1]), get_index(ms.meshes[2]), total_corresps,
                        sqrt(mean(res_norm.^2)), mean(res_norm), std(res_norm), maximum(res_norm),
                        sqrt(mean(res_norm_post.^2)), mean(res_norm_post), std(res_norm_post), maximum(res_norm_post),
                        avg_drifts[1], avg_drifts[2],
                        mean(r_maxs), std(r_maxs), minimum(r_maxs),
                        mean(sigmas95), std(sigmas95), maximum(sigmas95),
                        mean(sigmas75), std(sigmas75), maximum(sigmas75),
                        mean(sigmas50), std(sigmas50), maximum(sigmas50),
                        is_flagged(ms))
  end
  stats["summary"] = m

  path = get_path("stats", nextstage(get_index(ms.meshes[1])))
  println("Writing stats to $path")
  f = open(path, "w")
  write(f, JSON.json(stats))
  close(f)
  return stats
end

@inbounds function stats(meshset::MeshSet, first_ind = 1, last_ind = count_matches(meshset); flagged_only::Bool = false, summary::Bool = false)

  println("Computing statistics... sigma is computed at beta = 0.95")

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

  header = ""
  header = string(header, "index   ")
	header = string(header, "index   ")
	header = string(header, "src_index      ")
	header = string(header, "dst_index      ")
	header = string(header, "corrs")
	header = string(header, " | ")
	header = string(header, "rms_pre  ")
	header = string(header, "avg_pre  ")
	header = string(header, "std_pre  ")
	header = string(header, "max_pre")
	header = string(header, " | ")
	header = string(header, "rms_post ")
	header = string(header, "avg_post ")
	header = string(header, "std_post ")
	header = string(header, "max_post ")
	header = string(header, "drift_i ")
	header = string(header, "drift_j")
	header = string(header, " | ")
	header = string(header, "avg_r  ")
	header = string(header, "std_r  ")
	header = string(header, "min_r   ")
	header = string(header, "avg_σ  ")
	header = string(header, "std_σ  ")
	header = string(header, "max_σ ")
	header = string(header, "  |  ")
	header = string(header, "flags")
  header = string(header, "\n")
  print(header)

  divider = string(join(fill('-', 190)), "\n")
  print(divider)

  for match in meshset.matches[first_ind:last_ind]
    if flagged_only && !is_flagged(match) continue; end
  	src_mesh = meshes[get_src_index(match)]
  	dst_mesh = meshes[get_dst_index(match)]

  	# handle for empty match
  	if count_filtered_correspondences(match) < 2
      body = ""
  		body = string(body, @sprintf("%4i", findfirst(this -> meshset.matches[this] == match, 1:count_matches(meshset))));
  		body = string(body, @sprintf("%14s", string(get_index(src_mesh))))
  		body = string(body, "->")
  		body = string(body, @sprintf("%14s", string(get_index(dst_mesh))))
  		body = string(body, @sprintf("%6i", count_filtered_correspondences(match)))
      body = string(body, "\n")
  		print(body)
  		continue;
  	end

  	g_src_pts, g_dst_pts, filtered = get_correspondences(match; globalized = true)
  	g_src_pts_after, g_dst_pts_after, filtered_after = get_correspondences(match; src_mesh = src_mesh, dst_mesh = dst_mesh, globalized = true, use_post = true)

  	props = get_filtered_correspondence_properties(match);
  	r_maxs_match = Array{Float64}(map(get_dfs, props, repeated("r_max")));
  	sigs_match = Array{Float64}(map(get_dfs, props, repeated(Symbol("xcorr_sigma_0.95"))));

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

        body = ""
      	body = string(body, @sprintf("%4i", find_match_index(meshset, match)));
      	body = string(body, " ")
      	body = string(body, @sprintf("%13s", string(get_index(src_mesh))))
      	body = string(body, "->")
      	body = string(body, @sprintf("%13s", string(get_index(dst_mesh))))
      	body = string(body, @sprintf("%10i", count_filtered_correspondences(match)))
      	body = string(body, " ")
      	body = string(body, rms_pre_s)
      	body = string(body, avg_pre_s)
      	body = string(body, std_pre_s)
      	body = string(body, max_pre_s)
      	body = string(body, " ")
      	body = string(body, rms_post_s)
      	body = string(body, avg_post_s)
      	body = string(body, std_post_s)
      	body = string(body, max_post_s)
      	body = string(body, " ")
      	body = string(body, avg_drift_di_s)
      	body = string(body, avg_drift_dj_s)
      	body = string(body, "  ")
      	body = string(body, avg_r_s)
      	body = string(body, std_r_s)
      	body = string(body, min_r_s)
      	body = string(body, " ")
      	body = string(body, avg_sig_s)
      	body = string(body, std_sig_s)
      	body = string(body, max_sig_s)
      	body = string(body, "        ")

      	# FLAG PARAMETERS
      	if is_flagged(match)
      		body = string(body, "*")
      		push!(matches_to_review, find_match_index(meshset, match))
      	end
      	body = string(body, "\n")
        print(body)

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

  divider = string(join(fill('=', 190)), "\n")
  print(divider)

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

  footer = ""
	footer = string(footer, " ALL")
	footer = string(footer, " ")
	footer = string(footer, @sprintf("%13s", string(get_index(meshset.meshes[1]))))
	footer = string(footer, "->")
	footer = string(footer, @sprintf("%13s", string(get_index(last(meshset.meshes)))))
	footer = string(footer, @sprintf("%10i", total_corresps))
	footer = string(footer, " ")
	footer = string(footer, rms_pre_s)
	footer = string(footer, avg_pre_s)
	footer = string(footer, std_pre_s)
	footer = string(footer, max_pre_s)
	footer = string(footer, " ")
	footer = string(footer, rms_post_s)
	footer = string(footer, avg_post_s)
	footer = string(footer, std_post_s)
	footer = string(footer, max_post_s)
	footer = string(footer, " ")
	footer = string(footer, avg_drift_di_s)
	footer = string(footer, avg_drift_dj_s)
	footer = string(footer, "  ")
	footer = string(footer, avg_r_s)
	footer = string(footer, std_r_s)
	footer = string(footer, min_r_s)
	footer = string(footer, " ")
	footer = string(footer, avg_sig_s)
	footer = string(footer, std_sig_s)
	footer = string(footer, max_sig_s)
	footer = string(footer, "        ")

	# FLAG PARAMETERS
	if is_flagged(meshset)
		footer = string(footer, "*")
	end
	footer = string(footer, "\n")
  print(footer)

  divider = string(join(fill('-', 190)), "\n")
  print(divider)

  postscript = ""
  postscript = string(postscript, "Statistics on $(length(first_ind:last_ind)) matches from $first_ind -> $last_ind")
  if flagged_only 
    postscript = string(postscript, " --- only the flagged matches are included in the summary") 
  end
  postscript = string(postscript, "\n")
  postscript = string(postscript, "$(length(matches_to_review)) matches flagged for review: $matches_to_review")
  postscript = string(postscript, "\n")
  print(postscript)
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

"""
Solve script for loading, solving, & saving meshes
"""
function solve(z_start, z_stop; fix_list=[], reset_list=[], from_current=false, manual_only=true)
  if z_start > z_stop; z_start, z_stop = z_stop, z_start; end
  println("Solving from $z_start to $z_stop, fixing $fix_list, reseting $reset_list.")
  PARAMS[:mesh][:z_start] = z_start
  PARAMS[:mesh][:z_stop] = z_stop
  ms = compile_meshset()
  removed_ms = prune!(ms, manual_only=manual_only); 
  unfix!(ms)
  fix!(ms, fix_list)
  map(reset!, Iterators.repeated(ms), reset_list)
  elastic_solve!(ms, from_current=from_current, manual_only=manual_only, batch_size=2048);
  b = get_bbox(ms, globalized=true, use_post=true, scale=get_scale(:dst_image)/get_scale(:match_image)) 
  name = string((z_start, z_stop))
  write(joinpath(homedir(), "$(name).txt"), string(b)) 
  unfix!(ms)
  split_meshset(ms); 
end

"""
Solve a meshset in blocks

Args:
  * z_start: starting z index
  * z_stop: final z index (can be <= or > z_start)
  * z_increment: signed increment
  * z_overlap: unsigned amount of overlap between pieces
  * fixed_indices: 0-indexed range for sections in each piece to be fixed
  * reset_indices: 0-indexed range for sections in each piece to be reset

Loops through range in blocks of size z_increment. Solves block, fixing first
section of the block, and letting the tail float. Next block overlaps with 
prior block by z_overlap sections.
"""
function piecewise_solve(z_start, z_stop; z_increment=20, z_overlap=2, fixed_indices=1, reset_indices=4:19, manual_only=true)
  increment = z_increment > 0 ? 1 : -1
  z_overlap *= increment
  fixed_indices *= increment
  reset_indices *= increment
  assert(length(intersect(fixed_indices, reset_indices)) == 0)
  assert(length(intersect(0:increment:z_increment+z_overlap-increment, fixed_indices)) == length(fixed_indices))
  assert(length(intersect(0:increment:z_increment+z_overlap-increment, reset_indices)) == length(reset_indices))
  z_starts = z_start:z_increment:z_stop-increment
  z_stops = collect(z_start+z_increment+z_overlap-increment:z_increment:z_stop)
  if z_stops[end] != z_stop
    push!(z_stops, z_stop)
  end
  z_ranges = zip(z_starts, z_stops)
  for (k, (range_start, range_stop)) in enumerate(z_ranges)
    full_range = range_start:increment:range_stop
    reset_range = intersect(full_range, range_start + collect(reset_indices))
    fixed_range = intersect(full_range, range_start + collect(fixed_indices))
    # println((range_start, range_stop, fixed_range, reset_range))
    @time solve(range_start, range_stop, fix_list=fixed_range, reset_list=reset_range, from_current=true, manual_only=manual_only)
  end
end

"""
Filter a meshset in pieces (for large meshsets)

Hard-coded filter operations for basil SOA_v3.02
"""
function piecewise_filter(z_start, z_stop; z_increment=20, z_overlap=2)
  increment = z_increment > 0 ? 1 : -1
  z_overlap *= increment
  z_starts = z_start:z_increment:z_stop-increment
  z_stops = collect(z_start+z_increment+z_overlap-increment:z_increment:z_stop)
  if z_stops[end] != z_stop
    push!(z_stops, z_stop)
  end
  z_ranges = zip(z_starts, z_stops)
  for (k, (range_start, range_stop)) in enumerate(z_ranges)
    if range_start > range_stop; range_start, range_stop = range_stop, range_start; end
    println("Filtering $range_start to $range_stop")
    PARAMS[:mesh][:z_start] = range_start
    PARAMS[:mesh][:z_stop] = range_stop
    ms = compile_meshset()
    filter!(ms, (1, :get_correspondence_properties, :<=, 0.07, Symbol("xcorr_delta_15")))
    filter!(ms, (2, :get_correspondence_properties, :>=, 120, Symbol("vects_norm")))
    split_meshset(ms)
  end
end