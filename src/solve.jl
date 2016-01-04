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
function affine_solve_meshset!(Ms)
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



"""
Elastic solve
"""
function solve_meshset!(meshset)
  match_coeff = (meshset.properties["params"])["match_coeff"]
  mesh_coeff = (meshset.properties["params"])["mesh_coeff"]
  eta_gradient = (meshset.properties["params"])["eta_gradient"]
  ftol_gradient = (meshset.properties["params"])["ftol_gradient"]
  eta_newton = (meshset.properties["params"])["eta_newton"]
  ftol_newton = (meshset.properties["params"])["ftol_newton"]

  nodes = Array{Float64, 2}(2, 0);
  nodes_fixed = BinaryProperty(0)
  edges = spzeros(Float64, count_nodes(meshset), count_edges(meshset) + count_filtered_correspondences(meshset));
  edge_lengths = FloatProperty(0)
  edge_coeffs = FloatProperty(0)

  blockranges = Dict{Any, Any}();
  meshes = Dict{Any, Any}();
  cum_nodes = 0; cum_edges = 0;
  for mesh in meshset.meshes
	blockranges[mesh.index] = (cum_nodes + (1:count_nodes(mesh)), cum_edges + (1:count_edges(mesh)))
	meshes[mesh.index] = mesh
	cum_nodes = cum_nodes + count_nodes(mesh);
	cum_edges = cum_edges + count_edges(mesh);
  end

  for match in meshset.matches
	blockranges[match] = cum_edges + (1:count_filtered_correspondences(match));
	cum_edges = cum_edges + count_filtered_correspondences(match);
  end

  for mesh in meshset.meshes
    nodes = hcat(nodes, mesh.src_nodes...);
    num_nodes = count_nodes(mesh);
    append!(nodes_fixed, fill(false, num_nodes));
    append!(edge_lengths, get_edge_lengths(mesh));
    append!(edge_coeffs, fill(mesh_coeff, count_edges(mesh)));
    edges[blockranges[mesh.index]...] = mesh.edges;
  end

  for match in meshset.matches
    	src_mesh = meshes[match.src_index];
    	dst_mesh = meshes[match.dst_index];
	src_pts, dst_pts = get_filtered_correspondences(match);
	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	src_pt_weights = map(get_triangle_weights, repeated(src_mesh), src_pts, src_pt_triangles);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	dst_pt_weights = map(get_triangle_weights, repeated(dst_mesh), dst_pts, dst_pt_triangles);

	for ind in 1:count_filtered_correspondences(match)
		edges[((blockranges[match.src_index])[1])[collect(src_pt_triangles[ind])], ind] = (-1) * collect(src_pt_weights[ind]);
		edges[((blockranges[match.dst_index])[1])[collect(dst_pt_triangles[ind])], ind] = collect(dst_pt_weights[ind]);
	end
    	append!(edge_lengths, fill(0, count_filtered_correspondences(match)));
    	append!(edge_coeffs, fill(match_coeff, count_filtered_correspondences(match)));
  end

	println("nodes, ", size(nodes));
	println("nodes_fixed, ", size(nodes_fixed));
	println("edges, ", size(edges));
	println("edge_coeffs, ", size(edge_coeffs));
	println("edge_lengths, ", size(edge_lengths));

  SolveMesh!(nodes, nodes_fixed, edges, edge_coeffs, edge_lengths, eta_gradient, ftol_gradient, eta_newton, ftol_newton)
  dst_nodes = Points(0)
  for i in 1:size(nodes, 2)
          push!(dst_nodes, vec(nodes[:, i]))
        end

  for mesh in meshset.meshes
	mesh.dst_nodes = dst_nodes[blockranges[mesh.index][1]];
  end
  println("DONE!");
end

function stats(Ms::MeshSet)
  residuals = Points(0)
  residuals_t = Points(0)
  movement_src = Points(0)
  movement_dst = Points(0)
  for k in 1:Ms.M
    for i in 1:Ms.matches[k].n
      w = Ms.matches[k].dst_weights[i]
      t = Ms.matches[k].dst_triangles[i]
      p = Ms.matches[k].src_points_indices[i]
      src = Ms.meshes[find_index(Ms, Ms.matches[k].src_index)]
      dst = Ms.meshes[find_index(Ms, Ms.matches[k].dst_index)]
      p1 = src.nodes[p]
      p2 = dst.nodes[t[1]] * w[1] + dst.nodes[t[2]] * w[2] + dst.nodes[t[3]] * w[3]
      p1_t = src.nodes_t[p]
      p2_t = dst.nodes_t[t[1]] * w[1] + dst.nodes_t[t[2]] * w[2] + dst.nodes_t[t[3]] * w[3]
      push!(residuals, p2-p1)
      push!(residuals_t, p2_t-p1_t)
      push!(movement_src, p1_t-p1)
      push!(movement_dst, p2_t-p2)
    end
  end


   res_norm = map(norm, residuals)
   rms = sqrt(mean(res_norm.^2))
   avg = mean(res_norm)
   sig = std(res_norm)
   max = maximum(res_norm)


   res_norm_t = map(norm, residuals_t)
   rms_t = sqrt(mean(res_norm_t.^2))
   avg_t = mean(res_norm_t)
   sig_t = std(res_norm_t)
   max_t = maximum(res_norm_t)

   println("Statistics: ###")
   println("Total number of matched points: $(Ms.m_e) across $(Ms.M) Match between $(Ms.N) Meshes\n")
   println("Residuals before solving: rms: $rms,  mean: $avg, sigma = $sig, max = $max\n")
   println("Residuals after solving: rms: $rms_t,  mean: $avg_t, sigma = $sig_t, max = $max_t\n")
   println("###")
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
