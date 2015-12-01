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

function get_homogeneous_correspondences(match::Match)
  return homogenize_points(match.src_points), homogenize_points(match.dst_points);
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
end


function solve_meshset(meshset)

end

"""
Elastic solve
"""
function solve_meshset!(Ms)
  match_coeff = Ms.params["match_coeff"]
  eta_gradient = Ms.params["eta_gradient"]
  ftol_gradient = Ms.params["ftol_gradient"]
  eta_newton = Ms.params["eta_newton"]
  ftol_newton = Ms.params["ftol_newton"]

  nodes = Points(0)
  nodes_fixed = BinaryProperty(0)
  edges = spzeros(Float64, Ms.n, 0)
  edge_lengths = FloatProperty(0)
  edge_coeffs = FloatProperty(0)

  for i in 1:Ms.N
    cur_mesh = Ms.meshes[i]
    if i == 1 nodes = hcat(cur_mesh.nodes...)
    else nodes = hcat(nodes, hcat(cur_mesh.nodes...)); end
    append!(nodes_fixed, cur_mesh.nodes_fixed)
    append!(edge_lengths, cur_mesh.edge_lengths)
    append!(edge_coeffs, cur_mesh.edge_coeffs)
    if (i == Ms.N)  
      edges = hcat(edges, vcat(spzeros(Float64, Ms.nodes_indices[i], cur_mesh.m), 
             cur_mesh.edges))
    else
      edges = hcat(edges, vcat(spzeros(Float64, Ms.nodes_indices[i], cur_mesh.m), 
             cur_mesh.edges, 
             spzeros(Float64, Ms.n - Ms.nodes_indices[i] - cur_mesh.n, cur_mesh.m)))
    end
  end

  for i in 1:Ms.M
    cur_matches = Ms.matches[i]
    append!(edge_lengths, fill(0.0, cur_matches.n))
    append!(edge_coeffs, fill(convert(Float64, match_coeff), cur_matches.n))
    edges_padded = spzeros(Float64, Ms.n, cur_matches.n)

    for j in 1:Ms.matches[i].n
      edges_padded[find_node(Ms, Ms.matches_pairs[i][1], cur_matches.src_points_indices[j]), j] = -1
      edges_padded[find_node(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][1]), j] = cur_matches.dst_weights[j][1]
      edges_padded[find_node(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][2]), j] = cur_matches.dst_weights[j][2]
      edges_padded[find_node(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][3]), j] = cur_matches.dst_weights[j][3]
    end
    edges = hcat(edges, edges_padded)
  end

  SolveMesh!(nodes, nodes_fixed, edges, edge_coeffs, edge_lengths, eta_gradient, ftol_gradient, eta_newton, ftol_newton)
  nodes_t = Points(0)
  for i in 1:size(nodes, 2)
          push!(nodes_t, vec(nodes[:, i]))
        end
  for i in 1:Ms.N
    cur_mesh = Ms.meshes[i]
    cur_mesh.nodes_t = nodes_t[Ms.nodes_indices[i] + (1:cur_mesh.n)]
  end

    print(Ms.params)
    stats(Ms)
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
