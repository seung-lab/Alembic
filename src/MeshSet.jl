#import PyPlot

type MeshSet
  params::Dict

  N::Int64            # number of meshes in the set
  M::Int64            # number of matches in the set - (a -> b) and (b -> a) are distinct

  indices::Array{Index, 1}        # wafer, section, row, column as a tuple - if tileIndex is 0 then denotes entire section

  n::Int64            # number of nodes in the set across the whole set
  m::Int64            # number of edges in the set across the whole set
  m_i::Int64            # number of internal edges in the set across the whole set
  m_e::Int64            # number of edges between meshes
  
  meshes::Array{Mesh, 1}          # vector of meshes in the set
  nodes_indices::Array{Int64, 1}        # vector of number of nodes before the n-th mesh. To get index at i-th node of n-th mesh, nodes_indices[n] + i.

  matches::Array{Matches, 1}        # vector of matches in the set
  matches_pairs::Pairings       # vector of index (in meshes) - (a, b) means the match is between (meshes[a] -> meshes[b])
end


function find_node(Ms, mesh_ind, node_ind)
  return Ms.nodes_indices[mesh_ind] + node_ind
end

function find_index(Ms, mesh_index_tuple::Index)
  return findfirst(this -> mesh_index_tuple == this.index, Ms.meshes)
end

function find_section(Ms, section_num)
  return findfirst(this -> section_num == this.index[2], Ms.meshes)
end
function MeshSet(params::Dict)
  N = 0
  M = 0

  indices = Array{Index, 1}(0)
  
  n = 0
  m = 0
  m_i = 0
  m_e = 0

  meshes = Array{Mesh, 1}(0)
  nodes_indices = Array{Int64, 1}(0)

  matches = Array{Matches, 1}(0)
  matches_pairs = Array{Pair, 1}(0)

  return MeshSet(params, N, M, indices, n, m, m_i, m_e, meshes, nodes_indices, matches, matches_pairs)
end


function add_mesh(Am, Ms)
  push!(Ms.indices, Am.index)
  push!(Ms.meshes, Am)
  if length(Ms.nodes_indices) == 0 push!(Ms.nodes_indices, 0)
  else push!(Ms.nodes_indices, Ms.n)
  end
  Ms.N += 1
  Ms.m_i += Am.m
  Ms.m += Am.m
  Ms.n += Am.n
  return
end

function add_matches(M, Ms)
  if (typeof(M) == Void || M == Void) return; end
  push!(Ms.matches, M)
  push!(Ms.matches_pairs, (find_index(Ms, M.src_index), find_index(Ms, M.dst_index)))
  Ms.M += 1
  Ms.n
  Ms.m += M.n
  Ms.m_e += M.n
  return
end
#=
function consolidate(Ms::MeshSet)
  nodes_t = Points(0)
  
  for i in 1:Ms.N
    cur_mesh = Ms.meshes[i]
    if i == 1 nodes_t = hcat(cur_mesh.nodes_t...)
    else nodes_t = hcat(nodes_t, hcat(cur_mesh.nodes_t...)); end
  end
  
  nodes = Points(0)
  
  for i in 1:Ms.N
    cur_mesh = Ms.meshes[i]
    if i == 1 nodes = hcat(cur_mesh.nodes...)
    else nodes = hcat(nodes_t, hcat(cur_mesh.nodes...)); end
  end
end
=#

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

function save(filename::String, Ms::MeshSet)
  println("Saving meshset to ", filename)
  jldopen(filename, "w") do file
    write(file, "MeshSet", Ms)
  end
end

function save(Ms::MeshSet)
  firstindex = Ms.meshes[1].index
  lastindex = Ms.meshes[Ms.N].index


  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex))
    filename = joinpath(PREALIGNED_DIR, string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_prealigned.jld"))
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_prealigned(lastindex))
    filename = joinpath(ALIGNED_DIR, string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned.jld"))
  else 
    filename = joinpath(MONTAGED_DIR, string(join(firstindex[1:2], ","), "_montaged.jld"))
  end

  println("Saving meshset to ", filename)
  jldopen(filename, "w") do file
    write(file, "MeshSet", Ms)
  end
end

"""
Load montaged meshset for given wafer and section

`load_montaged(wafer_num, sec_num)`
"""
function load_montaged(wafer_num, sec_num)
  index = (wafer_num, sec_num, 1, 1)
	return load(index, index)
end

"""
Load prealigned meshset for given wafer and section

`load_prealigned(wafer_num, sec_num)`
"""
function load_prealigned(wafer_num, sec_num)
  lastindex = (wafer_num, sec_num, MONTAGED_INDEX, MONTAGED_INDEX)
  if sec_num == 1
    if wafer_num == 1 println("Error loading 1,1-prealigned.jld - the first section is the identity"); return Void
  else firstindex = MONTAGED_OFFSETS[findlast(i->MONTAGED_OFFSETS[i,2][1] == wafer_num -1, 1:size(MONTAGED_OFFSETS, 1)), 2]
  end
  else firstindex = (wafer_num, sec_num-1, MONTAGED_INDEX, MONTAGED_INDEX)
  end
	return load(firstindex, lastindex)
end

function load(firstindex::Index, lastindex::Index)
  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex))
    filename = joinpath(PREALIGNED_DIR, string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_prealigned.jld"))
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_prealigned(lastindex))
    filename = joinpath(ALIGNED_DIR, string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned.jld"))
  else 
    filename = joinpath(MONTAGED_DIR, string(join(firstindex[1:2], ","), "_montaged.jld"))
  end

  println("Loading meshset from", filename)
  return load(filename)
  end

function load(filename::String)
  Ms = JLD.load(filename, "MeshSet"); 
  return Ms
end

function get_all_overlaps(Ms)
adjacent_pairs = Pairings(0)
diagonal_pairs = Pairings(0)

  for i in 1:Ms.N, j in 1:Ms.N
    if is_adjacent(Ms.meshes[i], Ms.meshes[j]) push!(adjacent_pairs, (i, j)); end
    if is_diagonal(Ms.meshes[i], Ms.meshes[j]) push!(diagonal_pairs, (i, j)); end
  end

  pairs = vcat(adjacent_pairs, diagonal_pairs)

  return pairs
end

function add_pair_matches_reflexive!(Ms, src_index, dst_index, images = Void)
  if images == Void
  images = load_section_pair(Ms, src_index, dst_index)
  end
  matches_src_dst = Matches(images[1], Ms.meshes[find_section(Ms,src_index)], 
                              images[2], Ms.meshes[find_section(Ms,dst_index)], 
                              Ms.params)
  matches_dst_src = Matches(images[2], Ms.meshes[find_section(Ms,dst_index)], 
                              images[1], Ms.meshes[find_section(Ms,src_index)], 
                              Ms.params)
  add_matches(matches_src_dst, Ms)
  add_matches(matches_dst_src, Ms)
  return Ms
end

function add_pair_matches!(Ms, src_index, dst_index, images = Void)
  if images == Void
  images = load_section_pair(Ms, src_index, dst_index)
  end
  matches = Matches(images[1], Ms.meshes[find_section(Ms,src_index)], 
                              images[2], Ms.meshes[find_section(Ms,dst_index)], 
                              Ms.params)
  add_matches(matches, Ms)
  return Ms
end

"""
Include prealignment review image with prealignment process for faster review
"""
function add_pair_matches_with_thumbnails!(meshset, src_index, dst_index, images = Void)
  if images == Void
  images = load_section_pair(meshset, src_index, dst_index)
  end
  src_mesh = meshset.meshes[find_section(meshset, src_index)]
  dst_mesh = meshset.meshes[find_section(meshset, dst_index)]
  matches = calculate_matches(images..., src_mesh, dst_mesh, meshset.params)
  add_matches(matches, meshset)  
  write_prealignment_thumbnail(images..., meshset)
end

function calculate_matches(src_img, dst_img, src_mesh, dst_mesh, params)
  return Matches(src_img, src_mesh, dst_img, dst_mesh, params)
end

function add_all_matches!(Ms, images)

pairs = get_all_overlaps(Ms)
n = length(pairs)
i = 1
nextidx() = (idx=i; i+=1; idx)
matches_array = cell(n)

        while true
          idx = nextidx()
            if idx > n
              break
            end
          (a, b) = pairs[idx]
          matches_array[idx] = Matches(images[a], Ms.meshes[a], images[b], Ms.meshes[b], Ms.params)
        end
for k in 1:n
    M = matches_array[k]
    if typeof(M) == Void || M == Void continue; end
    add_matches(M, Ms)
end
  return Ms
end

function get_matched_points(Ms::MeshSet, k)
  src_mesh = Ms.meshes[find_index(Ms, Ms.matches[k].src_index)]
  src_indices = Ms.matches[k].src_points_indices
  src_pts = src_mesh.nodes[src_indices]
  dst_pts = Ms.matches[k].dst_points
  return src_pts, dst_pts
end

function get_matched_points(Ms::MeshSet)
  src_points = Points(0); 
  dst_points = Points(0); 
  for k in 1:Ms.M
    pts = get_matched_points(Ms, k)
    src_points = vcat(src_points,pts[1])
    dst_points = vcat(dst_points,pts[2])
  end
  return src_points, dst_points
end

function get_matched_points_t(Ms::MeshSet, k)
  src_pts = Points(0)
  dst_pts = Points(0)

  for i in 1:Ms.matches[k].n
      w = Ms.matches[k].dst_weights[i]
      t = Ms.matches[k].dst_triangles[i]
      p = Ms.matches[k].src_points_indices[i]
      src = Ms.meshes[find_index(Ms, Ms.matches[k].src_index)]
      dst = Ms.meshes[find_index(Ms, Ms.matches[k].dst_index)]
      p1 = src.nodes_t[p]
      p2 = dst.nodes_t[t[1]] * w[1] + dst.nodes_t[t[2]] * w[2] + dst.nodes_t[t[3]] * w[3]
      push!(src_pts, p1)
      push!(dst_pts, p2)
  end
  return src_pts, dst_pts
end

function get_matched_points_t(Ms::MeshSet)
  src_points = Points(0); 
  dst_points = Points(0); 
  for k in 1:Ms.M
    pts = get_matched_points_t(Ms, k)
    src_points = vcat(src_points,pts[1])
    dst_points = vcat(dst_points,pts[2])
   end
  return src_points, dst_points
end

function load_section(offsets, wafer_num, section_num)
  indices = find(i -> offsets[i,2][1:2] == (wafer_num, section_num), 1:size(offsets, 1))
  Ms = MeshSet(PARAMS_MONTAGE)
  num_tiles = length(indices)
  paths = Array{String, 1}(num_tiles)

# images = Array{SharedArray{UInt8, 2}, 1}(0)
  images = Array{Array{UInt8, 2}, 1}(0)


  for i in indices
    name = offsets[i, 1];
    index = offsets[i, 2];
    dy = offsets[i, 3] #/ 0.07; ##################################
    dx = offsets[i, 4] #/ 0.07; ##################################
    image = get_image(get_path(name));
    add_mesh(Mesh(name, image, index, dy, dx, false, PARAMS_MONTAGE), Ms);
    #image_shared = SharedArray(UInt8, size(image, 1), size(image, 2));
    #image_shared[:, :] = image[:, :];
    push!(images, image)
  end

  return Ms, images
end

"""
`[method]_APPROXIMATE` - apply transform to the Mesh type
`[method]_SOLVE` - apply transform to the Matches type
"""

"""
Retrieve nodes and nodes_t of mesh, make homogenous, and transpose to Nx3
"""
function get_homogeneous_correspondences(M::Mesh)
  pts = M.nodes
  pts_t = M.nodes_t
  num_pts = size(pts, 1)
  hpts = Array{Float64, 2}(3, num_pts)
  hpts_t = Array{Float64, 2}(3, num_pts)
  for i in 1:num_pts
    hpts[:, i] = [pts[i]; 1]
    hpts_t[:, i] = [pts_t[i]; 1]
  end
  return hpts', hpts_t'
end

"""
Return right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""
function rigid_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_correspondences(M)
  return find_rigid(pts_src, pts_dst)
end

"""
Return the right-hand matrix for the mesh
`nodes*T ~= nodes_T`
"""
function affine_approximate(M::Mesh)
  pts_src, pts_dst = get_homogeneous_correspondences(M)
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
Retrieve nodes and nodes_t of mesh, make homogenous, and transpose to Nx3
"""
function get_homogeneous_correspondences(Ms, k)
  pts_src, pts_dst = get_matched_points(Ms, k)
  num_pts = size(pts_src, 1)

  hpts_src = Array{Float64, 2}(3, num_pts)
  hpts_dst = Array{Float64, 2}(3, num_pts)

  for i in 1:num_pts
    hpts_src[:, i] = [pts_src[i]; 1]
    hpts_dst[:, i] = [pts_dst[i]; 1]
  end
  return hpts_src', hpts_dst'
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

function make_stack(offsets, wafer_num, batch::UnitRange{Int64})

  indices = find(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] in batch, 1:size(offsets, 1))
  Ms = MeshSet(PARAMS_ALIGNMENT)
#=
  index_aligned = (wafer_num, a, ALIGNED_INDEX, ALIGNED_INDEX)
  name_aligned = get_name(index_aligned);  
  dy_aligned = 0
  dx_aligned = 0
  size_i = 36000; 
  size_j = 36000

  add_mesh(Mesh(name_aligned, size_i, size_j, index_aligned, dy_aligned, dx_aligned, true, PARAMS_ALIGNMENT), Ms)
=#
  dy = 0
  dx = 0

  for i in indices
    name = offsets[i, 1]
    index = offsets[i, 2]
    dy += offsets[i, 3]
    dx += offsets[i, 4]
    size_i = offsets[i, 5]
    size_j = offsets[i, 6]
    add_mesh(Mesh(name, size_i, size_j, index, dy, dx, false, PARAMS_ALIGNMENT), Ms)
  end
  optimize_all_cores(Ms.params)
  return Ms
end

function affine_make_stack(offsets, wafer_num, a::Int64, b::Int64, optimize = true)
  i_dst = findfirst(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] == a, 1:size(offsets, 1))
  i_src = findfirst(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] == b, 1:size(offsets, 1))
  Ms = MeshSet(PARAMS_PREALIGNMENT)

  name_dst = offsets[i_dst, 1];
  index_dst = offsets[i_dst, 2];
  dy_dst = 0;#offsets[i_dst, 3];
  dx_dst = 0;#offsets[i_dst, 4];
  size_i = offsets[i_dst, 5]; 
  size_j = offsets[i_dst, 6]

  add_mesh(Mesh(name_dst, size_i, size_j, index_dst, dy_dst, dx_dst, true, PARAMS_PREALIGNMENT), Ms)

  name = offsets[i_src, 1];
  index = offsets[i_src, 2];
  dy = offsets[i_src, 3];
  dx = offsets[i_src, 4]; 
  size_i = offsets[i_src, 5];
  size_j = offsets[i_src, 6];

  add_mesh(Mesh(name, size_i, size_j, index, dy, dx, false, PARAMS_PREALIGNMENT), Ms)
  if optimize == true
  optimize_all_cores(Ms.params)
  end

  return Ms
end
function make_stack(offsets, wafer_num, a::Int64, b::Int64)
  i = findfirst(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] == b, 1:size(offsets, 1))
  Ms = MeshSet(PARAMS_ALIGNMENT)

  index_aligned = (wafer_num, a, ALIGNED_INDEX, ALIGNED_INDEX)
  name_aligned = get_name(index_aligned);  
  dy_aligned = 0
  dx_aligned = 0
  size_i = 36000; 
  size_j = 36000

  add_mesh(Mesh(name_aligned, size_i, size_j, index_aligned, dy_aligned, dx_aligned, true, PARAMS_ALIGNMENT), Ms)
  name = offsets[i, 1]
  index = offsets[i, 2]
  dy = offsets[i, 3]
  dx = offsets[i, 4]
  size_i = offsets[i, 5]
  size_j = offsets[i, 6]

  add_mesh(Mesh(name, size_i, size_j, index, dy, dx, false, PARAMS_ALIGNMENT), Ms)
  optimize_all_cores(Ms.params)

  return Ms
end

function make_stack(offsets, wafer_num, section_range, fixed_interval)
  indices = find(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] in section_range, 1:size(offsets, 1))
  Ms = MeshSet(PARAMS_ALIGNMENT)

  for i in indices
    name = offsets[i, 1];
    index = offsets[i, 2];
    dy = offsets[i, 3];
    dx = offsets[i, 4];
    size_i = offsets[i, 5]
    size_j = offsets[i, 6]
    is_fixed = false
    if findfirst(indices, i) in 1:fixed_interval:length(indices)
      is_fixed = true; println("$index is fixed")
    end
    add_mesh(Mesh(name, size_i, size_j, index, dy, dx, is_fixed, PARAMS_ALIGNMENT), Ms)
  end

  optimize_all_cores(Ms.params)

  return Ms
end

function crop_center(image, rad_ratio)
	size_i, size_j = size(image);
	rad = round(Int64, rad_ratio * (minimum([size_i, size_j]) / 2));
	range_i = round(Int64, size_i / 2) + (-rad:rad);
	range_j = round(Int64, size_j / 2) + (-rad:rad);
	return image[range_i, range_j];
end

function affine_load_section_pair(offsets, wafer_num, src, dst)
  i_dst = findfirst(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] == dst, 1:size(offsets, 1))
  i_src = findfirst(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] == src, 1:size(offsets, 1))

  name_dst = offsets[i_dst, 1];
  name_src = offsets[i_src, 1];

  @time dst_image = get_image(get_path(name_dst))
  @time src_image = get_image(get_path(name_src))
 
  dst_scaled = imwarp(dst_image, SCALING_FACTOR_TRANSLATE)[1]; 
  src_scaled = imwarp(src_image, SCALING_FACTOR_TRANSLATE)[1]; 
 
  src_cropped = crop_center(src_scaled, 0.5);

  offset_vect, xc = get_max_xc_vector(src_cropped, dst_scaled);

  offset_unscaled = round(Int64, offset_vect[1:2] / SCALING_FACTOR_TRANSLATE);

  println("Offsets from scaled blockmatches: $offset_unscaled");
  println("r: $(offset_vect[3])");
  update_offsets(name_src, offset_unscaled);
  return src_image, dst_image;
end

function load_section_pair(Ms, a, b)
  @time A_image = get_image(get_path(Ms.meshes[find_section(Ms,a)].name))
  @time B_image = get_image(get_path(Ms.meshes[find_section(Ms,b)].name))
      return A_image, B_image; 
end

function load_stack(offsets, wafer_num, section_range)
  indices = find(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] in section_range, 1:size(offsets, 1))
  Ms = MeshSet(PARAMS_ALIGNMENT)
  images = Array{SharedArray{UInt8, 2}, 1}(0)

  for i in indices
  name = offsets[i, 1]
  index = offsets[i, 2]
  dx = offsets[i, 4]
  dy = offsets[i, 3]
  image = get_image(get_path(name))
  add_mesh(Mesh(name, image, index, dy, dx, false, PARAMS_ALIGNMENT), Ms)
  
  image_shared = SharedArray(UInt8, size(image, 1), size(image, 2))
  image_shared[:, :] = image[:, :]
  push!(images, image_shared)
  end

  optimize_all_cores(Ms.params)

  return Ms, images
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

   println("Residuals before solving: rms: $rms,  mean: $avg, sigma = $sig, max = $max\n")
   println("Residuals after solving: rms: $rms_t,  mean: $avg_t, sigma = $sig_t, max = $max_t\n")
 # decomp_affine(affine_solve(Ms))
end
