"""

"""
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
#= JLD SAVE
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

  save(filename, Ms);
end
=#

# JLS SAVE
# function save(filename::String, Ms::MeshSet)
#   println("Saving meshset to ", filename)
#   open(filename, "w") do file
#     serialize(file, Ms)
#   end
# end

function save(filename::String, Ms::MeshSet)
  println("Saving meshset to ", filename)
  jldopen(filename, "w") do file
    write(file, "MeshSet", Ms)
  end
end

function save(Ms::MeshSet)
  firstindex = Ms.meshes[1].index
  lastindex = Ms.meshes[Ms.N].index
  filename = get_name(firstindex, lastindex)
  save(filename, Ms);
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

function get_name(firstindex::Index, lastindex::Index, ext="jld")
  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex))
    filename = joinpath(PREALIGNED_DIR, string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_prealigned.", ext))
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_aligned(lastindex))
    filename = joinpath(ALIGNED_DIR, string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned.", ext))
  else 
    filename = joinpath(MONTAGED_DIR, string(join(firstindex[1:2], ","), "_montaged.", ext))
  end
  return filename
end

function load(firstindex::Index, lastindex::Index)
  filename = get_name(firstindex, lastindex)
  println("Loading meshset from ", filename)
  return load(filename)
  # return open(deserialize, filename)
end

function load_aligned(firstindex::Index, lastindex::Index)
  filename = get_name(firstindex, lastindex, "jls")
  println("Loading meshset from ", filename)
  # return load(filename)
  return open(deserialize, filename)
end

function load(filename::String)
  return JLD.load(filename, "MeshSet"); 
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

function add_pair_matches_reflexive!(Ms, src, dst, images = Void)
  if images == Void
  images = load_section_pair(Ms, src, dst)
  end
  matches_src_dst = Matches(images[1], Ms.meshes[find_index(Ms,src)], 
                              images[2], Ms.meshes[find_index(Ms,dst)], 
                              Ms.params)
  matches_dst_src = Matches(images[2], Ms.meshes[find_index(Ms,dst)], 
                              images[1], Ms.meshes[find_index(Ms,src)], 
                              Ms.params)
  add_matches(matches_src_dst, Ms)
  add_matches(matches_dst_src, Ms)
  return Ms
end

function add_pair_matches!(Ms, src, dst, images = Void)
  if images == Void
  images = load_section_pair(Ms, src, dst)
  end
  matches = Matches(images[1], Ms.meshes[find_index(Ms,src)], 
                              images[2], Ms.meshes[find_index(Ms,dst)], 
                              Ms.params)
  add_matches(matches, Ms)
  return Ms
end

"""
Include prealignment review image with prealignment process for faster review
"""
function add_pair_matches_with_thumbnails!(meshset, src, dst, images = Void)
  if images == Void
  images = load_section_pair(meshset, src, dst)
  end
  src_mesh = meshset.meshes[find_index(meshset, src)]
  dst_mesh = meshset.meshes[find_index(meshset, dst)]
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

function make_stack(first_index, last_index, fixed_interval = 0)
	if first_index[3] != last_index[3] || first_index[4] != last_index[4]
		println("The indices are from different pipeline stages. Aborting.");
		return Void;
	end

	registry = get_registry(first_index); params = get_params(first_index); indices = get_range_in_registry(first_index, last_index);

	Ms = MeshSet(params);
	dy = 0; dx = 0;

  for i in indices
    name = registry[i, 1]
    index = registry[i, 2]
    if params["global_offsets"] == true
    	dy = registry[i, 3]
    	dx = registry[i, 4]
    elseif i == maximum(indices)
    	dy += registry[i, 3]
    	dx += registry[i, 4]
    end
    size_i = registry[i, 5]
    size_j = registry[i, 6]
    add_mesh(Mesh(name, size_i, size_j, index, dy, dx, false, params), Ms)
  end
  optimize_all_cores(Ms.params)
  return Ms;
end

function crop_center(image, rad_ratio)
	size_i, size_j = size(image);
	rad = round(Int64, rad_ratio * (minimum([size_i, size_j]) / 2));
	range_i = round(Int64, size_i / 2) + (-rad:rad);
	range_j = round(Int64, size_j / 2) + (-rad:rad);
	return image[range_i, range_j];
end

function affine_load_section_pair(src_index, dst_index)
  i_src = find_in_registry(src_index); 
  i_dst = find_in_registry(dst_index); 

  registry = get_registry(src_index);

  name_dst = registry[i_dst, 1];
  name_src = registry[i_src, 1];

  @time dst_image = get_uint8_image(get_path(name_dst))
  @time src_image = get_uint8_image(get_path(name_src))
 
  dst_scaled = imscale(dst_image, SCALING_FACTOR_TRANSLATE)[1]; 
  src_scaled = imscale(src_image, SCALING_FACTOR_TRANSLATE)[1]; 

  src_cropped = crop_center(src_scaled, 0.33);
  dst_cropped = crop_center(dst_scaled, 0.66);
  offset_vect, xc = get_max_xc_vector(src_cropped, dst_cropped);

  offset_unscaled = round(Int64, offset_vect[1:2] / SCALING_FACTOR_TRANSLATE);

  view(xc * 40);

  println(offset_vect[1:2]);
  println("Offsets from scaled blockmatches: $offset_unscaled");
  println("r: $(offset_vect[3])");
  update_offsets(name_src, offset_unscaled); 
  return src_image, dst_image;
end

function load_section_pair(Ms, a, b)
  @time A_image = get_h5_image(get_h5_path(Ms.meshes[find_index(Ms,a)].index))
  @time B_image = get_h5_image(get_h5_path(Ms.meshes[find_index(Ms,b)].index))
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

"""
Calculate the maximum bounding box of all the meshes in a meshset
"""
function get_global_bb(meshset)
  bbs = []
  println("Calculating global bounding box")
  for mesh in meshset.meshes
      nodes = hcat(mesh.nodes_t...)'
      push!(bbs, snap_bb(find_mesh_bb(nodes)))
  end
  global_bb = sum(bbs)
  global_bb.h += 1
  global_bb.w += 1
  println(global_bb)
  return global_bb
end  
