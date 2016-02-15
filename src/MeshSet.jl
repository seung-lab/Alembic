type MeshSet
 	meshes::Array{Mesh, 1}			# vector of meshes in the set
 	matches::Array{Match, 1}		# vector of matches in the set
  	properties::Dict{Any, Any}		# in the same array order, contains parameters and filters
end

### IO.jl
function get_path(meshset::MeshSet)			return get_path(meshset.meshes);		end

function get_images(meshset::MeshSet, dtype = UInt8)	return map(get_image, meshset.meshes, repeated(dtype));	end

function get_correspondence_patches(meshset::MeshSet, match_ind, corr_ind) return get_correspondence_patches(meshset.matches[match_ind], corr_ind)	end

function get_params(meshset::MeshSet)			return meshset.properties["params"];		end

function get_fixed(meshset::MeshSet)			return meshset.properties["fixed"];		end

### counting
function count_meshes(meshset::MeshSet)			return length(meshset.meshes);		end

function count_matches(meshset::MeshSet)		return length(meshset.matches);		end

function count_nodes(meshset::MeshSet)			return sum(map(count_nodes, meshset.meshes));		end

function count_edges(meshset::MeshSet)			return sum(map(count_edges, meshset.meshes));		end

function count_correspondences(meshset::MeshSet)	return sum(map(count_correspondences, meshset.matches));		end

function count_filtered_correspondences(meshset::MeshSet)	return sum(map(count_filtered_correspondences, meshset.matches));		end

### finding
function find_mesh_index(meshset::MeshSet, index)
  return findfirst(this -> index == this.index, meshset.meshes)
  end

function find_match_index(meshset::MeshSet, src_index, dst_index)
  return findfirst(this -> (src_index == this.src_index) && (dst_index == this.dst_index), meshset.matches)
  end

function find_node(Ms, mesh_ind, node_ind)
  return Ms.nodes_indices[mesh_ind] + node_ind
end

### adding
function add_mesh!(meshset::MeshSet, mesh::Mesh)
  push!(meshset.meshes, mesh);
end

function add_match!(meshset::MeshSet, match::Match)
  push!(meshset.matches, match);
end

### adding
function filter!(meshset::MeshSet, property_name, compare, threshold)
	total = 0;
	for match in meshset.matches
		total = total + filter!(match, property_name, compare, threshold)
	end
	println("$total / $(count_correspondences(meshset)) matches filtered")
end

### initialise
function MeshSet()
 	meshes = Array{Mesh, 1}(0)
 	matches = Array{Match, 1}(0)		
  	properties = Dict{Any, Any}()

	return MeshSet(meshes, matches, properties)
end

function prealign(index; params=get_params(index))
	src_index = index;
	dst_index = get_preceding(src_index)
	meshset = MeshSet();
	meshset.properties["params"] = params;
	push!(meshset.meshes, Mesh(src_index, params))
	push!(meshset.meshes, Mesh(dst_index, params))
	push!(meshset.matches, Match(meshset.meshes[1], meshset.meshes[2], params))
	solve!(meshset, method=params["solve"]["method"]);
	save(meshset);
	return meshset;
end

function MeshSet(index; params=get_params(index))
	if is_premontaged(index) return MeshSet(index, index); end
	if is_montaged(index) return MeshSet(index, index); end
end

function MeshSet(first_index, last_index; params=get_params(first_index), solve_method="elastic", fix_first=false)
	ind_range = get_index_range(first_index, last_index);
	if length(ind_range) == 0 return 0; end
	fixed = Array{Any, 1}(0);
	if fix_first push!(fixed, first_index); end
	meshes = map(Mesh, ind_range, repeated(params))
 	matches = Array{Match, 1}(0)		
	properties = Dict(	"params"  => params,
				"by"	  => ENV["USER"],
				"machine" => gethostname(),
				"timestamp" => string(now()),
				"fixed" => fixed
				)
	meshset = MeshSet(meshes, matches, properties);
	match!(meshset);
	solve!(meshset, method=solve_method);
	save(meshset);
	return meshset;
end


### match
function get_all_overlaps(meshset::MeshSet)	return get_all_overlaps(meshset.meshes);	end;
function get_all_overlaps(meshes::Array{Mesh, 1})
adjacent_pairs = Pairings(0)
diagonal_pairs = Pairings(0)
preceding_pairs = Pairings(0)
succeeding_pairs = Pairings(0)

  for i in 1:length(meshes), j in 1:length(meshes)
    if is_adjacent(meshes[i], meshes[j]) push!(adjacent_pairs, (i, j)); end
    if is_diagonal(meshes[i], meshes[j]) push!(diagonal_pairs, (i, j)); end
    if is_preceding(meshes[i], meshes[j]) 
    	push!(preceding_pairs, (i, j)); 
    	push!(succeeding_pairs, (j, i)); 
    end
  end

  pairs = vcat(adjacent_pairs, diagonal_pairs, preceding_pairs, succeeding_pairs)

  return pairs
end

function match!(meshset::MeshSet; prefetch_all = false)
	params = get_params(meshset);
	pairs = get_all_overlaps(meshset);
	if prefetch_all
	imgdict = Dict{Mesh, Array}();
	for mesh in meshset.meshes
		imgdict[mesh] = get_image(mesh);
	end
	for i in 1:length(pairs)
		add_match!(meshset, Match(meshset.meshes[pairs[i][1]], meshset.meshes[pairs[i][2]], src_image=imgdict[meshset.meshes[pairs[i][1]]], dst_image=imgdict[meshset.meshes[pairs[i][2]]], params));
	end
	else
	for i in 1:length(pairs)
		add_match!(meshset, Match(meshset.meshes[pairs[i][1]], meshset.meshes[pairs[i][2]], params));
	end
	end
end

function sanitize!(meshset::MeshSet)
  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
  	meshes[mesh.index] = mesh;
  end

  for match in meshset.matches
    	src_mesh = meshes[match.src_index];
    	dst_mesh = meshes[match.dst_index];
	src_pts, dst_pts = get_filtered_correspondences(match);
	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	invalids = union(find(ind -> src_pt_triangles[ind] == NO_TRIANGLE, 1:count_filtered_correspondences(match)), find(ind -> dst_pt_triangles[ind] == NO_TRIANGLE, 1:count_filtered_correspondences(match)))
	if length(invalids) !=0
	filter!(match; inds = invalids, filtertype = "sanitization");
	end
  end
end

# JLS SAVE
function save(filename::String, meshset::MeshSet)
  println("Saving meshset to ", filename)
  open(filename, "w") do file
    serialize(file, meshset)
  end
end

function save(meshset::MeshSet)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[count_meshes(meshset)].index

  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex))
    filename = joinpath(PREALIGNED_DIR, string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_prealigned.jls"))
    #update_offset(prealigned(firstindex), [0, 0]);
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_prealigned(lastindex))
    filename = joinpath(ALIGNED_DIR, string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned.jls"))
    #update_offset(aligned(firstindex), [0, 0]);
  else 
    filename = joinpath(MONTAGED_DIR, string(join(firstindex[1:2], ","), "_montaged.jls"))
    #update_offset(montaged(firstindex), [0, 0]);
  end
  save(filename, meshset);
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

function load(firstindex, lastindex)
  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex))
    filename = joinpath(PREALIGNED_DIR, string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_prealigned.jls"))
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_prealigned(lastindex))
    filename = joinpath(ALIGNED_DIR, string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned.jls"))
  else 
    filename = joinpath(MONTAGED_DIR, string(join(firstindex[1:2], ","), "_montaged.jls"))
  end

  println("Loading meshset from ", filename)
  return load(filename)
end

function load(index)
  	if is_montaged(index)
	  return load_montaged(index[1:2]...)
	end
end

function load(filename::String)
  return open(deserialize, filename)
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
