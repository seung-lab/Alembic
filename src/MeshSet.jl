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

#function get_fixed(meshset::MeshSet)			return meshset.properties["fixed"];		end

### counting
function count_meshes(meshset::MeshSet)			return length(meshset.meshes);		end

function count_matches(meshset::MeshSet)		return length(meshset.matches);		end

function count_nodes(meshset::MeshSet)			return sum(map(count_nodes, meshset.meshes));		end

function count_edges(meshset::MeshSet)			return sum(map(count_edges, meshset.meshes));		end

function count_correspondences(meshset::MeshSet)	return sum(map(count_correspondences, meshset.matches));		end

function count_filtered_correspondences(meshset::MeshSet)	return sum(map(count_filtered_correspondences, meshset.matches));		end

### finding
function find_mesh_index(meshset::MeshSet, index)
  return findfirst(this -> index == get_index(this), meshset.meshes)
  end

function find_match_index(meshset::MeshSet, src_index, dst_index)
  return findfirst(this -> (src_index == get_src_index(this)) && (dst_index == get_dst_index(this)), meshset.matches)
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



### reviewing
function set_reviewed!(meshset::MeshSet, match_ind, flag = false)
	set_reviewed!(meshset.matches[match_ind], flag);
end

function flag!(meshset::MeshSet, match_ind)
	flag!(meshset.matches[match_ind])
end

function unflag!(meshset::MeshSet, match_ind)
	unflag!(meshset.matches[match_ind])
end

function is_reviewed(meshset::MeshSet, match_ind)
  return is_reviewed(meshset.matches[match_ind])
end

function is_reviewed(meshset::MeshSet)
  return (&)(map(is_reviewed, meshset.matches)...)
end

function is_flagged(meshset::MeshSet, match_ind)
	return is_flagged(meshset.matches[match_ind])
end

function is_flagged(meshset::MeshSet)
	return |(map(is_flagged, meshset.matches)...)
end

function count_flags(meshset::MeshSet)
  return sum(map(is_flagged, meshset.matches))
end

function filter!(meshset::MeshSet, filters = Base.values(meshset.properties["params"]["filter"]))
  	for filter in filters
		filter!(meshset, filter...)
	end
end

function filter!(meshset::MeshSet, property_name, compare, threshold)
	total = 0;
	for match in meshset.matches
		total = total + filter!(match, property_name, compare, threshold)
	end
	println("$total / $(count_correspondences(meshset)) correspondences filtered on $property_name")
end

function check!(meshset::MeshSet, crits = Base.values(meshset.properties["params"]["review"])) 
  return |(map(check!, meshset.matches, repeated(crits))...)
end

function check_and_fix!(meshset::MeshSet, crits = Base.values(meshset.properties["params"]["review"]), filters = Base.values(meshset.properties["params"]["filter"])) 
  if |(map(check!, meshset.matches, repeated(crits))...)
    for match in meshset.matches
      if is_flagged(match)
      	clear_filters!(match);
      	filter!(match, filters)
      end
    end
    fixed = !check!(meshset, crits);
    if fixed
      println("fixed successfully");
    else println("failed to fix meshset")
    for match in meshset.matches
      if is_flagged(match)
      	clear_filters!(match);
      end
    end
    end
  end
end

### splitting
function split_meshset(meshset::MeshSet)
	parent_name = get_name(meshset)

	for i in 1:count_matches(meshset)
		properties = deepcopy(meshset.properties);
		properties["meta"]["parent"] = parent_name;
		properties["meta"]["split_index"] = i;
		matches = Array{Match, 1}();
		push!(matches, meshset.matches[i])
		meshes = Array{Mesh, 1}();
		src_mesh = meshset.meshes[find_mesh_index(meshset, get_src_index(meshset.matches[i]))]
		dst_mesh = meshset.meshes[find_mesh_index(meshset, get_dst_index(meshset.matches[i]))]
		push!(meshes, src_mesh, dst_mesh)
		save(MeshSet(meshes, matches, properties));
		println("Child ", i, " / ", count_matches(meshset), " saved");
	end
end

function concat_meshset(parent_name)
	ms = MeshSet();
	for i in 1:count_children(parent_name)
		cms = load_split(parent_name, i)
		append!(ms.matches, cms.matches)
		for cmesh in cms.meshes
			ind = find_mesh_index(ms,get_index(cmesh))
			if ind == 0 push!(ms.meshes, cmesh) end
		end

		if i == 1
			ms.properties = deepcopy(cms.properties);
			ms.properties["meta"]["parent"] = nothing;
			ms.properties["meta"]["split_index"] = 0;

		end
		println("Child ", i, " / ", count_children(parent_name), " concatanated");
	end
	return ms;
end

function concat!(meshset_one::MeshSet, meshset_two::MeshSet)
		append!(meshset_one.matches, meshset_two.matches)
		for mesh_two in meshset_two.meshes
			ind = find_mesh_index(ds, get_index(mesh_two))
			if ind == 0 push!(meshset_one.meshes, mesh_two) end
		end
	return meshset_one;
end

function flag!(meshset::MeshSet, match_ind)
	flag!(meshset.matches[match_ind])
end

function unflag!(meshset::MeshSet, match_ind)
	unflag!(meshset.matches[match_ind])
end
function is_flagged(meshset::MeshSet, match_ind)
	return is_flagged(meshset.matches[match_ind])
end

function is_flagged(meshset::MeshSet)
	return |(map(is_flagged, meshset.matches)...)
end

### initialise
function MeshSet()
 	meshes = Array{Mesh, 1}(0)
 	matches = Array{Match, 1}(0)		
  	properties = Dict{Any, Any}()

	return MeshSet(meshes, matches, properties)
end

function prealign(index; params=get_params(index), to_fixed=false)
	src_index = index;
	dst_index = get_preceding(src_index)
	if to_fixed
	dst_index = aligned(dst_index);
	end
	meshset = MeshSet();
	meshset.properties["params"] = params;
	push!(meshset.meshes, Mesh(src_index, params))
	push!(meshset.meshes, Mesh(dst_index, params, to_fixed))
	push!(meshset.matches, Match(meshset.meshes[1], meshset.meshes[2], params))
	filter!(meshset);
	check_and_fix!(meshset);
	solve!(meshset, method=params["solve"]["method"]);
	save(meshset);
	return meshset;
end

function MeshSet(index; params=get_params(index))
	if is_premontaged(index) return MeshSet(index, index); end
	if is_montaged(index) return MeshSet(premontaged(index), premontaged(index); prefetch_all=true); end
end

function MeshSet(first_index, last_index; params=get_params(first_index), solve=true, solve_method="elastic", fix_first=false, prefetch_all = false)
	ind_range = get_index_range(first_index, last_index);
	if length(ind_range) == 0 return nothing; end
	fixed_inds = Array{Any, 1}(0);
	if fix_first push!(fixed_inds, first_index); 
		ind_range[1] = aligned(ind_range[1])
	end
	meshes = map(Mesh, ind_range, repeated(params), map(in, ind_range, repeated(fixed_inds)))
 	matches = Array{Match, 1}(0)		
	properties = Dict{Any, Any}(	"params"  => params,
					"author" => author(),
					"meta" => Dict{Any, Any}(
					"parent" => nothing,
					"split_index" => 0)
					)
	meshset = MeshSet(meshes, matches, properties);
	match!(meshset, params["match"]["depth"]; prefetch_all=prefetch_all);

	filter!(meshset);
	check_and_fix!(meshset);
#=	
	if check!(meshset)
		save(meshset); return meshset;
	end =#

	if solve == true
	solve!(meshset, method=solve_method);
	end

	save(meshset);
	return meshset;
end




### match
function get_all_overlaps(meshset::MeshSet, within = 1)	return get_all_overlaps(meshset.meshes, within);	end;
function get_all_overlaps(meshes::Array{Mesh, 1}, within = 1)
adjacent_pairs = Pairings(0)
diagonal_pairs = Pairings(0)
preceding_pairs = Pairings(0)
succeeding_pairs = Pairings(0)

  for i in 1:length(meshes), j in 1:length(meshes)
    if is_adjacent(get_index(meshes[i]), get_index(meshes[j])) push!(adjacent_pairs, (i, j)); end
    if is_diagonal(get_index(meshes[i]), get_index(meshes[j])) push!(diagonal_pairs, (i, j)); end
    if is_preceding(get_index(meshes[i]), get_index(meshes[j]), within) 
    	push!(preceding_pairs, (i, j)); 
    	push!(succeeding_pairs, (j, i)); 
    end
  end

  pairs = vcat(adjacent_pairs, diagonal_pairs, preceding_pairs, succeeding_pairs)
  println("$(length(pairs)) pairs found")
  return pairs
end

function match!(meshset::MeshSet, within = 1; prefetch_all = false)
	params = get_params(meshset);
	pairs = get_all_overlaps(meshset, within);
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

function rematch!(meshset::MeshSet, match_ind, params)
  src_index = get_src_index(meshset.matches[match_ind])
  dst_index = get_dst_index(meshset.matches[match_ind])
  src_mesh = meshset.meshes[find_mesh_index(meshset, src_index)]
  dst_mesh = meshset.meshes[find_mesh_index(meshset, dst_index)]
  deleteat!(meshset.matches, match_ind)
  add_match!(meshset, Match(src_mesh, dst_mesh, params))
  meshset.properties["params"] = params
  meshset.properties["author"] = author()
  save(meshset)
end

function sanitize!(meshset::MeshSet)
  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
  	meshes[get_index(mesh)] = mesh;
  end

  for match in meshset.matches
    	src_mesh = meshes[match.src_index];
    	dst_mesh = meshes[match.dst_index];
	src_pts, dst_pts = get_correspondences(match);
	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), src_pts);
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), dst_pts);
	invalids = union(find(ind -> src_pt_triangles[ind] == NO_TRIANGLE, 1:count_correspondences(match)), find(ind -> dst_pt_triangles[ind] == NO_TRIANGLE, 1:count_correspondences(match)))
	if length(invalids) !=0
	clear_filters!(match; filtertype = "sanitization");
	filter_manual!(match, invalids; filtertype = "sanitization");
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
  filename = get_filename(meshset)
  save(filename, meshset);
end

function get_filename(meshset::MeshSet)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[count_meshes(meshset)].index
  return get_filename(firstindex, lastindex)
end

function get_filename(firstindex::Index, lastindex::Index)
  filename = string(get_name(firstindex, lastindex), ".jls")
  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_aligned(lastindex))
    filepath = PREALIGNED_DIR
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_prealigned(lastindex))
    filepath = ALIGNED_DIR
  else 
    filepath = MONTAGED_DIR
  end

  return joinpath(filepath, filename)
end

#### ONLY AS PARENT!!!!
function get_name(meshset::MeshSet)
  if has_parent(meshset)
    return get_split_index(meshset);
  else
    firstindex = meshset.meshes[1].index
    lastindex = meshset.meshes[count_meshes(meshset)].index
    return get_name(firstindex, lastindex)
  end
end

function get_name(firstindex::Index, lastindex::Index)
  name = ""
  if (is_prealigned(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_montaged(lastindex)) || (is_montaged(firstindex) && is_aligned(lastindex))
    name = string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_prealigned")
  elseif (is_prealigned(firstindex) && is_prealigned(lastindex)) || (is_aligned(firstindex) && is_prealigned(lastindex))
    name = string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned")
  else 
    name = string(join(firstindex[1:2], ","), "_montaged")
  end
	return name
end

function load_split(parent_name, split_index)
    filename = joinpath(ALIGNED_DIR, parent_name, join(split_index, ".jls"))
    println("Loading meshset for ", parent_name, ": child ", split_index, " / ", count_children(parent_name));
    return load(filename);
end

function has_parent(meshset::MeshSet)
	if !haskey(meshset.properties, "meta") || !haskey(meshset.properties["meta"], "parent") || meshset.properties["meta"]["parent"] == nothing return false end
	return true
end

function get_split_index(meshset::MeshSet)
	if !haskey(meshset.properties, "meta") || !haskey(meshset.properties["meta"], "parent") || meshset.properties["meta"]["parent"] == nothing return 0 end
	return meshset.properties["meta"]["split_index"];
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
  filename = get_filename(firstindex, lastindex)
  println("Loading meshset from ", filename)
  return load(filename)
end

function count_children(parent_name)
    dirname = joinpath(ALIGNED_DIR, parent_name)
    if !isdir(dirname) return 0 end
    dircontents = readdir(dirname);
    return length(find(i -> contains(dircontents[i], ".jls"), 1:length(dircontents)))
end

function load_split(parent_name, split_index)
    filename = joinpath(ALIGNED_DIR, parent_name, "$split_index.jls")
    println("Loading meshset for ", parent_name, ": child ", split_index, " / ", count_children(parent_name));
    return load(filename);
end

function load(index)
  if is_montaged(index)
	  return load_montaged(index[1:2]...)
  elseif is_prealigned(index)
	  return load(montaged(index), (get_preceding(montaged(index)))) 
  end
end

function load(filename::String)
  if !isfile(filename) return nothing end
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

#=
function get_param(meshset::MeshSet, property_name)
  property = nothing;
  if haskey(meshset.properties, "params")
    if haskey(meshset.properties["params"], "match")
      if haskey(meshset.properties["params"]["match"], property_name)
        property = meshset.properties[property_name]
      end
    end
  end
  end  
  return property
end =#

"""
Cycle through index range and return list of flagged meshset indices
"""
function view_flags(firstindex::Index, lastindex::Index)
  flagged_indices = []
  indexB = firstindex
  if is_montaged(indexB)
    indexB = premontaged(indexB)
  end

  for indexA in get_index_range(firstindex, lastindex)
    if is_montaged(indexA)
      indexA = premontaged(indexA)
      indexB = indexA
    end
    meshset = load(indexA, indexB)
    if is_flagged(meshset)
      push!(flagged_indices, (indexA, indexB))
    end
    indexB = indexA
  end

  return flagged_indices
end
