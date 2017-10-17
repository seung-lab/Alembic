type MeshSet
 	meshes::Array{Mesh, 1}			    # vector of meshes in the set
 	matches::Array{Match, 1}		    # vector of matches in the set
  properties::Dict{Symbol, Any}		# in the same array order, contains parameters and filters
end

### IO.jl
function get_index(meshset::MeshSet) 	
  return get_index(meshset.meshes[1]), get_index(meshset.meshes[end]) 
end

function get_images(meshset::MeshSet, dtype = UInt8)	return map(get_image, meshset.meshes, repeated(dtype));	end
function get_correspondence_patches(meshset::MeshSet, match_ind, corr_ind) return get_correspondence_patches(meshset.matches[match_ind], corr_ind)	end
function get_params(meshset::MeshSet)			return meshset.properties[:params];		end

function Base.string(ms::MeshSet)
  s = Dict()
  s["meshes"] = [get_name(m) for m in ms.meshes]
  s["matches"] = [get_name(m) for m in ms.matches]
  s["properties"] = ms.properties
  return JSON.dump(s)
end

#function get_fixed(meshset::MeshSet)			return meshset.properties["fixed"];		end

### counting
function count_meshes(meshset::MeshSet)			return length(meshset.meshes);		end
function count_matches(meshset::MeshSet)		return length(meshset.matches);		end
function count_nodes(meshset::MeshSet)			return sum(map(count_nodes, meshset.meshes));		end
function count_edges(meshset::MeshSet)			return sum(map(count_edges, meshset.meshes));		end
function count_correspondences(meshset::MeshSet)	return sum(map(count_correspondences, meshset.matches));		end
function count_filtered_correspondences(meshset::MeshSet)	return sum(map(count_filtered_correspondences, meshset.matches));		end

function get_mesh(meshset::MeshSet, index)
  return meshset.meshes[find_mesh_index(meshset, index)]
end

function get_matches(meshset::MeshSet, index)
  k = find(i -> (index == get_src_index(i)) || (index == get_dst_index(i)), meshset.matches)
  return meshset.matches[k]
end

### finding
function find_mesh_index(meshset::MeshSet, index)  	return findfirst(this -> index == get_index(this), meshset.meshes)	end
function find_mesh_index(meshset::MeshSet, mesh::Mesh)	return find_mesh_index(meshset, get_index(mesh))			end
function contains_mesh(meshset::MeshSet, index)		return !(find_mesh_index(meshset, index) == 0)				end
function contains_mesh(meshset::MeshSet, mesh::Mesh)	return contains_mesh(meshset, get_index(mesh))				end
function find_match_index(meshset::MeshSet, src_index, dst_index)
  	return findfirst(this -> (src_index == get_src_index(this)) && (dst_index == get_dst_index(this)), meshset.matches)
end
function find_match_index(meshset::MeshSet, match::Match)
	return find_match_index(meshset, get_src_and_dst_indices(match)...);
end
function find_match_indices(meshset::MeshSet, index)
  indices = []
  indices = [indices, find(this->index==get_src_index(this), meshset.matches)]
  indices = [indices, find(this->index==get_dst_index(this), meshset.matches)]
  return indices
end

### finding subs for purposes of making submeshsets
function find_matches_subarray(meshset::MeshSet, meshes::Array{Mesh, 1})
	matches_sub = Array{Match, 1}(0);
	meshes_inds = [get_index(mesh) for mesh in meshes]
	for match in meshset.matches
	  if issubset(get_src_and_dst_indices(match), meshes_inds)
	    push!(matches_sub, match);
	  end
	end
	return matches_sub;
end

function find_meshes_subarrays(meshset::MeshSet)
  meshes_subs = Array{Array{Mesh, 1}, 1}(0);
  current_sub = Array{Mesh, 1}(0)
  for (index, mesh) in enumerate(meshset.meshes)
      push!(current_sub, mesh);
      if (is_fixed(mesh) && index != 1) || index == count_meshes(meshset)
	if !is_fixed(current_sub[end-1])
		push!(meshes_subs, current_sub);
        end
  	current_sub = Array{Mesh, 1}(0)
        push!(current_sub, mesh);
      end
  end
  return meshes_subs;
end

function make_submeshsets(meshset::MeshSet)
   meshes_subs = find_meshes_subarrays(meshset);
   matches_subs = [find_matches_subarray(meshset, meshes) for meshes in meshes_subs];
   meshsets = Array{MeshSet, 1}(map(MeshSet, meshes_subs, matches_subs, repeated(meshset.properties)));
   return meshsets;
end

# Mesh.jl extensions
function fix!(meshset::MeshSet, mesh_ind::Int64) 	
  fix!(meshset.meshes[mesh_ind]); 	
end

function fix!(meshset::MeshSet, mesh_ind_first::Int64, mesh_ind_last::Int64)
	for mesh_ind in mesh_ind_first:mesh_ind_last	
    fix!(meshset.meshes[mesh_ind]);		
  end
end

function fix!(meshset::MeshSet)
	for mesh in meshset.meshes			
    fix!(mesh);			 	
  end
end

function fix_ends!(meshset::MeshSet)
   unfix!(meshset); fix!(meshset, 1); fix!(meshset, length(meshset.meshes));
end

function unfix!(meshset::MeshSet, mesh_ind::Int64)	
  unfix!(meshset.meshes[mesh_ind]);	
end

function unfix!(meshset::MeshSet, mesh_ind_first::Int64, mesh_ind_last::Int64)
	for mesh_ind in mesh_ind_first:mesh_ind_last
  	unfix!(meshset.meshes[mesh_ind]);
  end
end

function unfix!(meshset::MeshSet)
  for mesh in meshset.meshes
  	unfix!(mesh);
  end
end

### adding
function add_mesh!(meshset::MeshSet, mesh::Mesh)
  push!(meshset.meshes, mesh);
  sort!(meshset.meshes; by=get_index)
end

function add_match!(meshset::MeshSet, match::Match)
  push!(meshset.matches, match);
  sort!(meshset.matches; by=get_src_and_dst_indices)
end

function add_matches!(meshset::MeshSet, matches::Array{Match,1})
  append!(meshset.matches, matches);
  sort!(meshset.matches; by=get_src_and_dst_indices)
end

function remove_match!(meshset::MeshSet, src_index, dst_index)
  i = find_match_index(meshset, src_index, dst_index)
  deleteat!(meshset.matches, i)
end

function remove_mesh!(meshset::MeshSet, index)
  i = find_mesh_index(meshset, index)
  if i != 0
    deleteat!(meshset.meshes, i)
    match_indices = find_match_indices(meshset, index)
    mask = collect(setdiff(Set(1:length(meshset.matches)), match_indices))
    meshset.matches = meshset.matches[mask]
  end
end

### reviewing
function set_reviewed!(meshset::MeshSet, match_ind, flag = false)
	set_reviewed!(meshset.matches[match_ind], flag);
end

function is_reviewed(meshset::MeshSet, match_ind)
  return is_reviewed(meshset.matches[match_ind])
end

function is_reviewed(meshset::MeshSet)
  return (&)(map(is_reviewed, meshset.matches)...)
end

function count_flags(meshset::MeshSet)
  return sum(map(is_flagged, meshset.matches))
end

function flag!(meshset::MeshSet, match_ind)
  flag!(meshset.matches[match_ind])
end

function unflag!(meshset::MeshSet, match_ind)
  unflag!(meshset.matches[match_ind])
end

function unflag!(meshset::MeshSet)
  map(unflag!, meshset.matches)
end

function is_flagged(meshset::MeshSet, match_ind)
  return is_flagged(meshset.matches[match_ind])
end

function is_flagged(meshset::MeshSet)
  return |(map(is_flagged, meshset.matches)...)
end

function check!(meshset::MeshSet, crits=meshset.properties[:params][:review])
  unflag!(meshset)
  meshset.properties[:params][:review] = crits
  return |(map(check!, meshset.matches, repeated(Base.values(crits)))...)
end

function filter!(meshset::MeshSet, filters::Dict=meshset.properties[:params][:filter])
  filters = collect(Base.values(filters))
  filters = filters[sortperm(map(getindex, filters, repeated(1)))]
  for filter in filters
    filter!(meshset, filter)
  end
  passed = !check!(meshset)
  passed ? println("check passed") : println("check failed")
end

function filter!(meshset::MeshSet, filter::Tuple)
  total = sum(map(filter!, meshset.matches, repeated(filter)))
  println("$total / $(count_correspondences(meshset)) correspondences filtered on $filter")
end

function clear_filters!(meshset::MeshSet)
  for match in meshset.matches
    clear_filters!(match)
  end
end

function refilter!(meshset::MeshSet, filters=meshset.properties[:params][:filter])
  clear_filters!(meshset)
  filter!(meshset, filters)
end

function check_and_fix!(meshset::MeshSet, 
                  crits = Base.values(meshset.properties[:params][:review]), 
                  filters = Base.values(meshset.properties[:params][:filter])) 
  if check!(meshset, crits)
    for match in meshset.matches
      if is_flagged(match)
      	clear_filters!(match);
      	map(filter!, repeated(match), filters)
      end
    end
    fixed = !check!(meshset, crits);
    if fixed
      println("fixed successfully");
    else 
      println("failed to fix meshset")
      for match in meshset.matches
        if is_flagged(match)
        	# clear_filters!(match);
          flag!(match)
        end
      end
    end
  end
end

function get_parent(meshset::MeshSet)
  parent = nothing
  if haskey(meshset.properties, :meta)
    if haskey(meshset.properties[:meta], :parent)
      parent = meshset.properties[:meta][:parent]
    end
  end
  return parent
end

function split_meshset(meshset::MeshSet)
	for mesh in meshset.meshes
	  save(mesh);
	end
	for match in meshset.matches
	  save(match)
	end
	return (vcat(map(get_path, meshset.meshes), map(get_path, meshset.matches)));
end

function compile_meshset(first_index, last_index)
  inds = get_index_range(first_index, last_index);
  println("Compiling MeshSet between $first_index and $last_index - $(length(inds)) meshes expected")
  mesh_inds = Array{Index, 1}();
  match_inds = Array{Tuple{Index, Index}, 1}();

  for ind in inds
    if isfile(get_path(Mesh, ind)) push!(mesh_inds, ind) end
  end
  for ind in inds
    for dst_ind in inds
    if isfile(get_path(Match, (ind, dst_ind))) push!(match_inds, (ind, dst_ind)) end
  end
  end

  meshes = pmap(load, repeated(Mesh), mesh_inds);
  matches = pmap(load, repeated(Match), match_inds);
  #=
#  @sync begin
  for ind in inds
    if isfile(get_path(Mesh, ind)) @async push!(meshes, remotecall_fetch(WORKER_PROCS[ind[2]%length(WORKER_PROCS) + 1], load, Mesh, ind)); end
  end
  for ind in inds
    for dst_ind in inds
      if isfile(get_path(Match, (ind, dst_ind))) @async push!(matches, remotecall_fetch(WORKER_PROCS[ind[2]%length(WORKER_PROCS) + 1], load, Match, (ind, dst_ind))) end
    end
  end
#end #sync=#

  println("$(length(meshes)) meshes found...")
  if length(meshes) == 0 println("no meshes exist between the requested indices - aborting...."); return nothing end
  if length(meshes) != length(inds) println("not all meshes exist between requested indices - continuing anyway...") end
  if |(map(is_fixed, meshes)...) println("there are fixed meshes among the requested meshes - continuing anyway...") end
  println("$(length(matches)) matches found...")

  if length(matches) == 0 println("no matches exist between requested indices - continuing anyway...") end

  sort!(meshes; by=get_index)
  sort!(matches; by=get_dst_index)
  sort!(matches; by=get_src_index)

  ms = MeshSet();
  ms.meshes = meshes;
  ms.matches = matches;
  ms.properties = Dict{Symbol, Any}(:author => author(), :meta => Dict{Symbol, Any}(:parent => nothing, :split_index => 0), :params => deepcopy(meshes[1].properties[:params]))
  
  @everywhere gc();

  return ms;
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
			ms.properties[:meta][:parent] = nothing;
			ms.properties[:meta][:split_index] = 0;
		end
		println("Child ", i, " / ", count_children(parent_name), " concatanated");
	end
  sort!(ms.meshes; by=get_index)
	return ms;
end

"""
Create partial meshset of a meshset split into children
"""
function compile_partial_meshset(parent_name, firstindex, lastindex)
  ms = MeshSet()
  indices = get_index_range(prealigned(firstindex), prealigned(lastindex))
  ind = []
  ms.properties = deepcopy(load_split(parent_name, 1).properties)
  delete!(ms.properties, :meta)
  for i = 1:count_children(parent_name)
    child_ms = load_split(parent_name, i)
    if length(intersect(map(get_index, child_ms.meshes), indices)) > 0
      concat!(ms, child_ms)
    end
  end
  return ms
end

function find_mesh_indices(firstindex, lastindex, i::Int64)
  indices = get_index_range(firstindex, lastindex)
  return find_mesh_indices(firstindex, lastindex, indices[i])
end

function find_mesh_indices(firstindex, lastindex, index)
  parent_name = get_name(firstindex, lastindex)
  ind = []
  for i in 1:count_children(parent_name)
    ms = load_split(parent_name, i)
    if contains_mesh(ms, index)
      push!(ind, i)
    end
  end
  return ind
end

function concat_meshsets(filenames::Array)
  println(filenames[1])
  ms = load(filenames[1])
  for fn in filenames[2:end]
    println(fn)
    ms2 = load(fn)
    concat!(ms, ms2)
  end
  save(ms)
  return ms
end

"""
Merge meshes from list A into list B if they don't already exist
"""
function merge_meshes!(meshesA, meshesB)
  indices = map(get_index, meshesB)
  for mesh in meshesA
    index = get_index(mesh)
    if !(index in indices)
      push!(meshesB, mesh) 
    end
  end
  sort!(meshesB; by=get_index)
end

"""
Combine unique matches and meshes from meshset_two into meshset_one
"""
function concat!(meshset_one::MeshSet, meshset_two::MeshSet)
  for match in meshset_two.matches
    src_index = get_src_index(match)
    dst_index = get_dst_index(match)
    ind = find_match_index(meshset_one, src_index, dst_index)
    if ind == 0
      push!(meshset_one.matches, match)
    end
  end
  merge_meshes!(meshset_two.meshes, meshset_one.meshes)
  sort!(meshset_one.matches; by=get_dst_index)
  sort!(meshset_one.matches; by=get_src_index)
	return meshset_one;
end

### initialise
function MeshSet()
 	meshes = Array{Mesh, 1}(0)
 	matches = Array{Match, 1}(0)		
  properties = Dict{Symbol, Any}(  
    :params  => PARAMS,
        :author => author(),
        :meta => Dict{Symbol, Any}()
        )
	return MeshSet(meshes, matches, properties)
end

function MeshSet(indices::Array; solve=true, solve_method=:elastic)
  meshes = map(Mesh, indices)
  sort!(meshes; by=get_index)
  matches = Array{Match, 1}(0)    
  properties = Dict{Symbol, Any}(  
  	  :params  => PARAMS,
          :author => author(),
          :meta => Dict{Symbol, Any}(
          :parent => nothing,
          :split_index => 0)
          )
  meshset = MeshSet(meshes, matches, properties);
  match!(meshset, PARAMS[:match][:depth]; reflexive = PARAMS[:match][:reflexive]);

  # save(meshset)
  filter!(meshset);
  check!(meshset);
  # save(meshset);

  if solve == true
    solve!(meshset, method=solve_method);
    # save(meshset);
  end

  return meshset;
end

function mark_solved!(meshset::MeshSet)
  meshset.properties[:meta]["solved"] = author()
end

function reset!(meshset::MeshSet)
  m = string("Are you sure you want to reset all the meshes in ", get_name(meshset), "?")
  if user_approves(m)
    for mesh in meshset.meshes
      reset!(mesh)
    end
    meshset.properties[:meta]["reset"] = author()
  end
end

### match
function get_all_overlaps(meshset::MeshSet, within = 1; reflexive = reflexive)	
  return get_all_overlaps(meshset.meshes, within; reflexive = reflexive);	
end

function get_all_overlaps(meshes::Array{Mesh, 1}, within = 1; reflexive = true)
  preceding_pairs = Pairings(0)
  succeeding_pairs = Pairings(0)

  for i in 1:length(meshes), j in 1:length(meshes)
    if is_preceding(get_index(meshes[i]), get_index(meshes[j]), within) 
    	push!(preceding_pairs, (i, j)); 
    	if reflexive 
        push!(succeeding_pairs, (j, i)); 
      end
    end
  end

  pairs = vcat(preceding_pairs, succeeding_pairs)
  sort!(pairs)
  pairs = unique(pairs)
  println("$(length(pairs)) pairs found")
  return pairs
end

"""
Generate matches between all meshes (& components) that overlap

  See `get_all_overlaps` to see how meshes are determined to overlap.
  Uses `get_matches` which matches taking crack & fold masks into account.
    Meshes will be duplicated so that each mask component has its own mesh.
"""
function match!(ms::MeshSet, within=1; reflexive=true)
	pairs = get_all_overlaps(ms, within; reflexive=reflexive);
	for pair in pairs
		add_matches!(ms, get_matches(ms.meshes[pair[1]], ms.meshes[pair[2]]));
	end
  # add meshes for sections with multiple mask components
  for m in ms.matches
    if is_subsection(m.src_index) & contains_mesh(ms, m.src_index)
      mesh = ms.meshes[find_mesh_index(ms, get_z(m.src_index))]
      new_mesh = deepcopy(mesh, index=m.src_index)
      add_mesh!(ms, new_mesh)
    end
    if is_subsection(m.dst_index) & contains_mesh(ms, m.dst_index)
      mesh = ms.meshes[find_mesh_index(ms, get_z(m.dst_index))]
      new_mesh = deepcopy(mesh, index=m.dst_index)
      add_mesh!(ms, new_mesh)
    end
  end
end

# filters any correspondence that cannot be triangulated
function sanitize!(meshset::MeshSet)
  meshes = Dict{Any, Any}();
  for mesh in meshset.meshes
  	meshes[get_index(mesh)] = mesh;
  end

  for match in meshset.matches
    	src_mesh = meshes[match.src_index];
    	dst_mesh = meshes[match.dst_index];
	src_pts, dst_pts = get_correspondences(match);
	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), columnviews(src_pts));
	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), columnviews(dst_pts));
	invalids = union(find(ind -> src_pt_triangles[ind] == NO_TRIANGLE, 1:count_correspondences(match)), find(ind -> dst_pt_triangles[ind] == NO_TRIANGLE, 1:count_correspondences(match)))
	if length(invalids) !=0
	clear_filters!(match; filtertype = :sanitization);
	filter_manual!(match, invalids; filtertype = :sanitization);
	end
  end
end

function parse_meshset_filename(name::AbstractString)
    indexA, indexB = NO_INDEX, NO_INDEX
    # aligned-section
    m = Base.match(r"(\d*),(\d*)-(\d*),(\d*)_aligned", name)
    if typeof(m) != Void
    indexA = aligned(parse(Int, m[1]), parse(Int, m[2]))
    indexB = aligned(parse(Int, m[3]), parse(Int, m[4]))
    end
    return indexA, indexB
end

function make_stack()
	ms = MeshSet();
  for i in get_z_range()
    add_mesh!(ms, Mesh(i))
  end
  return ms
end

function save_stack(ms::MeshSet)
  for m in ms.meshes
    save(get_name(m), m)
  end
  for m in ms.matches
    save(get_name(m), m)
  end

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
  
  image_shared = SharedArray{UInt8}(size(image, 1), size(image, 2))
  image_shared[:, :] = image[:, :]
  push!(images, image_shared)
  end

  optimize_all_cores(Ms.params)

  return Ms, images
end

function get_correspondence_properties(meshset::MeshSet, key)
  return map(get_correspondence_properties, meshset.matches, repeated(key))
end

function get_correspondence_stats(ms::MeshSet, key)
  s = get_correspondence_properties(ms, key)
  indices = map(get_src_index, ms.matches)
  l = map(length, s)
  mn = map(round, map(minimum, s), repeated(2))
  m = map(round, map(median, s), repeated(2))
  mx = map(round, map(maximum, s), repeated(2))
  sd = map(round, map(std, s), repeated(2))
  return Dict("index"=>indices, "count"=>l, "min"=>mn, "med"=>m, "max"=>mx, "std"=>sd)
end

function compile_correspondence_stats(firstindex, lastindex, key)
  s = []
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = load(premontaged(index), premontaged(index))
    push!(s, get_correspondence_stats(ms, key))
  end
  return s
end

function flatten_correspondence_stats(stats)
  stats_keys = keys(stats[1])
  s = Dict()
  for k in stats_keys
    f = []
    for d in stats
      push!(f, d[k])
    end
    s[k] = vcat(f...)
    s[k] = convert(Array{typeof(stats[1][k][1]), 1}, s[k])
  end
  return s
end

function flag_correspondence_stats(stats, k="med", sigma=1.5)
  s = flatten_correspondence_stats(stats)
  threshold = median(s[k])
  e = round(std(s[k]), 2)
  println("Flagging correspondence stats for $k outside $threshold +/- $sigma*$e")
  flagged = []
  for (i, m) in enumerate(s[k])
    if m > (threshold+sigma*e) || m < (threshold-sigma*e)
      push!(flagged, [s["index"][i], m])
    end
  end
  return flagged
end
