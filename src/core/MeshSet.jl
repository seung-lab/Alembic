type MeshSet
 	meshes::Array{Mesh, 1}			    # vector of meshes in the set
 	matches::Array{Match, 1}		    # vector of matches in the set
  properties::Dict{Symbol, Any}		# in the same array order, contains parameters and filters
end

### IO.jl
function get_index(ms::MeshSet) 	
  return get_index(ms.meshes[1]), get_index(ms.meshes[end]) 
end

function get_images(ms::MeshSet, dtype = UInt8)	return map(get_image, ms.meshes, repeated(dtype));	end
function get_correspondence_patches(ms::MeshSet, match_ind, corr_ind) return get_correspondence_patches(ms.matches[match_ind], corr_ind)	end
function get_params(ms::MeshSet)			return ms.properties[:params];		end

function Base.string(ms::MeshSet)
  s = Dict()
  s["meshes"] = [get_name(m) for m in ms.meshes]
  s["matches"] = [get_name(m) for m in ms.matches]
  s["properties"] = ms.properties
  return JSON.dump(s)
end

# IO
function to_csv(ms::MeshSet)
  for match in ms.matches
    to_csv(match)
  end
end

function load_manual_filters!(ms::MeshSet)
  for match in ms.matches
    load_manual_filters!(match)
  end
end

### counting
function count_meshes(ms::MeshSet)			return length(ms.meshes);		end
function count_matches(ms::MeshSet)		return length(ms.matches);		end
function count_nodes(ms::MeshSet)			return sum(map(count_nodes, ms.meshes));		end
function count_edges(ms::MeshSet)			return sum(map(count_edges, ms.meshes));		end
function count_correspondences(ms::MeshSet)	return sum(map(count_correspondences, ms.matches));		end
function count_filtered_correspondences(ms::MeshSet)	return sum(map(count_filtered_correspondences, ms.matches));		end

### finding
function find_mesh_index(ms::MeshSet, index)  	return findfirst(this -> index == get_index(this), ms.meshes)	end
function find_mesh_index(ms::MeshSet, mesh::Mesh)	return find_mesh_index(ms, get_index(mesh))			end
function find_match_index(ms::MeshSet, src_index, dst_index)
    return findfirst(this -> (src_index == get_src_index(this)) && (dst_index == get_dst_index(this)), ms.matches)
end
function find_match_index(ms::MeshSet, match::Match)
  return find_match_index(ms, get_src_and_dst_indices(match)...);
end
function find_match_indices(ms::MeshSet, index)
  indices = []
  indices = [indices, find(this->index==get_src_index(this), ms.matches)]
  indices = [indices, find(this->index==get_dst_index(this), ms.matches)]
  return indices
end

### getting
function collect_z(ms::MeshSet)
  return map(get_z, ms.meshes)
end

function collect_mesh_indices(ms::MeshSet)
  return map(get_index, ms.meshes)
end

function Base.in(index, ms::MeshSet)
  return index in collect_mesh_indices(ms)
end

function Base.in(mesh::Mesh, ms::MeshSet)
  return get_index(mesh) in collect_mesh_indices(ms)
end

function get_subsections(ms::MeshSet, index)
  return ms.meshes[collect_z(ms) .== get_z(index)]
end

function get_subsections(ms::MeshSet, mesh::Mesh)
  return ms.meshes[collect_z(ms) .== get_z(mesh)]
end

function get_mesh(ms::MeshSet, index)
  return ms.meshes[find_mesh_index(ms, index)]
end

function get_meshes(ms::MeshSet)
  return ms.meshes
end

function get_match(ms::MeshSet, index)
  k = find(i -> (index == get_src_index(i)) || (index == get_dst_index(i)), ms.matches)
  return ms.matches[k]
end

function get_bbox(ms::MeshSet; globalized::Bool = false, use_post::Bool=false, scale=1.0)
  sum(map(x -> get_bbox(x, globalized=globalized, use_post=use_post, scale=scale), ms.meshes))
end

# Mesh.jl extensions
function fix!(ms::MeshSet, mesh_index::Int64) 	
  fix!(ms.meshes[mesh_index]); 	
end

function fix!(ms::MeshSet, mesh_indices)
	for i in mesh_indices	
    fix!(ms.meshes[i]);		
  end
end

function fix!(ms::MeshSet)
	for mesh in ms.meshes			
    fix!(mesh);			 	
  end
end

function fix_ends!(ms::MeshSet)
   unfix!(ms); fix!(ms, 1); fix!(ms, length(ms.meshes));
end

function unfix!(ms::MeshSet, mesh_index::Int64)	
  unfix!(ms.meshes[mesh_index]);	
end

function unfix!(ms::MeshSet, mesh_indices)
	for mesh_index in mesh_indices
  	unfix!(ms.meshes[mesh_index]);
  end
end

function unfix!(ms::MeshSet)
  for mesh in ms.meshes
  	unfix!(mesh);
  end
end

### adding
function add_mesh!(ms::MeshSet, mesh::Mesh)
  push!(ms.meshes, mesh);
  sort!(ms.meshes; by=get_index)
end

function add_match!(ms::MeshSet, match::Match)
  push!(ms.matches, match);
  sort!(ms.matches; by=get_src_and_dst_indices)
end

function add_matches!(ms::MeshSet, matches::Array{Match,1})
  append!(ms.matches, matches);
  sort!(ms.matches; by=get_src_and_dst_indices)
end

function remove_match!(ms::MeshSet, src_index, dst_index)
  i = find_match_index(ms, src_index, dst_index)
  deleteat!(ms.matches, i)
end

function remove_mesh!(ms::MeshSet, index)
  i = find_mesh_index(ms, index)
  if i != 0
    deleteat!(ms.meshes, i)
    match_indices = find_match_indices(ms, index)
    mask = collect(setdiff(Set(1:length(ms.matches)), match_indices))
    ms.matches = ms.matches[mask]
  end
end

### reviewing
function set_reviewed!(ms::MeshSet, match_ind, flag = false)
	set_reviewed!(ms.matches[match_ind], flag);
end

function is_reviewed(ms::MeshSet, match_ind)
  return is_reviewed(ms.matches[match_ind])
end

function is_reviewed(ms::MeshSet)
  return (&)(map(is_reviewed, ms.matches)...)
end

function count_flags(ms::MeshSet)
  return sum(map(is_flagged, ms.matches))
end

function flag!(ms::MeshSet, match_ind)
  flag!(ms.matches[match_ind])
end

function unflag!(ms::MeshSet, match_ind)
  unflag!(ms.matches[match_ind])
end

function unflag!(ms::MeshSet)
  map(unflag!, ms.matches)
end

function is_flagged(ms::MeshSet, match_ind)
  return is_flagged(ms.matches[match_ind])
end

function is_flagged(ms::MeshSet)
  return |(map(is_flagged, ms.matches)...)
end

function check!(ms::MeshSet, crits=ms.properties[:params][:review])
  unflag!(ms)
  ms.properties[:params][:review] = crits
  return |(map(check!, ms.matches, repeated(Base.values(crits)))...)
end

function Base.filter!(ms::MeshSet, filters::Dict=ms.properties[:params][:filter])
  filters = collect(Base.values(filters))
  filters = filters[sortperm(map(getindex, filters, repeated(1)))]
  for filter in filters
    filter!(ms, filter)
  end
  passed = !check!(ms)
  passed ? println("check passed") : println("check failed")
end

function Base.filter!(ms::MeshSet, filter::Tuple)
  total = sum(map(filter!, ms.matches, repeated(filter)))
  println("$total / $(count_correspondences(ms)) correspondences filtered on $filter")
end

function clear_filters!(ms::MeshSet)
  for match in ms.matches
    clear_filters!(match)
  end
end

function refilter!(ms::MeshSet, filters=ms.properties[:params][:filter])
  clear_filters!(ms)
  filter!(ms, filters)
end

function check_and_fix!(ms::MeshSet, 
                  crits = Base.values(ms.properties[:params][:review]), 
                  filters = Base.values(ms.properties[:params][:filter])) 
  if check!(ms, crits)
    for match in ms.matches
      if is_flagged(match)
      	clear_filters!(match);
      	map(filter!, repeated(match), filters)
      end
    end
    fixed = !check!(ms, crits);
    if fixed
      println("fixed successfully");
    else 
      println("failed to fix ms")
      for match in ms.matches
        if is_flagged(match)
        	# clear_filters!(match);
          flag!(match)
        end
      end
    end
  end
end

function get_parent(ms::MeshSet)
  parent = nothing
  if haskey(ms.properties, :meta)
    if haskey(ms.properties[:meta], :parent)
      parent = ms.properties[:meta][:parent]
    end
  end
  return parent
end

function split_ms(ms::MeshSet)
	for mesh in ms.meshes
	  save(mesh);
	end
	for match in ms.matches
	  save(match)
	end
	return (vcat(map(get_path, ms.meshes), map(get_path, ms.matches)));
end

function compile_ms(first_index, last_index)
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
  ms.properties = Dict{Symbol, Any}(:author => author(), :meta => Dict{Symbol, Any}(), :params => deepcopy(meshes[1].properties[:params]))
  
  @everywhere gc();

  return ms;
end

"""
Create partial ms of a ms split into children
"""
function compile_partial_ms(parent_name, firstindex, lastindex)
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
    if index in ms
      push!(ind, i)
    end
  end
  return ind
end

function concat_mss(filenames::Array)
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
Combine unique matches and meshes from ms_two into ms_one
"""
function concat!(ms_one::MeshSet, ms_two::MeshSet)
  for match in ms_two.matches
    src_index = get_src_index(match)
    dst_index = get_dst_index(match)
    ind = find_match_index(ms_one, src_index, dst_index)
    if ind == 0
      push!(ms_one.matches, match)
    end
  end
  merge_meshes!(ms_two.meshes, ms_one.meshes)
  sort!(ms_one.matches; by=get_dst_index)
  sort!(ms_one.matches; by=get_src_index)
	return ms_one;
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
          )
	  )
  ms = MeshSet(meshes, matches, properties);
  match!(ms, PARAMS[:match][:depth]; symmetric = PARAMS[:match][:symmetric]);

  # save(ms)
  filter!(ms);
  check!(ms);
  # save(ms);

  if solve == true
    solve!(ms, method=solve_method);
    # save(ms);
  end

  return ms;
end

function mark_solved!(ms::MeshSet)
  ms.properties[:meta][:solved] = author()
end

function reset!(ms::MeshSet)
  m = string("Are you sure you want to reset all the meshes in ", get_name(ms), "?")
  if user_approves(m)
    for mesh in ms.meshes
      reset!(mesh)
    end
    ms.properties[:meta][:reset] = author()
  end
end

### match
function get_all_overlaps(ms::MeshSet, within = 1; symmetric = true)	
  return get_all_overlaps(ms.meshes, within; symmetric = symmetric);	
end

function get_all_overlaps(meshes::Array{Mesh, 1}, within = 1; symmetric = true)
  return get_all_overlaps(1:length(meshes), within=within; symmetric=symmetric)
end

function get_all_overlaps(z_range, within = 1; symmetric = true)
  preceding_pairs = Pairings(0)
  succeeding_pairs = Pairings(0)

  for i in z_range, j in z_range
    if is_preceding(i, j, within) 
      push!(preceding_pairs, (i, j)); 
      if symmetric 
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
function match!(ms::MeshSet, within=1; symmetric=PARAMS[:match][:symmetric])
	pairs = get_all_overlaps(ms, within; symmetric=symmetric);
	for pair in pairs
		add_matches!(ms, get_matches(ms.meshes[pair[1]], ms.meshes[pair[2]]));
	end

  # add meshes for sections with multiple mask components
  function add_subsection!(ms, index)
    if !(index in ms) & is_subsection(index) & (get_z(index) in ms)
      mesh = ms.meshes[find_mesh_index(ms, get_z(index))]
      new_mesh = deepcopy(mesh, index=index)
      add_mesh!(ms, new_mesh)
    end
  end

  for m in ms.matches
    add_subsection!(ms, m.src_index)
    add_subsection!(ms, m.dst_index)
  end
end

# filters any correspondence that cannot be triangulated
function sanitize!(ms::MeshSet)
  meshes = Dict{Any, Any}();
  for mesh in ms.meshes
  	meshes[get_index(mesh)] = mesh;
  end

  for match in ms.matches
  	src_mesh = meshes[match.src_index];
  	dst_mesh = meshes[match.dst_index];
  	src_pts, dst_pts = get_correspondences(match);
  	src_pt_triangles = map(find_mesh_triangle, repeated(src_mesh), columnviews(src_pts));
  	dst_pt_triangles = map(find_mesh_triangle, repeated(dst_mesh), columnviews(dst_pts));
  	invalids = union(find(ind -> src_pt_triangles[ind] == NO_TRIANGLE, 1:count_correspondences(match)), find(ind -> dst_pt_triangles[ind] == NO_TRIANGLE, 1:count_correspondences(match)))
  	if length(invalids) !=0
    	clear_filters!(match; filtertype = :sanitization);
    	# filter_manual!(match, invalids; filtertype = :sanitization);
  	end
  end
end

function parse_ms_filename(name::AbstractString)
    indexA, indexB = NO_INDEX, NO_INDEX
    # aligned-section
    m = Base.match(r"(\d*),(\d*)-(\d*),(\d*)_aligned", name)
    if typeof(m) != Void
    indexA = aligned(parse(Int, m[1]), parse(Int, m[2]))
    indexB = aligned(parse(Int, m[3]), parse(Int, m[4]))
    end
    return indexA, indexB
end

function make_stack(start_index, end_index)
	return make_stack(intersect(start_index:end_index, get_z_range()))
end

function make_stack(indices = get_z_range())
	ms = MeshSet();
  for i in indices
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
  ms = MeshSet(PARAMS_ALIGNMENT)
  images = Array{SharedArray{UInt8, 2}, 1}(0)

  for i in indices
  name = offsets[i, 1]
  index = offsets[i, 2]
  dx = offsets[i, 4]
  dy = offsets[i, 3]
  image = get_image(get_path(name))
  add_mesh(Mesh(name, image, index, dy, dx, false, PARAMS_ALIGNMENT), ms)
  
  image_shared = SharedArray{UInt8}(size(image, 1), size(image, 2))
  image_shared[:, :] = image[:, :]
  push!(images, image_shared)
  end

  optimize_all_cores(ms.params)

  return ms, images
end

function get_correspondence_properties(ms::MeshSet, key)
  return map(get_correspondence_properties, ms.matches, repeated(key))
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
