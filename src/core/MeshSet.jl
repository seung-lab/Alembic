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
function to_csv(ms::MeshSet; manual_only=true)
  for match in ms.matches
    to_csv(match, manual_only=manual_only)
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
function count_filtered_correspondences(ms::MeshSet; manual_only=false)	return sum(map(x -> count_filtered_correspondences(x, manual_only=manual_only), ms.matches));		end

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
  fix!(get_mesh(ms, mesh_index)); 	
end

function fix!(ms::MeshSet, mesh_indices)
	for mesh_index in mesh_indices	
    fix!(ms, mesh_index);		
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
  unfix!(get_mesh(ms, mesh_index));	
end

function unfix!(ms::MeshSet, mesh_indices)
	for mesh_index in mesh_indices
  	unfix!(ms, mesh_index);
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

function save_meshes(ms::MeshSet)
  for z in unique(collect_z(ms))
    meshes = get_subsections(ms, z)
    meshset = MeshSet();
    append!(meshset.meshes, meshes)
    save(meshset, :mesh, string(z))
  end
end

function save_matches(ms::MeshSet)
  for match in ms.matches
    save(match, :match)
  end
end

function split_meshset(ms::MeshSet)
  save_meshes(ms)
  save_matches(ms)
end

function compile_meshset()
  ms = compile_meshes()
  compile_matches!(ms)
  return ms
end

function compile_meshset(z_range)
  ms = compile_meshes(z_range)
  compile_matches!(ms)
  return ms
end

function compile_meshset(z_range, pairs)
  ms = compile_meshes(z_range)
  compile_matches!(ms, pairs)
  return ms
end

function compile_meshes(z_range=get_z_range())
  ms = MeshSet()
  files = [string(z) for z in z_range]
  results, empties, errors = load(:mesh, files)
  for r in results
    println("Appending mesh $(r["filename"])")
    m = r["content"]
    append!(ms.meshes, m.meshes)
  end
  return ms
end

function get_pairs(ms::MeshSet)
  sections = []
  for m in ms.meshes
    z = get_z(m)
    if :subsections in keys(m.properties)
      append!(sections, m.properties[:subsections] + z)
    else
      push!(sections, z)
    end
  end
  return get_all_overlaps(sections, PARAMS[:match][:depth], 
                                        symmetric=PARAMS[:match][:symmetric])
end

function compile_matches!(ms::MeshSet, pairs=get_pairs(ms))
  files = [string(p) for p in pairs]
  results, empties, errors = load(:match, files)
  for r in results
    println("Appending match $(r["filename"])")
    m = r["content"]
    push!(ms.matches, m)
  end
end

"""
Merge meshes from list A into list B if they don't already exist
"""
function merge_meshes!(meshesA, meshesB)
  indices = map(get_index, meshesA)
  for mesh in meshesB
    index = get_index(mesh)
    if !(index in indices)
      push!(meshesA, mesh) 
    end
  end
  # sort!(meshesB; by=get_index)
end

"""
Combine unique matches and meshes from ms_two into ms_one
"""
function concat!(msA::MeshSet, msB::MeshSet)
  for match in msB.matches
    src_index = get_src_index(match)
    dst_index = get_dst_index(match)
    ind = find_match_index(msA, src_index, dst_index)
    if ind == 0
      push!(msA.matches, match)
    end
  end
  merge_meshes!(msA.meshes, msB.meshes)
  # sort!(ms_one.matches; by=get_dst_index)
  # sort!(ms_one.matches; by=get_src_index)
	return msA;
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

function mark_solved!(ms::MeshSet)
  ms.properties[:meta][:solved] = author()
end

function reset!(ms::MeshSet, index)
  try
    mi = find_mesh_index(ms, index)
    m = ms.meshes[mi]
    ms.meshes[mi] = deepcopy(m, dst_nodes=m.src_nodes)
  end
end

function reset!(ms::MeshSet)
  map(reset!, Iterators.repeated(ms), collect_mesh_indices(ms))
end

### match
function get_all_overlaps(ms::MeshSet, within=1; symmetric=true)	
  return get_all_overlaps(ms.meshes, within; symmetric=symmetric);	
end

function get_all_overlaps(meshes::Array{Mesh, 1}, within=1; symmetric=true)
  return get_all_overlaps(map(get_index, meshes), within; symmetric=symmetric)
end

function get_all_overlaps(z_range, within=1; symmetric=true)
  println("Finding all overlaps within $within")
  preceding_pairs = Array{Tuple{Number,Number},1}()
  succeeding_pairs = Array{Tuple{Number,Number},1}()

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
function match!(ms::MeshSet, within::Int64=PARAMS[:match][:depth]; 
                                          symmetric=PARAMS[:match][:symmetric])
	pairs = get_all_overlaps(ms, within; symmetric=symmetric);
  match!(ms, pairs)
end

"""
Add meshes for sections with multiple mask components
"""
function add_subsection!(ms::MeshSet, index)
  if !(index in ms) & is_subsection(index) & (get_z(index) in ms)
    println("Adding subsection $index")
    mesh = get_mesh(ms, get_z(index))
    new_mesh = deepcopy(mesh, index=index)
    add_mesh!(ms, new_mesh)
  end
end

"""
Make sure that all meshes required by matches are present
"""
function sync!(ms::MeshSet)
  for m in ms.matches
    src_index, dst_index = get_index(m)
    add_subsection!(ms, src_index)
    add_subsection!(ms, dst_index)
  end
end

function match!(ms::MeshSet, pairs::Array)
	for pair in pairs
    matches = get_matches(get_mesh(ms, pair[1]), get_mesh(ms, pair[2]))
    for m in matches
      println((get_index(m), length(m.src_points)))
      if length(m.src_points) > 0
    		push!(ms.matches, m)
      end
    end
	end
  sync!(ms)
end

"""
Remove unconnected meshes & matches below correspondence threshold. 
Return those meshes & matches in their own meshset.
"""
function prune!(ms::MeshSet; min_match_count=3, manual_only=true)
  meshes_match_count = Dict(get_index(mesh) => 0 for mesh in ms.meshes)
  mesh_indices = keys(meshes_match_count)
  matches_ret = fill(true, count_matches(ms))
  for (i, match) in enumerate(ms.matches)
    match_count = count_filtered_correspondences(match, manual_only=manual_only)
    if (size(match.filters, 1) > 0) & (match_count > min_match_count)
      src_index, dst_index = get_index(match)
      if (src_index in mesh_indices) & (dst_index in mesh_indices)
        match_count = count_filtered_correspondences(match, manual_only=manual_only)
        meshes_match_count[src_index] += match_count
        meshes_match_count[dst_index] += match_count
      end
    end
  end

  for (i, match) in enumerate(ms.matches)
    match_count = count_filtered_correspondences(match, manual_only=manual_only)
    if (size(match.filters, 1) > 0) & (match_count > min_match_count)
      src_index, dst_index = get_index(match)
      if (src_index in mesh_indices) & (dst_index in mesh_indices)
        src_count = meshes_match_count[src_index]
        dst_count = meshes_match_count[dst_index]
        if (src_count < min_match_count) | (dst_count < min_match_count)
          matches_ret[i] = false
        end
      else
        matches_ret[i] = false
      end
    else
      matches_ret[i] = false
    end
  end

  meshes_ret = [meshes_match_count[get_index(m)] >= min_match_count for m in ms.meshes]
  removed_ms = MeshSet()
  removed_ms.meshes = ms.meshes[.~meshes_ret]
  removed_ms.matches = ms.matches[.~matches_ret]
  println("Removed $(length(removed_ms.matches)) / $(length(ms.matches)) matches")
  println("Removed $(length(removed_ms.meshes)) / $(length(ms.meshes)) meshes")
  ms.meshes = ms.meshes[meshes_ret]
  ms.matches = ms.matches[matches_ret]
  return removed_ms
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
