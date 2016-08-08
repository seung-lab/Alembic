# functions to get indices in different pipeline stages
# index -> index
function overview(index)		return (index[1], index[2], OVERVIEW_INDEX, OVERVIEW_INDEX);		end
function premontaged(index)		return (index[1], index[2], PREMONTAGED_INDEX, PREMONTAGED_INDEX);	end
function montaged(index)		return (index[1], index[2], MONTAGED_INDEX, MONTAGED_INDEX);		end
function prealigned(index)		return (index[1], index[2], PREALIGNED_INDEX, PREALIGNED_INDEX);	end
function subsection(index, num)		return (index[1], index[2], PREALIGNED_INDEX, num);			end
function aligned(index)			return (index[1], index[2], ALIGNED_INDEX, ALIGNED_INDEX);		end
function finished(index)		return (index[1], index[2], FINISHED_INDEX, FINISHED_INDEX);		end
# wafer, sec -> index
function overview(wafer, section)	return (wafer, section, OVERVIEW_INDEX, OVERVIEW_INDEX);		end
function premontaged(wafer, section) 	return (wafer, section, PREMONTAGED_INDEX, PREMONTAGED_INDEX);		end
function montaged(wafer, section)	return (wafer, section, MONTAGED_INDEX, MONTAGED_INDEX);		end
function prealigned(wafer, section)	return (wafer, section, PREALIGNED_INDEX, PREALIGNED_INDEX);		end
function subsection(wafer, section, num) return (wafer, section, PREALIGNED_INDEX, num);			end
function aligned(wafer, section)	return (wafer, section, ALIGNED_INDEX, ALIGNED_INDEX);			end
function finished(wafer, section)	return (wafer, section, FINISHED_INDEX, FINISHED_INDEX);		end

# functions for checking whether an index belongs to a particular pipeline stage
function is_overview(index)	 	return index[3:4] == (OVERVIEW_INDEX, OVERVIEW_INDEX);			end
function is_premontaged(index)	 	return index[3] > 0 && index[4] > 0;					end
function is_montaged(index)	 	return index[3:4] == (MONTAGED_INDEX, MONTAGED_INDEX);			end
function is_prealigned(index)	 	return index[3] == PREALIGNED_INDEX;					end
function is_subsection(index)		return index[3] == PREALIGNED_INDEX && index[4] >= 0; 			end
function is_aligned(index)	 	return index[3:4] == (ALIGNED_INDEX, ALIGNED_INDEX);			end
function is_finished(index)	 	return index[3:4] == (FINISHED_INDEX, FINISHED_INDEX);			end

# function to get the next stage
function nextstage(index)
    if is_overview(index)		return premontaged(index)
    elseif is_premontaged(index) 	return montaged(index)
    elseif is_montaged(index) 		return prealigned(index)
    elseif is_prealigned(index) 	return aligned(index)
    elseif is_aligned(index) 		return finished(index)
    elseif is_finished(index) 		return NO_INDEX		end
end

function prevstage(index)
    if is_overview(index)		return NO_INDEX
    elseif is_premontaged(index) 	return overview(index)
    elseif is_montaged(index) 		return premontaged(index)
    elseif is_prealigned(index) 	return montaged(index)
    elseif is_aligned(index) 		return prealigned(index)
    elseif is_finished(index) 		return aligned(index)		end
end

# functions for checking whether two indices are next to each other
function is_adjacent(A_index, B_index)
  if !(is_premontaged(A_index)) || !(is_premontaged(B_index))	return false end
  if abs(A_index[3] - B_index[3]) + abs(A_index[4] - B_index[4])  == 1 return true; end
  return false;
end

function is_diagonal(A_index, B_index)
  if !(is_premontaged(A_index)) || !(is_premontaged(B_index))	return false end
  if abs(A_index[3] - B_index[3]) + abs(A_index[4] - B_index[4]) == 2 && A_index[3] != B_index[3] && A_index[4] != B_index[4] return true; end
  return false;
end

function is_preceding(A_index, B_index, within = 1)
  if (is_premontaged(A_index)) || (is_premontaged(B_index))	return false end
  for i in 1:within if A_index[1:2] == get_preceding(B_index, i)[1:2] return true; end end
  return false;
end

function get_neighbor(index::Index, direction=[0,1])
  if !(is_premontaged(index)) return NO_INDEX end
  neighbor_index = (index[1:2]..., index[3]+direction[1], index[4]+direction[2])
  if !in_registry(neighbor_index) return NO_INDEX end
  return neighbor_index
end

function get_above(index::Index)
  return get_neighbor(index, [-1,0])
end

function get_below(index::Index)
  return get_neighbor(index, [1,0])
end

function get_right(index::Index)
  return get_neighbor(index, [0,1])
end

function get_left(index::Index)
  return get_neighbor(index, [0,-1])
end

function get_cardinal_neighbors(index::Index)
  neighbors = [get_above(index), get_left(index), get_below(index), get_right(index)] 
  return filter(i->i!=NO_INDEX, neighbors)
end

function index_rank(index)
  return index[1]*10^3 + index[2]
end

function index_to_int(index)
  return index[1]*10^7 + index[2]*10^4 + index[3]*10^2 + index[4] 
end

# function isless(indexA, indexB)
#   return Base.isless(index_rank(indexA), index_rank(indexB))
# end

# function isequal(indexA, indexB)
#   return Base.isequal(index_rank(indexA), index_rank(indexB))
# end

function get_overview_index(index)
  return (index[1:2]..., OVERVIEW_INDEX, OVERVIEW_INDEX)
end

function get_registry(index)
  if is_premontaged(index) registry = REGISTRY_PREMONTAGED;
  elseif is_montaged(index) registry = REGISTRY_MONTAGED;
  elseif is_prealigned(index) registry = REGISTRY_PREALIGNED;
  elseif is_aligned(index) registry = REGISTRY_ALIGNED;
  elseif is_finished(index) registry = REGISTRY_FINISHED;
  else registry = Void; println("Index $index does not correspond to a pipeline stage."); end
  return registry; 
end

function get_indices(index)
  return get_registry(index)[:,2]
end

function get_params(index)
  if is_premontaged(index) params = PARAMS_MONTAGE;
  elseif is_montaged(index) params = PARAMS_PREALIGNMENT;
  elseif is_prealigned(index) params = PARAMS_ALIGNMENT;
  elseif is_aligned(index) params = PARAMS_ALIGNMENT;
  else params = Void; println("Index $index does not correspond to a pipeline stage."); end
  return params;
end

function find_in_registry(index)
  registry = get_registry(index);
  if length(registry) == 0 return 0 end
  return findfirst(registry[:,2], index);
end

function in_registry(index::Index)
  if index[3] == 0 || index[4] == 0 return false end
  return find_in_registry(index) != 0
end

"""
Remove index from registry file & reload that registry
"""
function purge_from_registry!(index)
  # assert(is_premontaged(index))
  registry_path = get_registry_path(index)
  registry = readdlm(registry_path)
  i = find_in_registry(index)
  println("Purging $index from $registry_path")
  registry = registry[1:size(registry,1).!=i, :]
  writedlm(registry_path, registry)
  reload_registry(index)
end

function get_metadata(index)
  registry = get_registry(index);
  if find_in_registry(index) == 0 return [get_name(index) index 0 0 0 0 0 false]
  else
    return registry[find_in_registry(index), :];
  end
end

# rotated offset is just the offset of the image once rotated around [0,0]
function get_offset(index; rotated = false)
#=function get_offset(index, get_from_master=false)
	if get_from_master =#
		if myid() != IO_PROC return remotecall_fetch(IO_PROC, get_offset, index) end
	#end
	metadata = get_metadata(index);
	if rotated
    	rotation = make_rotation_matrix(metadata[3])
	if rotation == 0 return Point([0,0]) end
    	rotation_bb = snap_bb(tform_bb(sz_to_bb(Point(metadata[6:7])), rotation))
    	return Point([rotation_bb.i, rotation_bb.j])
	else
	return Point(metadata[4:5]);
      end
end

function get_rotation(index)
		if myid() != IO_PROC return remotecall_fetch(IO_PROC, get_rotation, index) end
	metadata = get_metadata(index);
	return Float64(metadata[3])
end

function get_image_size(index; rotated = false)
		if myid() != IO_PROC return remotecall_fetch(IO_PROC, get_image_size, index) end
	metadata = get_metadata(index);
	size_raw = Array{Int64, 1}(metadata[6:7]);
	if rotated
	  tform = make_rotation_matrix(get_rotation(index));
	  size_rotated = bb_to_sz(tform_bb(sz_to_bb(size_raw), tform))
	  return collect(size_rotated)
	else return size_raw
	end
end

function is_rendered(index)
	metadata = get_metadata(index);
	return Bool(metadata[8])
end

function get_preceding(index, num = 1)
  registry = get_registry(index);
  loc_in_reg = find_in_registry(index);
  if loc_in_reg == 0 return NO_INDEX; end
  if loc_in_reg <= num return NO_INDEX; end
  cur_sec = registry[loc_in_reg, 2][2]
  ind_num = 0;
  for sec_num in 1:num
    while cur_sec == registry[loc_in_reg - ind_num, 2][2]
	ind_num += 1;
	if loc_in_reg == ind_num return NO_INDEX end
    end
    cur_sec = registry[loc_in_reg - ind_num, 2][2];
  end
  return registry[loc_in_reg - ind_num, 2];
end

function get_succeeding(index, num = 1)
  registry = get_registry(index);
  loc_in_reg = find_in_registry(index);
  if loc_in_reg == 0 return NO_INDEX; end
  if loc_in_reg > size(registry, 1) - num return NO_INDEX; end
  cur_sec = registry[loc_in_reg, 2][2]
  ind_num = 0;
  for sec_num in 1:num
    while cur_sec == registry[loc_in_reg + ind_num, 2][2]
	ind_num += 1;
	if loc_in_reg == ind_num return NO_INDEX end
    end
    cur_sec = registry[loc_in_reg + ind_num, 2][2];
  end
  return registry[loc_in_reg + ind_num, 2];
end

function in_same_wafer(indexA, indexB)
  return indexB[1] == indexA[1]
end

function get_succeeding_in_wafer(index)
  new_index = get_succeeding(index)
  if in_same_wafer(new_index, index)
    return new_index
  else
    return NO_INDEX  
  end
end

function get_preceding_in_wafer(index)
  new_index = get_preceding(index)
  if in_same_wafer(new_index, index)
    return new_index
  else
    return NO_INDEX  
  end
end

function get_index_range(firstindex, lastindex; exclude_wholes_when_subs_exist = true)
  #firstindex, lastindex = match_index_stages(firstindex, lastindex)
  if is_premontaged(firstindex)
    return get_registry(firstindex)[get_range_in_registry(firstindex, lastindex), 2]
  end
    inds = filter(i->index_rank(firstindex) <= index_rank(i) <= index_rank(lastindex), get_indices(firstindex));
    if exclude_wholes_when_subs_exist
      inds_to_exclude = similar(inds, 0)
      for ind in inds
	if is_subsection(ind)
		push!(inds_to_exclude, prealigned(ind))
	end
      end
      inds = setdiff(inds, inds_to_exclude)
    end
    return inds
end

function match_index_stages(indexA, indexB)
  if is_premontaged(indexA)
    indexB = premontaged(indexB)
  elseif is_montaged(indexA)
    indexB = montaged(indexB)
  elseif is_prealigned(indexA)
    indexB = prealigned(indexB)
  elseif is_aligned(indexA)
    indexB = aligned(indexB)
  elseif is_finished(indexA)
    indexA = aligned(indexA)
    indexB = aligned(indexB)
  end
  return indexA, indexB
end

function get_range_in_registry(firstindex, lastindex)
	firstindex, lastindex = match_index_stages(firstindex, lastindex)
	if firstindex[1:2] == lastindex[1:2] && is_premontaged(firstindex) && is_premontaged(lastindex)
			return find(ind -> ind[1:2] == firstindex[1:2], REGISTRY_PREMONTAGED[:, 2]);
	end
  indices = get_indices(firstindex)
  return eachindex(indices)[firstindex .<= indices .<= lastindex]
end

"""
Find first row of the offset file that matches index and return cumulative sum 
of previous offset arrays
"""
function find_cumulative_offset(offset_file, index)
  if findfirst(offset_file[:,2], index) != 0
    return collect(sum(offset_file[1:findfirst(offset_file[:,2], index), 3:4], 1))
  else
    return [0,0]
  end
end

"""
Return zip object of an index and the index that follows it
"""
function get_sequential_index_pairs(indexA, indexB)
  indices = get_index_range(indexA, indexB)
  return zip(indices[1:end-1], indices[2:end])
end

"""
Boolean if an index is the first section in the stack
"""
function is_first_section(index)
  return get_registry(index)[1,2] == index
end

function reset_offset(index)
  update_registry(index; offset = [0,0])
end

function update_registry(index; rotation::Union{Float64, Int64} = get_rotation(index), offset::Union{Point, Array{Int64, 1}} = get_offset(index), image_size::Union{Array{Int64, 1}, Tuple{Int64, Int64}} = get_image_size(index), rendered::Bool = is_rendered(index))
  if myid() != IO_PROC return remotecall_fetch(IO_PROC, () ->update_registry(index, rotation = rotation, offset = offset, image_size = image_size, rendered = rendered)) end
  image_fn = string(get_name(index));
  registry_fp = get_registry_path(index)
  println("Updating registry for ", image_fn, " in:\n", registry_fp, ": rotation: $rotation, offset: $offset, image_size: $image_size")

  if !isfile(registry_fp)
    f = open(registry_fp, "w")
    close(f)
    registry = [image_fn, rotation, offset..., image_size..., rendered]'
  else  
    registry = readdlm(registry_fp)
    registry_line = [image_fn, rotation, offset..., image_size..., rendered]
    idx = findfirst(registry[:,1], image_fn)
    if idx != 0
      registry[idx, :] = registry_line';
    else
      registry = vcat(registry, registry_line')
    end
  end
  registry = registry[sortperm(registry[:, 1], by=parse_name), :];
  writedlm(registry_fp, registry)
  reload_registry(index);
end

function update_offsets(indices, offsets, sizes, needs_render=falses(length(indices)))
  index = indices[1]
  registry_fp = get_registry_path(index)
  filenames = map(get_name, indices)
  if !isfile(registry_fp)
    f = open(registry_fp, "w")
    close(f)
    registry = [[f, 0, o..., s..., n] for (f,o,s,n) in zip(filenames, offsets, sizes, needs_render)]
    registry = hcat(registry...)'
    println("Block line updating registry with ", registry, " in:\n", registry_fp)
  else
    registry = readdlm(registry_fp)
    locations = map(findfirst, repeated(registry), filenames)
    for (k, i) in enumerate(locations)
      println("Updating registry for ", filenames[k], " in:\n", registry_fp, ": offset is now ", offsets[k])
      registry_line = [filenames[k], 0, offsets[k]..., sizes[k]..., needs_render[k]]
      if i > 0
        registry[i,:] = registry_line
      else
        registry = vcat(registry, registry_line')
      end
    end
  end
  registry = registry[sortperm(registry[:, 1], by=parse_name), :]
  writedlm(registry_fp, registry)
  reload_registry(index)
  remotecall_fetch(IO_PROC, reload_registry, index)
end

function reload_registry(index)
  registry_fp = get_registry_path(index)
  if is_montaged(index) global REGISTRY_MONTAGED = parse_registry(registry_fp);
  elseif is_prealigned(index) global REGISTRY_PREALIGNED = parse_registry(registry_fp);
  elseif is_aligned(index) global REGISTRY_ALIGNED = parse_registry(registry_fp);
  else global REGISTRY_PREMONTAGED = parse_registry(registry_fp);
  end
end

function globalize!(pts::Points, offset::Point)
  @simd for i in 1:length(pts) @fastmath @inbounds pts[i] = pts[i] + offset; end
end
