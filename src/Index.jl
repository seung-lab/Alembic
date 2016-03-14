function overview(index)
	return (index[1], index[2], OVERVIEW_INDEX, OVERVIEW_INDEX);
end

function premontaged(index)
	return (index[1], index[2], PREMONTAGED_INDEX, PREMONTAGED_INDEX);
end

function montaged(index)
	return (index[1], index[2], MONTAGED_INDEX, MONTAGED_INDEX);
end

function prealigned(index)
	return (index[1], index[2], PREALIGNED_INDEX, PREALIGNED_INDEX);
end

function aligned(index)
	return (index[1], index[2], ALIGNED_INDEX, ALIGNED_INDEX);
end

function finished(index)
  return (index[1], index[2], FINISHED_INDEX, FINISHED_INDEX);
end
###

function overview(wafer, section)
	return (wafer, section, OVERVIEW_INDEX, OVERVIEW_INDEX);
end

function premontaged(wafer, section)
	return (wafer, section, PREMONTAGED_INDEX, PREMONTAGED_INDEX);
end

function montaged(wafer, section)
	return (wafer, section, MONTAGED_INDEX, MONTAGED_INDEX);
end

function prealigned(wafer, section)
	return (wafer, section, PREALIGNED_INDEX, PREALIGNED_INDEX);
end

function aligned(wafer, section)
	return (wafer, section, ALIGNED_INDEX, ALIGNED_INDEX);
end

function finished(wafer, section)
  return (wafer, section, FINISHED_INDEX, FINISHED_INDEX);
end

###

function is_overview(index)
    if index[3:4] == (OVERVIEW_INDEX, OVERVIEW_INDEX)   return true else return false end
end

function is_premontaged(index)
    if index[3] > 0 && index[4] > 0 return true else return false end
end

function is_montaged(index)
    if index[3:4] == (MONTAGED_INDEX, MONTAGED_INDEX)   return true else return false end
end

function is_prealigned(index)
    if index[3:4] == (PREALIGNED_INDEX, PREALIGNED_INDEX)   return true else return false end
end

function is_aligned(index)
    if index[3:4] == (ALIGNED_INDEX, ALIGNED_INDEX) return true else return false end
end

function is_finished(index)
    if index[3:4] == (FINISHED_INDEX, FINISHED_INDEX) return true else return false end
end

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

function is_preceding(A_index, B_index)
  if (is_premontaged(A_index)) || (is_premontaged(B_index))	return false end
	if A_index[1:2] == get_preceding(B_index)[1:2] return true; end
	return false;
end


function get_overview_index(index)
  return (index[1:2]..., OVERVIEW_INDEX, OVERVIEW_INDEX)
end

function get_registry(index)
  if is_premontaged(index) registry = REGISTRY_PREMONTAGED;
  elseif is_montaged(index) registry = REGISTRY_MONTAGED;
  elseif is_prealigned(index) registry = REGISTRY_PREALIGNED;
  elseif is_aligned(index) registry = REGISTRY_ALIGNED;
  elseif is_registered(index) registry = REGISTRY_ALIGNED;
  else registry = Void; println("Index $index does not correspond to a pipeline stage."); end
  return registry; 
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
  return findfirst(registry[:,2], index);
end

function get_metadata(index)
  registry = get_registry(index);
  return registry[find_in_registry(index), :];
end

function get_offset(index)
#=function get_offset(index, get_from_master=false)
	if get_from_master =#
		if myid() != 1 return remotecall_fetch(1, get_offset, index) end
	#end
	metadata = get_metadata(index);
	return Point(metadata[3:4]);
end

function get_image_sizes(index)
	metadata = get_metadata(index);
	return metadata[5:6];
end

function get_preceding(index)
  registry = get_registry(index);
  loc_in_reg = find_in_registry(index);
  if loc_in_reg == 0 return NO_INDEX; end
  #if loc_in_reg == 1 return index; end
  if loc_in_reg == 1 return NO_INDEX; end
  return registry[loc_in_reg - 1, 2];
end

function get_succeeding(index)
  registry = get_registry(index);
  loc_in_reg = find_in_registry(index);
  if loc_in_reg == 0 return NO_INDEX; end
  if loc_in_reg == 1 return NO_INDEX; end
  return registry[loc_in_reg + 1, 2];
end

function get_index_range(first_index, last_index)
  if is_finished(first_index) || is_finished(last_index)
    first_index, last_index = aligned(first_index), aligned(last_index)
  end
	ran = get_registry(last_index)[get_range_in_registry(first_index, last_index), 2];
	ran[1] = first_index;
	return ran;
end

function get_range_in_registry(indexA, indexB)
	
	if indexA[1:2] == indexB[1:2] && is_premontaged(indexA) && is_premontaged(indexB)
			return find(ind -> ind[1:2] == indexA[1:2], REGISTRY_PREMONTAGED[:, 2]);
	end
	if indexA[3] != indexB[3] || indexA[4] != indexB[4]
		println("The indices are from different pipeline stages or from different sections for premontage. Defaulting to the latter's stage...");
	end
	return find_in_registry((indexA[1:2]..., indexB[3:4]...)):find_in_registry(indexB);
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
