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

function is_adjacent(A_index, B_index)
  if abs(A_index[3] - B_index[3]) + abs(A_index[4] - B_index[4])  == 1 return true; end
  return false;
end

function is_diagonal(A_index, B_index)
  if abs(A_index[3] - B_index[3]) + abs(A_index[4] - B_index[4]) == 2 && A_index[3] != B_index[3] && A_index[4] != B_index[4] return true; end
  return false;
end

function is_preceding(A_index, B_index)
	if A_index == get_preceding(B_index) && A_index[3:4] == B_index[3:4] return true; end
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
  else registry = Void; println("Index $index does not correspond to a pipeline stage."); end
  return registry; 
end

function get_params(index)
  if is_premontaged(index) params = PARAMS_MONTAGE;
  elseif is_montaged(index) params = PARAMS_PREALIGNMENT;
  elseif is_prealigned(index) params = PARAMS_ALIGNMENT;
  elseif is_aligned(index) params = Void;
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
  if loc_in_reg == 1 return index; end
  #if loc_in_reg == 1 return NO_INDEX; end
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
	return get_registry(first_index)[get_range_in_registry(first_index, last_index), 2];
end

function get_range_in_registry(indexA, indexB)
	
	if indexA[1:2] == indexB[1:2] && is_premontaged(indexA) && is_premontaged(indexB)
			return find(ind -> ind[1:2] == indexA[1:2], REGISTRY_PREMONTAGED[:, 2]);
	end
	if indexA[3] != indexB[3] || indexA[4] != indexB[4]
		println("The indices are from different pipeline stages or from different sections for premontage. Aborting.");
		return nothing;
	end
	return find_in_registry(indexA):find_in_registry(indexB);
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
Generate list of indices between indexA and indexB
"""
function get_index_range(indexA, indexB)
  return get_registry(indexA)[get_range_in_registry(indexA, indexB), 2]
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
