function is_adjacent(A::Index, B::Index, cross_wafer_guess = false)
	if sum(abs([A...] - [B...])) == 1
		return true;
	elseif cross_wafer_guess && (A[1]-B[1]==1 && A[2]==1) || (A[1]-B[1] == -1 && B[2]==1)
		return true;
	end
	return false;
end

function is_adjacent(Am::Mesh, Bm::Mesh)
  if abs(Am.index[3] - Bm.index[3]) + abs(Am.index[4] - Bm.index[4]) + abs(Am.index[2] - Bm.index[2]) == 1 return true; end
  return false;
end

function is_diagonal(Am::Mesh, Bm::Mesh)
  if abs(Am.index[3] - Bm.index[3]) + abs(Am.index[4] - Bm.index[4]) == 2 && Am.index[3] != Bm.index[3] && Am.index[4] != Bm.index[4] return true; end
  return false;
end

function get_overview_index(index::Index)
  return (index[1:2]..., OVERVIEW_INDEX, OVERVIEW_INDEX)
end

function get_registry(index::Index)
  if is_premontaged(index) registry = PREMONTAGED_OFFSETS;
  elseif is_montaged(index) registry = MONTAGED_OFFSETS;
  elseif is_prealigned(index) registry = PREALIGNED_OFFSETS;
  elseif is_aligned(index) registry = ALIGNED_OFFSETS;
  else registry = Void; println("Index $index does not correspond to a pipeline stage."); end
  return registry; 
end

function get_params(index::Index)
  if is_premontaged(index) params = PARAMS_MONTAGE;
  elseif is_montaged(index) params = PARAMS_PREALIGNMENT;
  elseif is_prealigned(index) params = PARAMS_ALIGNMENT;
  elseif is_aligned(index) params = Void;
  else params = Void; println("Index $index does not correspond to a pipeline stage."); end
  return params;
end

function find_in_registry(index::Index)
  registry = get_registry(index);
  return findfirst(registry[:,2], index);
end

function find_preceding(index::Index)
  registry = get_registry(index);
  loc_in_reg = find_in_registry(index);
  if loc_in_reg == 0 println("Index $index not found in the registry."); return NO_INDEX; end
  if loc_in_reg == 1 return ("Index $index is the first entry in the registry."); NO_INDEX; end
  return registry[loc_in_reg - 1, 2];
end

function get_range_in_registry(first_index, last_index)
	if first_index[3] != last_index[3] || first_index[4] != last_index[4]
		println("The indices are from different pipeline stages. Aborting.");
		return Void;
	end
	
	return find_in_registry(first_index):find_in_registry(last_index);
end

"""
Find first row of the offset file that matches index and return
of previous offset arrays
"""
function find_offset(offset_file, index)
  if findfirst(offset_file[:,2], index) != 0
    return collect(offset_file[findfirst(offset_file[:,2], index), 3:4])
  else
    return [0,0]
  end
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
Find appropriate offset file and pull out the offset array for the index
"""
function load_offset(index)
  if is_montaged(index)
    return find_offset(MONTAGED_OFFSETS, index)
  elseif is_aligned(index)
    return find_offset(PREALIGNED_OFFSETS, index)
  else
    return [0,0]
  end
end
"""
INCOMPLETE

Lazy function to generate list of indices between indexA and indexB

Fixes to include:
* switch between wafers
* check indexA[3:4] against indexB[3:4]
"""
function create_index_range(indexA, indexB)
  return [(indexA[1], i, indexA[3:4]...) for i in indexA[2]:indexB[2]]
end

"""
Return zip object of an index and the index that follows it
"""
function create_sequential_index_pairs(indexA, indexB)
  indices = create_index_range(indexA, indexB)
  return zip(indices[1:end-1], indices[2:end])
end

"""
Test if an index is the first section in the stack
"""
function is_first_section(index)
  return index[1] == 1 && index[2] == 1
end
