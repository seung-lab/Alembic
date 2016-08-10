"""
Create 2D integer spiral ordering from a midpoint
http://stackoverflow.com/questions/3706219/algorithm-for-iterating-over-an-outward-spiral-on-a-discrete-2d-grid-from-the-or

If NUM_POINTS not specified, assumes: 
* bounding box origin is [1,1]
* bounding box is square off the maximum dimension
* center is calculated with ceil (so, bounding width & height are even)
"""
function create_spiral(center::Array{Int64,1}, NUM_POINTS=(2*maximum(center))^2)
  assert(*((center .> [0,0])...))
  spiral = Array{Array{Int64,1}, 1}()
  i, j = center
  push!(spiral, center)
  di, dj = [1, 0]
  seg_length = 1
  seg_passed = 0
  for k in 1:NUM_POINTS-1
    i += di
    j += dj
    push!(spiral, [i,j])
    seg_passed += 1
    if seg_passed == seg_length
      seg_passed = 0
      di, dj = -dj, di # 90 degree counterclockwise
      if dj == 0
        seg_length += 1
      end
    end
  end
  return spiral
end

function spiral_sort(indices)
  rows = map(getindex, indices, repeated(3))
  cols = map(getindex, indices, repeated(4))
  row_range = maximum(rows) + minimum(rows)
  col_range = maximum(cols) + minimum(cols)
  center = [ceil(Int64, row_range/2), ceil(Int64, col_range/2)]
  spiral = create_spiral(center)
  section = indices[1][1:2]
  spiral_indices = [(section..., s...) for s in spiral]
  return filter(i->i in indices, spiral_indices)
end

function premontage(index::Index)
  scale = 0.5
  println("Premontaging tiles in $(premontaged(index))")
  tile_indices = spiral_sort(get_index_range(premontaged(index), premontaged(index)))
  bbs = map(get_bb, tile_indices)
  overlaps = []
  for k in 2:length(tile_indices)
    mindex = tile_indices[k]
    println("Finding overlap for $mindex")
    moving_bb = bbs[k]
    fixed_bbs = bbs[1:k-1]
    overlap_area = [get_area(fixed_bb - moving_bb) for fixed_bb in fixed_bbs]
    max_idx = selectperm(overlap_area, 1, rev=true, by=negative_nans)
    if overlap_area[max_idx] > 100*100
      findex = tile_indices[max_idx]
      fixed_bb = fixed_bbs[max_idx]
      push!(overlaps, (mindex, findex, moving_bb, fixed_bb, scale))
    else
      println("No overlap found")
    end
  end
  
  translations = pmap(find_translation, overlaps)
  offsets = map(get_offset, tile_indices)
  sizes = map(get_image_size, tile_indices)
  moving_indices = [i[1] for i in translations]
  for (i, index) in enumerate(tile_indices)
    k = findfirst(i->i==index, moving_indices)
    if k > 0
      mindex, findex, moving_bb, fixed_bb, dv = translations[k]
      j = findfirst(i->i==findex, tile_indices)
      fixed_offset = offsets[j]
      moving_offset = offsets[i]
      offsets[i] = fixed_offset - get_offset(fixed_bb) + moving_offset + dv
      # update_offset(mindex, offset, get_image_size(mindex));
      # offset = get_offset(findex) - get_offset(fixed_bb) + get_offset(mindex) + dv
    end
  end
  update_offsets(tile_indices, offsets, sizes)
  save_premontage_review(tile_indices[1])
end

function premontage(firstindex::Index, lastindex::Index)
  indices = get_index_range(firstindex, lastindex)
  sections = unique([ind[1:2] for ind in indices])
  for sec in sections
    premontage(premontaged(sec...))
  end
end

function find_translation(overlap)
  mindex, findex, moving_bb, fixed_bb, scale = overlap
  dv = find_translation_offset(mindex, findex, moving_bb, fixed_bb, scale=scale)
  return mindex, findex, moving_bb, fixed_bb, dv
end

function find_translation(moving_index::Index, fixed_index::Index, mbb=get_bb(moving_index), fbb=get_bb(fixed_index); scale=0.5, highpass_sigma = 10)
  println("Translating $moving_index to $fixed_index")
  overlap_bb = get_bb(moving_index) - get_bb(fixed_index)
  println("Overlap area: ", get_area(overlap_bb), " px^2")
  moving = get_slice(moving_index, overlap_bb, scale, is_global=true)
  fixed = get_slice(fixed_index, overlap_bb, scale, is_global=true) 
	if highpass_sigma != 0
		highpass_sigma = highpass_sigma * scale
	  	moving = Array{Float64, 2}(moving);
	  	moving_g = copy(moving)
		@fastmath Images.imfilter_gaussian_no_nans!(moving_g, [highpass_sigma, highpass_sigma])
		elwise_sub!(moving, moving_g);

	  	fixed = Array{Float64, 2}(fixed);
	  	fixed_g = copy(fixed)
		@fastmath Images.imfilter_gaussian_no_nans!(fixed_g, [highpass_sigma, highpass_sigma])
		elwise_sub!(fixed, fixed_g);
      end
  xc = normxcorr2_preallocated(moving, fixed; shape="full")
  return xc, moving, fixed
end

function find_translation_offset(moving_index::Index, fixed_index::Index, mbb=get_bb(moving_index), fbb=get_bb(fixed_index); scale=0.5)
  xc, moving, fixed = find_translation(moving_index, fixed_index, mbb, fbb, scale=scale)
  rm, ind = findmax(xc)
  i_max, j_max = ind2sub(xc, ind);
  i = 1/scale * (i_max - median(1:size(xc, 1)))
  j = 1/scale * (j_max - median(1:size(xc, 2)))
  return [i,j]
end

function get_major_overlaps(indices)
  bbs = [i=>get_bb(i) for i in indices]
  neighbors = map(get_cardinal_neighbors, indices)
end

function check_premontage(index::Index)
  indices = get_index_range(premontaged(index), premontaged(index))
end

function set_relative_offset(moving_index::Index, fixed_index::Index, rel_offset)

  offset = get_offset(fixed_index) + rel_offset
  println("$moving_index, $fixed_index, $rel_offset")
  # update_offset(moving_index, offset, get_image_size(moving_index))
end

function snap(index::Index, overlap=[880, 690])
  h, w = get_image_size(index)
  offset = get_offset(index)
  vertical = 0
  horizontal = 0
  below = get_below(index)
  if below != NO_INDEX
    vertical = get_offset(below)[1] - (h-overlap[1])
  else
    above = get_above(index)
    if above != NO_INDEX
      vertical = get_offset(above)[1] + (h-overlap[1])
    end
  end
  left = get_left(index)
  if left != NO_INDEX
    horizontal = get_offset(left)[2] + (w-overlap[2])
  else
    right = get_right(index)
    if right != NO_INDEX
      horizontal = get_offset(right)[2] - (w-overlap[2])
    end
  end
  println("$index, $offset -> ($vertical, $horizontal)")
  update_registry(index; offset = [vertical, horizontal], image_size = [h,w])
  save_premontage_review(index)
end

#### ASSUMES HOMOGENEOUS TILE SIZE
#function premontage(wafer_range::UnitRange{Int64})
function premontage_to_overview(wafer, start, tile_size = [8000,8000], scale = 0.05)
#  for wafer in wafer_range
    dir = PREMONTAGED_DIR_PATH
    tiles = sort_dir(dir, ".h5");
    tile_indices = sort(map(parse_name, tiles))
    tile_indices = filter(x->x[1] == wafer, tile_indices)
    section_range = start:tile_indices[end][2]
  for sec in section_range
    dir = PREMONTAGED_DIR_PATH
    println("$wafer, $sec")
    tiles = sort_dir(dir, ".h5");
    tile_indices = sort(map(parse_name, tiles))
    tile_indices = filter(x->x[1:2] == (wafer,sec), tile_indices)

    if length(tile_indices) == 0 continue; end

    fixed_indices = similar(tile_indices, 0)
    images = map(get_image, tile_indices)
    #images = map(Images.restrict, images);
    #images = map(Images.restrict, images);
    #images = map(Images.restrict, images);

    sigma = [1,1]

    overview_image = get_image(get_path(overview(wafer, sec), ".tif"))
    function find_patch_locs(cur_image, overview_image)
	cur_image = imscale(cur_image, scale)[1]
	xc = normxcorr2(cur_image, overview_image)
	rm, ind = findmax(xc)
	i_max, j_max = ind2sub(xc, ind);
	return [round(Int64, i_max / scale), round(Int64, j_max / scale)]
    end

    offsets = pmap(find_patch_locs, images, repeated(overview_image))
    for i in 1:length(tile_indices)
      update_registry(tile_indices[i]; offset = offsets[i], image_size = tile_size)
    end
    
# end
 end
end
