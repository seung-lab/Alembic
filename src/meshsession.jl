function montage(firstindex::Index, lastindex::Index)
  ind_range = get_index_range(premontaged(firstindex), premontaged(lastindex))
  ind_range = unique([montaged(i[1:2]...) for i in ind_range])
  crop = [250,250]
  for index in ind_range
    ms = MeshSet(index)
    render_montaged(ms; render_full=true, render_review=false)
    # if is_flagged(ms)
    #   render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
    # else
    #   render_montaged(ms; render_full=true, render_review=false)
    # end
  end
end

function fix_montages(firstindex::Index, lastindex::Index)
  ind_range = get_index_range(premontaged(firstindex), premontaged(lastindex))
  ind_range = unique([montaged(i[1:2]...) for i in ind_range])
  for index in ind_range
    params = get_params(premontaged(index))
    ms = load(index)
    ms.properties["params"]["filter"] = params["filter"]
    ms.properties["params"]["review"] = params["review"]
    refilter!(ms);
    save(ms);
    if is_flagged(ms)
      render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
    else
      render_montaged(ms; render_full=true, render_review=false)
    end
  end
end

function remontage(index::Index)
  index = montaged(index)
  ms = MeshSet(index)
  if is_flagged(ms)
    render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
  else
    render_montaged(ms; render_full=true, render_review=false)
  end
  return ms
end

function render_montage_and_prealign(firstindex::Index, lastindex::Index)
  params = get_params(premontaged(firstindex))
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = load(index)
    ms.properties["params"]["solve"] = params["solve"]
    solve!(ms)
    save(ms)
    render_montaged(ms; render_full=true, render_review=false)
    if index > firstindex
      ms = prealign(index)
      if is_flagged(ms)
        render_prealigned(index; render_full=false, render_review=true, startindex=montaged(firstindex))
      end
    end
  end
end

function prealign(firstindex::Index, lastindex::Index; fix_start=false)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = MeshSet()
    if index==firstindex && fix_start
      println("Prealigning first section to FIXED aligned section")
      ms = prealign(index; to_fixed=fix_start)
    else 
      ms = prealign(index)
    end
    # if is_flagged(ms)
    render_prealigned(index; render_full=false, render_review=true)
    # end
  end
end

function reprealign(firstindex::Index, lastindex::Index, params)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = load(prealigned(index))
    if is_flagged(ms)
      try
        reset_offset(index)
        ms = prealign(montaged(index); params=params)
        if is_flagged(ms)
          render_prealigned(index; render_full=false, render_review=true)
        end
      catch e
        log_error(prealigned(index); fn="match_error_log", comment=e)
      end
    end
  end
end

function align(firstindex::Index, lastindex::Index)
  ms = MeshSet(firstindex, lastindex; solve=false)
  render_aligned_review(ms)
  split_meshset(ms)
end

function align(index_list)
  for (indexA, indexB) in index_list
    align(indexA, indexB)
  end
end

function align_over_missing_tiles(index_depth_list)
  for (firstindex, lastindex, depth) in index_depth_list
    params = get_params(firstindex);
    params["match"]["depth"] = depth;
    ms = MeshSet(firstindex, lastindex; solve=false, params=params);
    render_aligned_review(ms)
  end
end

function solve_align(firstindex::Index, lastindex::Index)
  parent_name = get_name(firstindex, lastindex)
  ms = concat_meshset(parent_name)
  save(ms)
  solve!(ms)
  save(ms)
  split_meshset(ms)
end

function refilter!(firstindex::Index, lastindex::Index, params=get_params(firstindex))
  parent_name = get_name(firstindex, lastindex)
  for i = 1:count_children(parent_name)
    refilter!(firstindex, lastindex, i, params)
  end
end

function refilter!(firstindex::Index, lastindex::Index, ind::Int64, params=get_params(firstindex))
  parent_name = get_name(firstindex, lastindex)
  ms = load_split(parent_name, ind)
  ms.properties["params"]["filter"] = params["filter"]
  ms.properties["params"]["review"] = params["review"]
  refilter!(ms)
  save(ms)
end

function copy_through_first_section(index::Index)
  img = get_image(montaged(index))

  function write_image(index, img)
    fn = string(get_name(index), ".h5")
    dir = get_dir(index)

    update_offset(index, [0,0], size(img))
    println("Writing image:\n\t", fn)
    # @time imwrite(stage["img"], joinpath(dir, fn))
    f = h5open(joinpath(dir, fn), "w")
    chunksize = min(1000, min(size(img)...))
    @time f["img", "chunk", (chunksize,chunksize)] = img
    f["offset"] = [0,0]
    f["scale"] = 1.0
    close(f)
  end

  write_image(prealigned(index), img)
  write_image(aligned(index), img)
end

function check_and_view_flags!(firstindex::Index, lastindex::Index)
  params = get_params(firstindex)
  review_params = params["review"]

  flagged_indices = []

  if is_montaged(firstindex) || is_prealigned(firstindex)
    if is_montaged(firstindex)
      func = montaged
    elseif is_prealigned(firstindex)
      func = prealigned
    end
    for index in get_index_range(montaged(firstindex), montaged(lastindex))
      index = func(index)
      ms = load(index)
      ms.properties["params"]["review"] = get_params(get_index(ms.meshes[1]))["review"]
      check!(ms)
      save(ms)
      if is_flagged(ms)
        push!(flagged_indices, index)
      end
    end
  elseif is_aligned(firstindex)
    parent_name = get_name(prealigned(firstindex), prealigned(lastindex))
    for k in 1:count_children(parent_name)
      ms = load_split(parent_name, k)
      ms.properties["params"]["review"] = review_params
      check!(ms)
      save(ms)
      if is_flagged(ms)
        push!(flagged_indices, k)
      end
    end
  end

  return flagged_indices

end

"""
Cycle through index range and return list of flagged meshset indices
"""
function view_flags(firstindex::Index, lastindex::Index)
  flagged_indices = []

  if is_montaged(firstindex) || is_prealigned(firstindex)
    if is_montaged(firstindex)
      func = montaged
    elseif is_prealigned(firstindex)
      func = prealigned
    end
    for index in get_index_range(montaged(firstindex), montaged(lastindex))
      meshset = load(func(index))
      # check!(meshset)
      if is_flagged(meshset)
        push!(flagged_indices, index)
      end
    end
  elseif is_aligned(firstindex)
    parent_name = get_name(prealigned(firstindex), prealigned(lastindex))
    for k in 1:count_children(parent_name)
      ms = load_split(parent_name, k)
      if is_flagged(ms)
        push!(flagged_indices, k)
      end
    end
  end

  return flagged_indices
end

function write_reviews_as_needed(firstindex::Index, lastindex::Index)
  if is_montaged(firstindex)
    for index in get_index_range(montaged(firstindex), montaged(lastindex))
      ms = load(index)
      if is_flagged(ms)
        render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
      end
    end
  elseif is_prealigned(firstindex)
    for index in get_index_range(montaged(firstindex), montaged(lastindex))
      ms = load(prealigned(index))
      if is_flagged(ms)
        render_prealigned(index; render_full=false, render_review=true)
      end
    end
  end
end

function negative_nans(a)
  return isnan(a) ? -Inf : a
end

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

function premontage(firstindex::Index, lastindex::Index)
  scale = 0.5

  indices = get_index_range(firstindex, lastindex)
  sections = unique([ind[1:2] for ind in indices])

  for sec in sections
    println("Premontaging tiles in section $sec")
    tile_indices = spiral_sort(get_index_range(premontaged(sec...), premontaged(sec...)))
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
    moving_indices = [i[1] for i in translations]
    for index in tile_indices
      k = findfirst(i -> i==index, moving_indices)
      if k > 0
        mindex, findex, moving_bb, fixed_bb, dv = translations[k]
        offset = get_offset(findex) - get_offset(fixed_bb) + get_offset(mindex) + dv
        update_offset(mindex, offset, get_image_size(mindex));
      end
    end
    save_premontage_review(tile_indices[1])
  end
end

function find_translation(overlap)
  mindex, findex, moving_bb, fixed_bb, scale = overlap
  dv = find_translation_offset(mindex, findex, moving_bb, fixed_bb, scale=scale)
  return mindex, findex, moving_bb, fixed_bb, dv
end

function find_translation(moving_index::Index, fixed_index::Index, mbb=get_bb(moving_index), fbb=get_bb(fixed_index); scale=0.5)
  println("Translating $moving_index to $fixed_index")
  overlap_bb = get_bb(moving_index) - get_bb(fixed_index)
  println("Overlap area: ", get_area(overlap_bb), " px^2")
  moving = get_slice(moving_index, overlap_bb, scale, is_global=true)
  fixed = get_slice(fixed_index, overlap_bb, scale, is_global=true) 
  xc = normxcorr2(moving, fixed; shape="full")
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
  update_offset(index, [vertical, horizontal], (h,w))
  save_premontage_review(index)
end

#### ASSUMES HOMOGENEOUS TILE SIZE
#function premontage(wafer_range::UnitRange{Int64})
function premontage_to_overview(wafer, start, tile_size = [8000,8000], scale = 0.05)
#  for wafer in wafer_range
    premontaged_path = get_path(premontaged(wafer, 1))
    dir,name = splitdir(premontaged_path)
    tiles = sort_dir(dir, ".h5");
    tiles = filter(x->contains(x,"Tile"), tiles)
    tile_indices = sort(map(parse_name, tiles))
    tile_indices = filter(x->x[1] == wafer, tile_indices)
    section_range = start:tile_indices[end][2]
  for sec in section_range
    premontaged_path = get_path(premontaged(wafer, sec))
    dir,name = splitdir(premontaged_path)
    println("$wafer, $sec")
    tiles = sort_dir(dir, ".h5");
    tiles = filter(x->contains(x,"Tile"), tiles)
    tile_indices = sort(map(parse_name, tiles))
    tile_indices = filter(x->x[1:2] == (wafer,sec), tile_indices)

    if length(tile_indices) == 0 continue; end

    fixed_indices = similar(tile_indices, 0)
    images = map(get_image, tile_indices)
    #images = map(Images.restrict, images);
    #images = map(Images.restrict, images);
    #images = map(Images.restrict, images);

    sigma = [1,1]

    overview_image = get_image(overview(wafer, sec))
    function find_patch_locs(cur_image, overview_image)
	cur_image = imscale(cur_image, scale)[1]
	xc = normxcorr2(cur_image, overview_image)
	rm, ind = findmax(xc)
	i_max, j_max = ind2sub(xc, ind);
	return [round(Int64, i_max / scale), round(Int64, j_max / scale)]
    end

    offsets = pmap(find_patch_locs, images, repeated(overview_image))
    for i in 1:length(tile_indices)
      update_offset(tile_indices[i], offsets[i], tile_size)
    end
    
# end
 end
end

"""
Write any errors to a log file
"""
function log_error(index::Index; fn="render_error_log", comment="")
  ts = parse(Dates.format(now(), "yymmddHHMMSS"))
  dir = get_dir(index)
  path = joinpath(dir, string(fn, ".txt"))
  new_row = [ts, index, comment]'
  if !isfile(path)
    f = open(path, "w")
    close(f)
    log = new_row
  else  
    log = readdlm(path)
    log = vcat(log, new_row)
  end
  log = log[sortperm(log[:, 1]), :]
  println("Logging error:\n", path)
  writedlm(path, log)
end
