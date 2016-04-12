function montage(firstindex::Index, lastindex::Index)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = MeshSet(index)
    if is_flagged(ms)
      render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
    else
      render_montaged(ms; render_full=true, render_review=false)
    end
  end
end

function render_flagged_montages(firstindex::Index, lastindex::Index)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = load(index)
    if is_flagged(ms)
      unflag!(ms)
      save(ms)
      render_montaged(ms; render_full=true, render_review=false)
    end
  end
end

function remontage(firstindex::Index, lastindex::Index, params)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = load(index)
    check!(ms)
    if is_flagged(ms)
      ms = MeshSet(index; params=params)
      if is_flagged(ms)
        render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
      else
        render_montaged(ms; render_full=true, render_review=false)
      end
    end
  end
end

function prealign(firstindex::Index, lastindex::Index; start_to_fixed=false)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = MeshSet()
    # try
      if index==firstindex
        ms = prealign(index; to_fixed=true)
      else 
        ms = prealign(index)
      end
      if is_flagged(ms)
        render_prealigned(index; render_full=false, render_review=true)
      end
    # catch e
    #   log_error(prealigned(index); fn="match_error_log", comment=e)
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

function align(firstindex::Index, lastindex::Index; fix_first=false)
  ms = MeshSet(firstindex, lastindex; solve=false, fix_first=fix_first)
  render_aligned_review(ms)
  split_meshset(ms)
end

function align(index_list)
  for (indexA, indexB) in index_list
    align(indexA, indexB)
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


#### ASSUMES HOMOGENEOUS TILE SIZE
#function premontage(wafer_range::UnitRange{Int64})
function premontage(wafer, start)
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
    images = map(Images.restrict, images);
    images = map(Images.restrict, images);

    #fix the first one
    update_offset(tile_indices[1], [0, 0], [8000,8000]);
    push!(fixed_indices, tile_indices[1])
    oset = Points(length(tile_indices))
    oset[1] = [0,0];

    sigma = [1,1]

    for index in tile_indices[2:end]
	target_index = fixed_indices[findfirst(this -> is_adjacent(index, this), fixed_indices)]

	xc = normxcorr2(images[findfirst(ind -> ind == index, tile_indices)], images[findfirst(ind -> ind == target_index, tile_indices)]; shape="full")
	xcd = xc - Images.imfilter_gaussian(xc, sigma);

	rm, ind = findmax(xcd)
	i_max, j_max = ind2sub(xcd, ind);
	i_diff = 4 * (i_max - median(1:size(xcd, 1)))
	j_diff = 4 * (j_max - median(1:size(xcd, 2)))

#    	update_offset(index, [i_diff, j_diff] + get_offset(target_index), [size(images[findfirst(ind -> ind == index, tile_indices)])...]);
    	update_offset(index, [i_diff, j_diff] + oset[findfirst(ind -> ind == target_index, tile_indices)], [8000,8000]);
    	oset[findfirst(ind -> ind == index, tile_indices)] = [i_diff, j_diff] + oset[findfirst(ind -> ind == target_index, tile_indices)]
    	push!(fixed_indices, index)
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