function montage(firstindex::Index, lastindex::Index)
  ind_range = get_index_range(premontaged(firstindex), premontaged(lastindex))
  ind_range = unique([montaged(i[1:2]...) for i in ind_range])
  for index in ind_range
    premontage(premontaged(index))
    ms = MeshSet(index)
    render_montaged(ms; render_full=true, render_review=true, flagged_only=true)
    calculate_stats(ms)
    reset_cache()
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
    ms.properties[:params][:filter] = params[:filter]
    ms.properties[:params][:review] = params[:review]
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
    ms.properties[:params][:solve] = params[:solve]
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

function prealign(firstindex::Index, lastindex::Index)
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    ms = prealign(index)
    calculate_stats(ms)
    try
      render(ms, review=true)
    catch
      render_prematch_review(index)
    end
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
    params["match"][:depth] = depth;
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
  ms.properties[:params][:filter] = params[:filter]
  ms.properties[:params][:review] = params[:review]
  refilter!(ms)
  save(ms)
end

function copy_through_first_section(index::Index)
  img = get_image(montaged(index))

  function write_image(index, img)
    fn = string(get_name(index), ".h5")
    dir = get_dir_path(index)

    update_registry(index; offset = [0,0], image_size = size(img))
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
  review_params = params[:review]

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
      ms.properties[:params][:review] = get_params(get_index(ms.meshes[1]))[:review]
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
      ms.properties[:params][:review] = review_params
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
Write any errors to a log file
"""
function log_error(index::Index; fn="render_error_log", comment="")
  ts = parse(Dates.format(now(), "yymmddHHMMSS"))
  dir = get_dir_path(index)
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


function d_montage(z, threshold=3)
  index = montaged(1,z); 
  ms = load(MeshSet, index);
  # for i in 1:length(ms.matches)
  #   calculate_post_statistics!(ms, i)
  #   flag!(ms.matches[i], "norm_post", >, 3)
  # end
  # save(ms)
  s = load("stats", index); 
  c = filter_stats(s, ["max_post"]); 
  println(c[c[:,2] .> 3, :])
  i = sort(map(parse, c[c[:,2] .> 3, 1]))[1]
  inspect(ms, i)
  return ms, s, c
end
