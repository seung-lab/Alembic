function inspect(index::Index, match_ind=1)
  meshset = load(index)
  inspect(meshset, match_ind)
end

function inspect(firstindex::Index, lastindex::Index, match_ind=1)
  meshset = load(firstindex, lastindex)
  inspect(meshset, match_ind)
end

function inspect(meshset::MeshSet, match_ind=1)
  src_index = meshset.matches[match_ind].src_index
  dst_index = meshset.matches[match_ind].dst_index
  println("\n", get_name(meshset), ": ", (src_index, dst_index), " @ ", match_ind, " / ", length(meshset.matches))
  imgc, img2, vectors, params = view_matches(meshset, match_ind)
  enable_inspection(imgc, img2, meshset, match_ind, vectors, params)
end

function get_next(meshset::MeshSet, match_ind=1)
end

function get_previous(meshset::MeshSet, match_ind=1)
end

"""
Return next match in the matches list that's flagged
"""
function get_next_flag(meshset::MeshSet, match_ind=1)
  matches = meshset.matches
  if match_ind > length(matches) & match_ind > 0
    return nothing
  end
  if is_flagged(matches[match_ind])
    return meshset, match_ind
  else
    return get_next_flag(meshset, match_ind+1)
  end
end

"""
Return previous meshset & match index in the matches list that's flagged
"""
function get_previous_flag(meshset::MeshSet, match_ind)
  matches = meshset.matches
  if match_ind > length(matches) & match_ind > 0
    if 
    return nothing
  end
  if is_flagged(matches[match_ind])
    return meshset, match_ind
  else
    return get_previous_flag(meshset, match_ind-1)
  end
end

"""
Display match index match_ind as overlayed images with points
"""
function view_matches(meshset::MeshSet, match_ind)
  match = meshset.matches[match_ind]
  indexA = match.src_index
  indexB = match.dst_index

  path = get_review_path(indexB, indexA)
  if !isfile(path)
    path = get_review_path(indexA, indexB)
  end
  img_orig = h5read(path, "img")
  offset = h5read(path, "offset")
  println("offset: ", offset)
  scale = h5read(path, "scale")

  # Add border to the image for vectors that extend
  pad = 400
  img = zeros(UInt32, size(img_orig,1)+pad*2, size(img_orig,2)+pad*2)
  img[pad:end-pad-1, pad:end-pad-1] = img_orig

  params = deepcopy(meshset.properties["params"]["match"])
  params["offset"] = offset - pad
  params["scale"] = scale
  params["match_index"] = match_ind
  params["vector_scale"] = 4
  params["post_matches"] = false
  params["dist"] = 90
  params["sigma"] = 3

  view_inspection_statistics(match, params["search_r"])

  imgc, img2 = view(img, pixelspacing=[1,1])
  vectors = make_vectors(meshset, match_ind, params)
  show_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1))
  update_annotations(imgc, img2, match, vectors)
  return imgc, img2, vectors, params
end

function make_vectors(meshset, match_ind, params)
  scale = params["scale"]
  offset = params["offset"]
  factor = params["vector_scale"]
  println("Vector scale: ", factor)
  src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences(meshset, match_ind)
  if params["post_matches"]
    src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences_post(meshset, match_ind)
  end
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vecs = offset_matches(vectorsA, vectorsB, offset)
  vectors = [hcat(vecs[1]...); hcat(vecs[2]...)]
  return change_vector_lengths([hcat(vecs[1]...); hcat(vecs[2]...)], factor)
end

"""
Display the images involved in the blockmatching at one point in a Match object
"""
function view_blockmatch(match, ind, params)
  println(match.correspondence_properties[ind])
  src_patch, src_pt, dst_patch, dst_pt, xc, offset = get_correspondence_patches(match, ind)
  block_r = params["block_r"]
  search_r = params["search_r"]
  beta = 0.5
  N=size(xc, 1)
  M=size(xc, 2)
  if !contains(gethostname(), "seunglab") && !ON_AWS
    fig = figure("correlogram")
    surf([i for i=1:N, j=1:M], [j for i=1:N, j=1:M], xc, cmap=get_cmap("hot"), 
                            rstride=10, cstride=10, linewidth=0, antialiased=false)
    grid("on")
    title("match $ind")
  end
  println("max r-value: ", maximum(xc))
  xc_image = xcorr2Image(xc)
  xc_beta_mask = xc .> beta*maximum(xc)
  xc_beta = xcorr2Image(xc .* xc_beta_mask)
  # xc_image = padimage(xc_image, block_r, block_r, block_r, block_r, 1)
  hot = create_hot_colormap()
  xc_color = apply_colormap(xc_image, hot)
  xc_beta_color = apply_colormap(xc_beta, hot)
  fused_img, _ = imfuse(dst_patch, [0,0], src_patch, offset)

  src = convert(Array{UInt32,2}, src_patch).<< 8
  dst = convert(Array{UInt32,2}, dst_patch).<< 16

  cgrid = canvasgrid(2,3; pad=10)
  opts = Dict(:pixelspacing => [1,1])

  imgc, img2 = view(cgrid[1,1], src; opts...)
  imgc, img2 = view(cgrid[2,1], dst; opts...)
  imgc, img2 = view(cgrid[2,2], fused_img; opts...)
  imgc, img2 = view(cgrid[1,2], xc_color'; opts...)
  imgc, img2 = view(cgrid[1,3], xc_beta_color'; opts...)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  return win
end

function xcorr2Image(xc)
  b = xc' / maximum(xc)
  b[b .> 1] = 1
  b[b .< 0] = 0
  # b = b.*b / maximum(b.*b) * 254 + 1
  b = b * 254 + 1
  # b = (b + 1) ./ 2 * 254 + 1
  b[isnan(b)] = 0
  return round(UInt8, b)
  # return round(UInt8, (b + 1) ./ 2 * 254 + 1)
end

function create_hot_colormap()
  return vcat(linspace(RGB(0,0,0), RGB(1,0,0), 100), 
              linspace(RGB(1,0,0), RGB(1,1,0), 100), 
              linspace(RGB(1,1,0), RGB(1,1,1), 55))
end

function apply_colormap{T}(img, colormap::Array{T})
  new_img = zeros(T, size(img)...)
  for i in eachindex(img, new_img)
    @inbounds new_img[i] = colormap[img[i]]
  end
  return new_img
end

"""
Extend ImageView to inspect & remove matches
"""
function enable_inspection(imgc::ImageView.ImageCanvas, 
                                img2::ImageView.ImageSlice2d, 
                                meshset, match_ind, vectors, params)
  println("Enable inspection")
  matches = meshset[match_ind]
  c = canvas(imgc)
  win = Tk.toplevel(c)
  c.mouse.button2press = (c, x, y) -> brushtool_start(c, x, y, (c, bb) -> remove_contained_points(imgc, img2, matches, vectors, bb))
  bind(c, "<Button-3>", (c, x, y)->inspect_match(imgc, img2, 
                                                        parse(Int, x), 
                                                        parse(Int, y), 
                                                        matches,
                                                        vectors, 
                                                        params))
  bind(c, "<Control-Button-3>", (c, x, y)->remove_match(imgc, img2, 
                                                                parse(Int, x), 
                                                                parse(Int, y), 
                                                                matches,
                                                                vectors))
  bind(win, "f", path->flag_inspection(imgc, img2, matches))
  bind(win, "<Control-z>", path->undo_match_filter(imgc, img2, matches, vectors))
  bind(win, "<Escape>", path->disable_inspection(imgc, img2))
  bind(win, "<Destroy>", path->disable_inspection(imgc, img2))
  # bind(win, "z", path->end_edit())
  bind(win, "c", path->compare_filter(imgc, img2, matches, vectors))
  bind(win, ",", path->decrease_distance_filter(imgc, img2, matches, vectors, params))
  bind(win, ".", path->increase_distance_filter(imgc, img2, matches, vectors, params))
  bind(win, "n", path->decrease_sigma_filter(imgc, img2, matches, vectors, params))
  bind(win, "m", path->increase_sigma_filter(imgc, img2, matches, vectors, params))
  bind(win, "=", path->increase_vectors(imgc, img2, meshset, matches, vectors, params))
  bind(win, "-", path->decrease_vectors(imgc, img2, meshset, matches, vectors, params))
  bind(win, "p", path->switch_pre_to_post(imgc, img2, meshset, matches, vectors, params))
  bind(win, "s", path->save_inspection(meshset, matches))
  bind(win, "<Return>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=true, flag=false))
  # bind(win, "<Return>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=false, flag=false))
  # bind(win, "<Return>", path->go_to_next_flag(imgc, img2, meshset, match_ind; forward=true, flag=true))
  # bind(win, "<Return>", path->go_to_previous_flag(imgc, img2, meshset, match_ind; forward=false, flag=true))
end

function disable_inspection(imgc::ImageView.ImageCanvas, img2::ImageView.ImageSlice2d)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  println("Disable inspection")
  bind(c, "<Button-3>", path->path)
  bind(c, "<Control-Button-3>", path->path)
  bind(win, "f", path->path)
  bind(win, "<Destroy>", path->path)
  bind(win, "<Escape>", path->path)
  bind(win, "<Control-z>", path->path)
  bind(win, "c", path->path)
  bind(win, ",", path->path)
  bind(win, ".", path->path)
  bind(win, "=", path->path)
  bind(win, "-", path->path)
  bind(win, "p", path->path)
  bind(win, "s", path->path)
  bind(win, "<Return>", path->path)
  destroy(win)
end

function save_inspection(meshset, matches)
  println("meshset saved")
  set_reviewed!(matches)
  save(meshset)
end

function go_to_next_inspection(imgc, img2, meshset, match_ind)
  disable_inspection(imgc, img2)

  if stage == "montage"
    index, match_ind = calls
    match_ind += 1
    if match_ind > length(meshset.matches)
      index = get_succeeding(index)
      match_ind = 1
    end
    inspect_montages(index, match_ind)
  elseif stage == "prealignment"
    src_index, dst_index = calls
    src_index, dst_index = dst_index, get_succeeding(dst_index)
    inspect_prealignments(src_index, dst_index)
  elseif stage == "alignment"
    firstindex, lastindex, src_index, dst_index = calls
    src_index, dst_index = dst_index, get_succeeding(dst_index)
    inspect_alignments(firstindex, lastindex, src_index, dst_index)
  elseif stage == "meshset"
    match_ind = calls
    match_ind += 1
    if match_ind <= length(meshset.matches)
      inspect_meshset(meshset, match_ind)
    end
  end
end

"""
Preliminary method to test filters in the GUI
"""
function compare_filter(imgc, img2, match, vectors)
  filter = (0.5, >, 5)
  inds_to_filter = Array{Any, 1}()
  attributes = get_properties(match, filter[1])
  push!(inds_to_filter, find(i -> filter[2](i, filter[3]), attributes))
  inds_to_filter = union(inds_to_filter...)
  rejected_inds = get_rejected_indices(match)

  manual_removed_sigma_did_not = setdiff(rejected_inds, inds_to_filter)
  reverse_mask = collect(setdiff(Set(1:length(attributes)), manual_removed_sigma_did_not))

  clear_filters!(match)
  filter_manual!(match, reverse_mask)
  update_annotations(imgc, img2, match, vectors)
end

"""
Invert to just the bad matches
"""
function show_removed(imgc, img2, match, vectors)
  all_inds = Set(1:length(match.src_points))
  rejected_inds = get_rejected_indices(match)

  accepted_inds = collect(setdiff(all_inds, rejected_inds))

  clear_filters!(match)
  filter_manual!(match, accepted_inds)
  update_annotations(imgc, img2, match, vectors)
end

function update_annotations(imgc, img2, match, vectors)
  mask = get_filtered_indices(match)
  for an in Base.values(imgc.annotations)
    if :pts in fieldnames(an.data)
      an.data.pts = vectors[1:2, mask]
    elseif :lines in fieldnames(an.data)
      an.data.lines = vectors[:, mask]
    end
  end
  ImageView.redraw(imgc)
end

function inspect_match(imgc, img2, x, y, matches, vectors, params, prox=0.0125)
  # prox: 0.0125 = 100/8000
  indices = get_filtered_indices(matches)
  lines = vectors[:, indices]

  xu, yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
  xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
  limit = (img2.zoombb.xmax - img2.zoombb.xmin) * prox 
  annidx = find_idx_of_nearest_pt(lines[1:2,:], [xi, yi], limit)
  if annidx > 0
    idx = indices[annidx]
    ptA = vectors[1:2,idx] # - params["src_offset"]
    ptB = vectors[3:4,idx] # - params["dst_offset"]
    println(idx, ": ", ptA, ", ", ptB)
    bm_win = view_blockmatch(matches, idx, params)
    detect_blockmatch_removal(imgc, img2, bm_win, matches, idx, vectors)
  end
end

"""
Find index of location in array with point closest to the given point (if there
exists such a unique point within a certain pixel limit).

Args:

* pts: 2xN array of coordinates
* pt: 2-element coordinate array of point in interest
* limit: number of pixels that nearest must be within

Returns:

* index of the nearest point in the pts array

  idx = find_idx_of_nearest_pt(pts, pt, limit)
"""
function find_idx_of_nearest_pt(pts, pt, limit)
    d = sum((pts.-pt).^2, 1).^(1/2)
    idx = eachindex(d)'[d .< limit]
    if length(idx) == 1
        return idx[1]
    else
        return 0
    end
end

function remove_contained_points(imgc, img2, matches, vectors, bb)
  indices = get_filtered_indices(matches)
  pts = vectors[1:2, indices]
  within_bb = vcat((bb.xmin .<= pts[1,:] .<= bb.xmax) & (bb.ymin .<= pts[2,:] .<= bb.ymax)...)
  contained_indices = indices[within_bb]
  println("Brushtool removed ", length(contained_indices), " points")
  filter_manual!(matches, contained_indices)
  update_annotations(imgc, img2, matches, vectors)
end

function detect_blockmatch_removal(imgc::ImageView.ImageCanvas, 
                                      img2::ImageView.ImageSlice2d, 
                                      win, matches, idx, vectors)
  bind(win, "<Delete>", path->remove_blockmatch_from_patch_window(win, imgc, 
                                                  img2, matches, idx, vectors))
end

function remove_blockmatch_from_patch_window(win, imgc, img2, matches, idx, vectors)
  bind(win, "<Delete>", path->path)
  destroy(win)
  filter_manual!(matches, idx)
  update_annotations(imgc, img2, matches, vectors)
end

function remove_match(imgc, img2, x, y, matches, vectors, prox=0.0125)
  indices = get_filtered_indices(matches)
  lines = vectors[:, indices]

  xu, yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
  xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
  limit = (img2.zoombb.xmax - img2.zoombb.xmin) * prox
  annidx = find_idx_of_nearest_pt(lines[1:2,:], [xi, yi], limit)
  if annidx > 0
    idx = indices[annidx]
    pt = vectors[1:2,idx]
    println("Manually removed ", idx, ": ", pt)
    filter_manual!(matches, idx)
    update_annotations(imgc, img2, matches, vectors)
  end
end

function undo_match_filter(imgc, img2, matches, vectors)
  println("Undo")
  idx = undo_filter!(matches)
  update_annotations(imgc, img2, matches, vectors)
end

function increase_distance_filter(imgc, img2, matches, vectors, params)
  params["dist"] += 10
  filter_match_distance(imgc, img2, matches, vectors, params["dist"])
end

function decrease_distance_filter(imgc, img2, matches, vectors, params)
  params["dist"] = max(params["dist"]-10, 0)
  filter_match_distance(imgc, img2, matches, vectors, params["dist"])
end

function increase_sigma_filter(imgc, img2, matches, vectors, params)
  params["sigma"] += 0.5
  filter_match_sigma(imgc, img2, matches, vectors, params["sigma"])
end

function decrease_sigma_filter(imgc, img2, matches, vectors, params)
  params["sigma"] = max(params["sigma"]-0.5, 0)
  filter_match_sigma(imgc, img2, matches, vectors, params["sigma"])
end

function filter_match_distance(imgc, img2, matches, vectors, dist)
  # hack to test if a match_distance filter was just implemented
  if length(matches.filters) > 0
    if matches.filters[end]["type"] == "norm"
      undo_filter!(matches)
    end
  end
  println("Distance filter @ ", dist)
  filter!(matches, "norm", >, dist)
  update_annotations(imgc, img2, matches, vectors)
end

function filter_match_sigma(imgc, img2, matches, vectors, sigma)
  if !ON_AWS
  # hack to test if a match_distance filter was just implemented
  if length(matches.filters) > 0
    if matches.filters[end]["type"] == "sigma_5"
      undo_filter!(matches)
    end
  end
  println("Sigma filter @ ", sigma)
  filter!(matches, "sigma_5", >, sigma)
  else
  if length(matches.filters) > 0
    if matches.filters[end]["type"] == .5
      undo_filter!(matches)
    end
  end
  println("Sigma filter @ ", sigma)
  filter!(matches, .5, >, sigma)
  end
  update_annotations(imgc, img2, matches, vectors)
end

function flag_inspection(imgc, img2, matches)
  println("Flag inspection!")
  flag!(matches)
end

function increase_vectors(imgc, img2, meshset, matches, vectors, params)
  params["vector_scale"] *= 1.1
  k = params["match_index"]
  vectors = make_vectors(meshset, k, params)
  update_annotations(imgc, img2, matches, vectors)
end

function decrease_vectors(imgc, img2, meshset, matches, vectors, params)
  params["vector_scale"] = max(params["vector_scale"]*0.9, 0.1)
  k = params["match_index"]
  vectors = make_vectors(meshset, k, params)
  update_annotations(imgc, img2, matches, vectors)
end

function switch_pre_to_post(imgc, img2, meshset, matches, vectors, params)
  params["post_matches"] = !params["post_matches"]
  println(params["post_matches"] ? "Using post matches" : "Using pre matches" )
  k = params["match_index"]
  vectors = make_vectors(meshset, k, params)
  update_annotations(imgc, img2, matches, vectors)
end

function show_filtered_points(imgc, img2, meshset, matches, vectors, params)
  # display filtered matches in different color/opacity
end

function get_inspection_path(username, stage_name)
  return joinpath(INSPECTION_DIR, string(stage_name, "_inspection_", username, ".txt"))
end

"""
Combine stack error log files of all tracers to create one log matrix
"""
function compile_review_logs(stage)
  tracers = ["hmcgowan", "bsilverman", "merlinm", "kpw3", "mmoore", "dih"]
  logs = []
  for tracer in tracers
    path = get_inspection_path(tracer, stage)
    if isfile(path)
      push!(logs, readdlm(path))
    end
  end
  return vcat(logs...)
end

function show_montage_inspection_seam_progress(index)
  meshset = load(index)
  seam_count = length(meshset.matches)
  seams = 1:seam_count
  reviewed = falses(seam_count)
  flagged = falses(seam_count)
  for (i, match) in enumerate(meshset.matches)
    reviewed[i] = is_reviewed(match)
    flagged[i] = is_flagged(match)
  end
  fig = figure("montage_inspection_seam_progress $index")
  subplot(211)
  title("is_reviewed")
  plot(seams, reviewed, ".")
  subplot(212)
  title("is_flagged")
  plot(seams, flagged, ".")
end

function show_montage_inspection_section_progress(start=0, finish=9999999)
  firstindex, lastindex = (2,1,-2,-2), (8,173,-2,-2)
  montages = get_index_range(firstindex, lastindex)
  section_count = length(montages)
  sections = 1:section_count
  reviewed = falses(section_count)
  flagged = falses(section_count)
  for (i, index) in enumerate(montages)
    if start < i < finish
      meshset = load(index) 
      reviewed[i] = is_reviewed(meshset)
      flagged[i] = is_flagged(meshset)
    end
  end
  fig = figure("montage_inspection_section_progress")
  subplot(211)
  title("is_reviewed")
  plot(sections, reviewed, ".")
  subplot(212)
  title("is_flagged")
  plot(sections, flagged, ".")
end

function show_alignment_inspection_progress(firstindex, lastindex)
  # firstindex, lastindex = (1,167,-3,-3), (2,149,-3,-3)
  parent_name = string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned")
  splits_count = count_children(parent_name)
  println(parent_name)
  matches = 1:splits_count
  reviewed = falses(splits_count)
  flagged = falses(splits_count)
  for i = 1:splits_count
    meshset = load_split(parent_name, i) 
    reviewed[i] = is_reviewed(meshset)
    flagged[i] = is_flagged(meshset)
  end
  fig = figure("$parent_name alignment_inspection_progress")
  subplot(211)
  title("is_reviewed")
  plot(matches, reviewed, ".")
  subplot(212)
  title("is_flagged")
  plot(matches, flagged, ".")
end

"""
Display time-based chart of montage point inspections for the tracers
"""
function show_montage_inspection_progress(show_stats=false)
  path = joinpath(MONTAGED_DIR, "review")
  montages = get_index_range((1,1,-2,-2), (8,173,-2,-2))
  factor = 100
  x = 1:(factor*length(montages))
  y = zeros(Int64, factor*length(montages))
  logs = compile_review_logs("montage")
  tracers = unique(logs[:, 2])
  fig, ax = PyPlot.subplots()
  for tracer in tracers
    log = logs[logs[:,2] .== tracer, :]
    log_time = round(Int64, (log[:,1] % 10^6) / 100)
    meshset_indices = [(parse_index(l)[1:2]..., -2, -2) for l in log[:,3]]
    log_k = map(x -> findfirst(montages, x)*factor, meshset_indices)
    log_k += log[:,5]
    y[log_k] = log_time
    ax[:plot](log_k, log_time, ".", label=tracer)
    if show_stats
      times = log[sortperm(log[:, 1]), 1]
      dt = [times[i] - times[i-1] for i in 2:length(times)]
      dt_med = median(dt)
      reviews = size(log,1)
      first_day = round(Int, minimum(times) / 1000000) % 10
      last_day = round(Int, maximum(times) / 1000000) % 10
      work_hrs = 8
      println(tracer, ": ", dt_med, " s, ", reviews, ", ", round(Int, reviews/(last_day-first_day+1)/8), " seams/hr")
    end
  end
  ax[:legend](loc="best")
  x_wo = x[y.==0]
  x_wo = x_wo[0 .< x_wo % factor .<= 84]
  PyPlot.plot(x_wo, zeros(Int64, length(x_wo)), ".")
  seams = [join([logs[i,3], logs[i,4]]) for i in 1:size(logs, 1)]
  println("Reviewed seams: ", length(unique(seams)), " / ", sum(y.>0) + length(x_wo))
end

"""
Turn index into a ten-based integer
"""
function create_index_id(ind::Index)
  return abs(ind[1])*10^7 + abs(ind[2])*10^4 + abs(ind[3])*10^2 + abs(ind[4])
end

"""
Create integer from two indices
"""
function create_index_id(indexA::Index, indexB::Index)
  return create_index_id(indexA)*10^8 + create_index_id(indexB)
end

"""
Filter logs by most recent entry for each seam
"""
function get_most_recent_logs(logs)
  logs = logs[sortperm(logs[:,1]), :]
  indicesA = [parse_index(l) for l in logs[:,3]]
  indicesB = [parse_index(l) for l in logs[:,4]]
  seam_id = [create_index_id(a,b) for (a,b) in zip(indicesA, indicesB)]
  logs = hcat(logs, seam_id)
  logs = logs[sortperm(logs[:,end]), :]
  last_id = vcat([logs[i,end] != logs[i+1,end] for i=1:(size(logs,1)-1)], true)
  return logs[last_id .== true, :]
end

function get_meshset_with_edits(meshset, ind, logs)
  meshset_indices = [(parse_index(l)[1:2]..., ind[3:4]...) for l in logs[:,3]]
  logs = logs[meshset_indices .== ind, :]
  for i in 1:size(logs,1)
    match = meshset.matches[logs[i,5]]
    inds_to_filter = [logs[i,6]]
    if typeof(inds_to_filter[1]) != Int64
      inds_to_filter = readdlm(IOBuffer(logs[i,6]), ',', Int)
    end
    if inds_to_filter[1] != 0
      filter_manual!(match, inds_to_filter)
    end
    set_reviewed!(match)
    match.properties["review"]["author"]["by"] = logs[i,2]
  end
  return meshset
end

function update_montage_meshsets(waferA, secA, waferB, secB)
  indices = get_index_range(montaged(waferA,secA), montaged(waferB,secB))
  logs = compile_review_logs("montage")
  logs = get_most_recent_logs(logs)
  for index in indices
    meshset = load(index)
    meshset = get_meshset_with_edits(meshset, index, logs)
    save(meshset)
    # solve!(meshset, method="elastic")
    # save(meshset)
  end
end

function update_prealignment_meshsets(waferA, secA, waferB, secB)
  logs = compile_review_logs("prealignment")
  logs = get_most_recent_logs(logs)
  index_pairs = get_sequential_index_pairs((waferA,secA,-2,-2), (waferB,secB,-2,-2))
  for (indexA, indexB) in index_pairs
    println(indexB, indexA)
    meshset = load(indexB, indexA)
    meshset = get_meshset_with_edits(meshset, indexB, logs)
    solve!(meshset, method="regularized")
    save(meshset)
  end
end

function get_review_path(src_index, dst_index=(0,0,0,0))
  prefix = "review"
  dir = ALIGNED_DIR
  ind = indices_to_string(src_index, dst_index)
  if is_premontaged(src_index) || is_premontaged(dst_index)
    dir = MONTAGED_DIR
    ind = string(join(src_index, ","), "-", join(dst_index, ","))
  elseif is_montaged(src_index) || is_montaged(dst_index)
    dir = PREALIGNED_DIR
  elseif is_prealigned(src_index) || is_prealigned(dst_index)
    dir = ALIGNED_DIR
  end
  fn = string(prefix, "_", ind, ".h5")
  return joinpath(dir, "review", fn)
end

function indices_to_string(indexA, indexB)
  if indexB[1] == 0
    return join(indexA[1:2], ",")
  end
  return string(join(indexA[1:2], ","), "-", join(indexB[1:2], ","))
end

function is_dv_near_search_r(match::Match, sr, factor=0.05)
  dv = hcat(get_filtered_properties(match, "dv")...)
  return sum(abs(dv) .> (1-factor)*sr) > 0
end

function mark_suspicious(meshset::MeshSet)
  search_r = get_dfs(meshset.properties["params"], "search_r")
  for (i, match) in enumerate(meshset.matches)
    print("\n", i, "\t", match.src_index, match.dst_index)
    num_matches = length(match.src_points)
    if num_matches > 20
      r_max = get_filtered_properties(match, "r_max")
      num_filtered = length(r_max)
      if num_filtered == 0
        print("\tall filtered")
        flag!(match)
      else
        dyn_range = get_filtered_properties(match, "src_normalized_dyn_range")
        dv = get_filtered_properties(match, "dv")
        
        if sum(abs(hcat(dv...)) .> (search_r*0.9)) > 0
          print("\tdv")
          flag!(match)
        else
          print("\t  ")
        end
        if sum(dyn_range .< 0.5) > 0
          print("\tdyn")
          flag!(match)
        else
          print("\t   ")
        end
        if (num_filtered / num_matches < 0.2) & (num_matches > 30)
          print("\tpts")
          flag!(match)
        else
          print("\t   ")
        end
      end
    end
  end
end

function view_dvs(index::Index)
  meshset = load(index)
  sr = get_param(meshset, "search_r")
  k = 1
  n = 0
  fig = figure("dv $index $n", figsize=(20,20))
  for (i, match) in enumerate(meshset.matches)
    if rem(i,9) == 0
      n += 1
      fig = figure("dv $index $n", figsize=(20,20))
      k = 1
    end
    dv_all = hcat(get_properties(match, "dv")...)
    dv = hcat(get_filtered_properties(match, "dv")...)
    dv_bounds = vcat([[i -sr] for i=-sr:5:sr]...,
                      [[i sr] for i=-sr:5:sr]...,
                      [[sr i] for i=-sr:5:sr]...,
                      [[-sr i] for i=-sr:5:sr]...)
    subplot(330+k)
    if length(dv_all) > 1
      plt[:scatter](dv_all[2,:], -dv_all[1,:], color="#990000")
      plt[:scatter](dv_bounds[:,1], dv_bounds[:,2], color="#0f0f0f", marker="+")
    end
    if length(dv) > 1
      plt[:scatter](dv[2,:], -dv[1,:], color="#009900")
    end
    grid("on")
    title("dv $i")
    k += 1
  end
end

function view_sigma_plots(meshset::MeshSet, range=0:5:100)
  props = sort(collect(keys(meshset.matches[1].correspondence_properties[1])))
  colors = ["#000055", "#005500", "#550000", "#0000dd", "#00dd00", "#dd0000"]
  c_ind = 0
  for k in props
    if contains(k, "sigma")
      c_ind += 1
      precision = []
      recall = []
      for r in range
        fp, fn, tp, total = eval_filters(meshset, [(k, >, r, 0)] ,:);
        push!(precision, (100 * tp / (fp + tp)))
        push!(recall, (100 * tp / (fn + tp)))
      end
      plt[:scatter](precision, recall, label=k, color=colors[c_ind], alpha=0.5)
    end
  end
  legend(loc="upper right",fancybox="true")
  grid("on")
end


function view_dv_dispersion(index::Index)
  meshset = load(index)
  fig = figure("dv $index", figsize=(20,20))
  for (i, match) in enumerate(meshset.matches)
    dv_all = hcat(get_properties(match, "dv")...)
    dv_all_ind = ones(size(dv_all, 2))*i
    dv = hcat(get_filtered_properties(match, "dv")...)
    dv_ind = ones(size(dv, 2))*i
    # println(i, " ", size(dv_all, 2), " ", size(dv, 2)) 
    if size(dv_all, 2) > 1
      subplot(211)
      plt[:scatter](dv_all_ind, dv_all[1,:], color="#bb0000", alpha=0.5)
      plt[:scatter](dv_ind, dv[1,:], color="#00bb00", alpha=0.5)
      grid("on")
      title("x")
    end
    if size(dv, 2) > 1
      subplot(212)
      plt[:scatter](dv_all_ind, dv_all[2,:], color="#bb0000", alpha=0.5)
      plt[:scatter](dv_ind, dv[2,:], color="#00bb00", alpha=0.5)
      grid("on")
      title("y")
    end
  end
end

function view_property_histogram(firstindex::Index, lastindex::Index, property_name, nbins=20)
  attr = []
  indrange = get_index_range(firstindex, lastindex)
  for index in indrange
    meshset = load(index)
    for match in meshset.matches
      attr = vcat(attr, get_properties(match, property_name))
    end
  end
  masknan = map(!, map(isnan, attr))
  attr = attr .* masknan
  fig = figure("histogram")
  p = plt[:hist](attr, nbins)
  grid("on")
  title(property_name)
  return p
end

function view_property_histogram(meshset::MeshSet, property_name, nbins=20)
  attr = []
  for match in meshset.matches
    attr = vcat(attr, get_properties(match, property_name))
  end
  fig = figure("histogram")
  p = plt[:hist](attr, nbins)
  grid("on")
  title(property_name)
  return p
end


function view_property_scatter(match, property_name)
  attr = get_properties(match, property_name)
  xy = Array{Int, 2}()
  if length(attr[1]) == 2
    xy = hcat(attr...)
  else
    xy = hcat(1:length(attr), attr)'
  end
  fig = figure("scatter")
  p = plt[:scatter](xy[1,:], xy[2,:])
  grid("on")
  title(property_name)
  return p
end

function view_property_histogram(match::Match, property_name; filtered=true, nbins=20)
  color="#990000"
  properties_func = get_properties
  if filtered
    properties_func = get_filtered_properties
    color="#009900"
  end
  attr = properties_func(match, property_name)
  if length(attr) > 1
    p = plt[:hist](attr, nbins, color=color)
    grid("on")
    title(property_name)
    return p
  end
end

function view_dv_dispersion(match, sr; filtered=true)
  color="#990000"
  properties_func = get_properties
  if filtered
    properties_func = get_filtered_properties
    color="#009900"
  end
  dv = hcat(properties_func(match, "dv")...)
  dv_bounds = vcat([[i -sr] for i=-sr:5:sr]...,
                    [[i sr] for i=-sr:5:sr]...,
                    [[sr i] for i=-sr:5:sr]...,
                    [[-sr i] for i=-sr:5:sr]...)
  if length(dv) > 1
    p = plt[:scatter](dv[2,:], -dv[1,:], color=color, alpha=0.5)
    p = plt[:scatter](dv_bounds[:,1], dv_bounds[:,2], color="#0f0f0f", marker="+", alpha=0.5)
    grid("on")
    title("dv")
    return p
  end
end

function view_property_spatial_scatter(match, property_name; filtered=true, factor=1)
  color="#990000"
  correspondences_func = get_correspondences
  properties_func = get_properties
  if filtered
    correspondences_func = get_filtered_correspondences
    properties_func = get_filtered_properties
    color="#009900"
  end
  pts = hcat(correspondences_func(match)[1]...)
  attr = properties_func(match, property_name)
  if length(attr) > 1
    p = plt[:scatter](pts[2,:], -pts[1,:], s=attr*factor, color=color, alpha=0.5)
    grid("on")
    title(property_name)
    return p
  end
end

function view_inspection_statistics(match, sr)
  fig = figure("image pair statistics", figsize=(20,20))
  PyPlot.clf()
  subplot(231)
  view_property_histogram(match, "r_max"; filtered=false, nbins=20)
  view_property_histogram(match, "r_max"; filtered=true, nbins=20)
  subplot(232)
  view_dv_dispersion(match, sr; filtered=false)
  view_dv_dispersion(match, sr; filtered=true)
  subplot(233)
  view_property_spatial_scatter(match, "r_max"; filtered=false, factor=100)
  view_property_spatial_scatter(match, "r_max"; filtered=true, factor=100)
  subplot(234)
  view_property_spatial_scatter(match, 0.5; filtered=false, factor=1)
  view_property_spatial_scatter(match, 0.5; filtered=true, factor=1)
  subplot(235)
  view_property_spatial_scatter(match, "src_kurtosis"; filtered=false, factor=10)
  view_property_spatial_scatter(match, "src_kurtosis"; filtered=true, factor=10)
  subplot(236)
  view_property_spatial_scatter(match, "src_normalized_dyn_range"; filtered=false, factor=100)
  view_property_spatial_scatter(match, "src_normalized_dyn_range"; filtered=true, factor=100)
  fig[:canvas][:draw]()
end
