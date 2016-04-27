function inspect(index::Index, match_ind=0)
  meshset = load(index)
  inspect(meshset, match_ind)
end

function inspect(firstindex::Index, lastindex::Index, match_ind=0)
  if is_prealigned(firstindex) && is_prealigned(lastindex)
    k = match_ind
    if match_ind == 0
      k = 1
    end 
    meshset = load_split(get_name(firstindex, lastindex), k)
  else
    meshset = load(firstindex, lastindex)
  end
  inspect(meshset, match_ind)
end

function inspect(meshset::MeshSet, match_ind=0)
  if match_ind == 0
    meshset, match_ind = get_next_flagged_match(meshset, 0)
  end
  k = match_ind
  name = get_name(meshset)
  num_string = string(" @ ", match_ind, " / ", length(meshset.matches))
  if is_prealigned(meshset.meshes[1].index) || is_prealigned(meshset.meshes[2].index)
    k = 1
    name = string("#", get_name(meshset), " ", get_parent(meshset))
    num_string = ""
  end
  src_index = meshset.matches[k].src_index
  dst_index = meshset.matches[k].dst_index
  println("\n", name, ": ", (src_index, dst_index), num_string)
  imgc, img2, vectors, params = view_match(meshset, k)
  enable_inspection(imgc, img2, meshset, k, vectors, params)
end

"""
Get next match in the wafer (could span meshsets)
"""
function get_next_match(meshset::MeshSet, match_ind=1)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[end].index
  if !is_aligned(meshset)
    matches = meshset.matches
    if 1 <= match_ind < length(matches)
      return meshset, match_ind+1
    elseif match_ind == length(matches)
      if is_premontaged(firstindex)
        indexA = premontaged(get_succeeding_in_wafer(montaged(lastindex)))
        indexB = firstindex
      elseif is_montaged(firstindex)
        indexA = get_succeeding_in_wafer(firstindex)
        indexB = firstindex
      end

      if indexA[1:2] == (0,0) || indexB[1:2] == (0,0)
        return nothing, nothing
      else
        meshset = load(indexA, indexB)
        match_ind = 1
        return meshset, match_ind
      end
    end

  elseif is_aligned(meshset)
    parent_name = get_parent(meshset)
    match_ind += 1
    if match_ind <= count_children(parent_name)
      return load_split(parent_name, match_ind), match_ind
    else
      return nothing, nothing
    end
  end
  return meshset, 1
end

"""
Get previous match in the wafer (could span meshsets)
"""
function get_previous_match(meshset::MeshSet, match_ind=1)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[end].index
  if !is_aligned(meshset)
    matches = meshset.matches
    if 1 < match_ind <= length(matches)
      return meshset, match_ind-1
    elseif match_ind == 1
      if is_premontaged(firstindex)
        indexA = premontaged(get_preceding_in_wafer(montaged(firstindex)))
        indexB = firstindex
      elseif is_montaged(firstindex)
        indexA = lastindex
        indexB = get_preceding_in_wafer(lastindex)
      end

      if indexA[1:2] == (0,0) || indexB[1:2] == (0,0)
        return nothing, nothing
      else
        meshset = load(indexA, indexB)
        match_ind = length(meshset.matches)
        return meshset, match_ind
      end
    end

  elseif is_aligned(meshset)
    parent_name = get_parent(meshset)
    match_ind -= 1
    if match_ind > 0
      return load_split(parent_name, match_ind), match_ind
    else
      return nothing, nothing
    end
  end
  return meshset, 1
end

"""
Return next match in the matches list that's flagged
"""
function get_next_flagged_match(meshset::MeshSet, match_ind=1)
  meshset, match_ind = get_next_match(meshset, match_ind)
  if meshset == nothing || match_ind == nothing
    return nothing, nothing
  end
  k = match_ind
  if is_aligned(meshset)
    k = 1
  end
  if is_flagged(meshset.matches[k])
    return meshset, match_ind
  else
    return get_next_flagged_match(meshset, match_ind)
  end
end

"""
Return previous meshset & match index in the matches list that's flagged
"""
function get_previous_flagged_match(meshset::MeshSet, match_ind)
  meshset, match_ind = get_previous_match(meshset, match_ind)
  if meshset == nothing || match_ind == nothing
    return nothing, nothing
  end
  k = match_ind
  if is_aligned(meshset)
    k = 1
  end
  if is_flagged(meshset.matches[k])
    return meshset, match_ind
  else
    return get_previous_flagged_match(meshset, match_ind)
  end
end

"""
Display match index match_ind as overlayed images with points
"""
function view_match(meshset::MeshSet, match_ind)
  match = meshset.matches[match_ind]
  indexA = match.src_index
  indexB = match.dst_index

  path = get_review_path(indexB, indexA)
  if !isfile(path)
    path = get_review_path(indexA, indexB)
    if !isfile(path)
      error("NO REVIEW IMAGE CREATED")
    end
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
  params["vector_scale"] = 40
  params["dist"] = 90
  params["sigma"] = 20
  if is_montaged(meshset)
    params["vector_scale"] = 4
    params["sigma"] = 7
  end
  params["post_matches"] = false # is_prealigned(indexA)

  if USE_PYPLOT
    view_inspection_statistics(match, params["search_r"])
  end

  println("make image")
  imgc, img2 = view(img, pixelspacing=[1,1])
  # resize(imgc, 400, 600)
  vectors = make_vectors(meshset, match_ind, params)
  show_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1))
  update_annotations(imgc, img2, match, vectors, params)

  c = canvas(imgc)
  win = Tk.toplevel(c)
  fnotify = ImageView.Frame(win)
  lastrow = 2
  ImageView.grid(fnotify, lastrow+=1, 1, sticky="ew")
  xypos = ImageView.Label(fnotify)
  imgc.handles[:pointerlabel] = xypos
  ImageView.grid(xypos, 1, 1, sticky="ne")
  ImageView.set_visible(win, true)
  c.mouse.motion = (path,x,y)-> updatexylabel(xypos, imgc, img2, x, y, params["offset"]..., scale)

  return imgc, img2, vectors, params
end

function resize(imgc, w, h)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  Tk.set_size(win, w, h)
end

"""
Make vectors for display, scaling and offsetting appropriately
"""
function make_vectors(meshset, match_ind, params)
  scale = params["scale"]
  offset = params["offset"]
  src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences(meshset, match_ind)
  if params["post_matches"]
    src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences_post(meshset, match_ind)
  end
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vecs = offset_matches(vectorsA, vectorsB, offset)
  vectors = [hcat(vecs[1]...); hcat(vecs[2]...)]
  return [vectors[2,:]; vectors[1,:]; vectors[4,:]; vectors[3,:]]
end

"""
Scale a list of 2-element points
"""
function scale_matches(pts, scale)
  return [x*scale for x in pts]
end

"""
Maintain vector start point, but adjust end point for more prominent visual
"""
function change_vector_lengths(vectors, k)
  v = []
  if length(vectors) > 0
    v = [vectors[1,:]; 
          vectors[2,:]; 
          (vectors[3,:]-vectors[1,:])*k + vectors[1,:]; 
          (vectors[4,:]-vectors[2,:])*k + vectors[2,:]]
  end
  return v
end

function offset_matches(src_pts, dst_pts, offset)
  src_pts = [x - offset for x in src_pts]
  dst_pts = [x - offset for x in dst_pts]
  return src_pts, dst_pts
end

function update_annotations(imgc, img2, match, vectors, params)
  v = change_vector_lengths(vectors, params["vector_scale"])
  mask = get_filtered_indices(match)
  for an in Base.values(imgc.annotations)
    if :pts in fieldnames(an.data)
      an.data.pts = v[1:2, mask]
    elseif :lines in fieldnames(an.data)
      an.data.lines = v[:, mask]
    end
  end
  ImageView.redraw(imgc)
end

"""
Display the images involved in the blockmatching at one point in a Match object
"""
function view_blockmatch(match, match_ind, params)
  for (k, v) in match.correspondence_properties[match_ind]
    println(k, ":\t", v)
  end
  src_patch, src_pt, dst_patch, dst_pt, xc, offset = get_correspondence_patches(match, match_ind)
  block_r = params["block_r"]
  search_r = params["search_r"]
  beta = 0.5
  N=size(xc, 1)
  M=size(xc, 2)
  if USE_PYPLOT
    fig = figure("correlogram")
    surf([i for i=1:N, j=1:M], [j for i=1:N, j=1:M], xc, cmap=get_cmap("hot"), 
                            rstride=10, cstride=10, linewidth=0, antialiased=false)
    grid("on")
    title("match $match_ind")
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
  matches = meshset.matches[match_ind]
  c = canvas(imgc)
  win = Tk.toplevel(c)
  c.mouse.button2press = (c, x, y) -> brushtool_start(c, x, y, (c, bb) -> remove_contained_points(imgc, img2, matches, vectors, bb, params))
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
                                                                vectors,
                                                                params))
  bind(win, "f", path->toggle_flag(imgc, img2, matches))
  bind(win, "<Control-z>", path->undo_match_filter(imgc, img2, matches, vectors, params))
  bind(win, "<Control-r>", path->refresh(imgc, img2, meshset, match_ind))
  bind(win, "<Escape>", path->disable_inspection(imgc, img2))
  bind(win, "<Destroy>", path->disable_inspection(imgc, img2))
  # bind(win, "z", path->end_edit())
  # bind(win, "c", path->compare_filter(imgc, img2, matches, vectors))
  bind(win, "i", path->show_removed(imgc, img2, matches, vectors, params))
  bind(win, ",", path->decrease_distance_filter(imgc, img2, matches, vectors, params))
  bind(win, ".", path->increase_distance_filter(imgc, img2, matches, vectors, params))
  bind(win, "n", path->decrease_sigma_filter(imgc, img2, matches, vectors, params))
  bind(win, "m", path->increase_sigma_filter(imgc, img2, matches, vectors, params))
  bind(win, "b", path->adjust_sigma_filter(imgc, img2, matches, vectors, params))
  bind(win, "=", path->increase_vectors(imgc, img2, meshset, matches, vectors, params))
  bind(win, "-", path->decrease_vectors(imgc, img2, meshset, matches, vectors, params))
  bind(win, "p", path->switch_pre_to_post(imgc, img2, meshset, matches, vectors, params))
  bind(win, "s", path->save_inspection(meshset, match_ind))
  # bind(win, "<Return>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=true, flag=true, save=true))
  bind(win, "<Control-Shift-Right>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=true, flag=true))
  bind(win, "<Control-Right>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=true, flag=false))
  bind(win, "<Control-Shift-Left>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=false, flag=true))
  bind(win, "<Control-Left>", path->go_to_next_inspection(imgc, img2, meshset, match_ind; forward=false, flag=false))
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
  bind(win, "i", path->path)
  bind(win, ",", path->path)
  bind(win, ".", path->path)
  bind(win, "n", path->path)
  bind(win, "m", path->path)
  bind(win, "b", path->path)
  bind(win, "=", path->path)
  bind(win, "-", path->path)
  bind(win, "p", path->path)
  bind(win, "s", path->path)
  # bind(win, "<Return>", path->path)
  bind(win, "<Control-Shift-Right>", path->path)
  bind(win, "<Control-Right>", path->path)
  bind(win, "<Control-Shift-Left>", path->path)
  bind(win, "<Control-Left>", path->path)
  destroy(win)
end

function save_inspection(meshset, match_ind)
  println("meshset saved")
  set_reviewed!(meshset.matches[match_ind])
  save(meshset)
end

function go_to_next_inspection(imgc, img2, meshset, match_ind; forward=true, flag=true, save=false)
  if save
    save_inspection(meshset, match_ind)
  end

  if is_aligned(meshset)
    match_ind = get_name(meshset)
  end

  next_ms, next_match_ind = nothing, nothing
  if forward && flag
    println("Go to next flagged match")
    next_ms, next_match_ind = get_next_flagged_match(meshset, match_ind)
  elseif forward && !flag
    println("Go to next match")
    next_ms, next_match_ind = get_next_match(meshset, match_ind)
  elseif !forward && flag
    println("Go to previous flagged match")
    next_ms, next_match_ind = get_previous_flagged_match(meshset, match_ind)
  elseif !forward && !flag
    println("Go to previous match")
    next_ms, next_match_ind = get_previous_match(meshset, match_ind)
  end
  if next_ms == nothing && next_match_ind == nothing
    println("End of the inspection sequence")
  else
    disable_inspection(imgc, img2)
    inspect(next_ms, next_match_ind)
  end
end

function refresh(imgc, img2, meshset, match_ind)
  disable_inspection(imgc, img2)
  inspect(meshset, match_ind)
end

"""
Preliminary method to test filters in the GUI
"""
function compare_filter(imgc, img2, match, vectors, params)
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
  update_annotations(imgc, img2, match, vectors, params)
end

"""
Invert to just the bad matches
"""
function show_removed(imgc, img2, match, vectors, params)
  all_inds = Set(1:length(match.src_points))
  rejected_inds = get_rejected_indices(match)

  accepted_inds = collect(setdiff(all_inds, rejected_inds))

  clear_filters!(match)
  filter_manual!(match, accepted_inds)
  update_annotations(imgc, img2, match, vectors, params)
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
    detect_blockmatch_removal(imgc, img2, bm_win, matches, idx, vectors, params)
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

function remove_contained_points(imgc, img2, matches, vectors, bb, params)
  indices = get_filtered_indices(matches)
  pts = vectors[1:2, indices]
  within_bb = vcat((bb.xmin .<= pts[1,:] .<= bb.xmax) & (bb.ymin .<= pts[2,:] .<= bb.ymax)...)
  contained_indices = indices[within_bb]
  println("Brushtool removed ", length(contained_indices), " points")
  filter_manual!(matches, contained_indices)
  update_annotations(imgc, img2, matches, vectors, params)
end

function detect_blockmatch_removal(imgc::ImageView.ImageCanvas, 
                                      img2::ImageView.ImageSlice2d, 
                                      win, matches, idx, vectors, params)
  bind(win, "<Delete>", path->remove_blockmatch_from_patch_window(win, imgc, 
                                                  img2, matches, idx, vectors, params))
end

function remove_blockmatch_from_patch_window(win, imgc, img2, matches, idx, vectors, params)
  bind(win, "<Delete>", path->path)
  destroy(win)
  filter_manual!(matches, idx)
  update_annotations(imgc, img2, matches, vectors, params)
end

function remove_match(imgc, img2, x, y, matches, vectors, params, prox=0.0125)
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
    update_annotations(imgc, img2, matches, vectors, params)
  end
end

function undo_match_filter(imgc, img2, matches, vectors, params)
  println("Undo")
  idx = undo_filter!(matches)
  update_annotations(imgc, img2, matches, vectors, params)
end

function increase_distance_filter(imgc, img2, matches, vectors, params)
  params["dist"] += 10
  filter_match_distance(imgc, img2, matches, vectors, params)
end

function decrease_distance_filter(imgc, img2, matches, vectors, params)
  params["dist"] = max(params["dist"]-10, 0)
  filter_match_distance(imgc, img2, matches, vectors, params)
end

function filter_match_distance(imgc, img2, matches, vectors, params)
  # hack to test if a match_distance filter was just implemented
  dist = params["dist"]
  if length(matches.filters) > 0
    if matches.filters[end]["type"] == "norm"
      undo_filter!(matches)
    end
  end
  println("Distance filter @ ", dist)
  filter = (:get_properties, >, dist, "norm")
  filter!(matches, filter...)
  update_annotations(imgc, img2, matches, vectors, params)
end

function adjust_sigma_filter(imgc, img2, matches, vectors, params)
  println("Enter sigma filter value:")
  val = chomp(readline())
  val = parse(Int, val)
  params["sigma"] = val
  filter_match_sigma(imgc, img2, matches, vectors, params)
end

function increase_sigma_filter(imgc, img2, matches, vectors, params)
  params["sigma"] += 1
  filter_match_sigma(imgc, img2, matches, vectors, params)
end

function decrease_sigma_filter(imgc, img2, matches, vectors, params)
  params["sigma"] = max(params["sigma"]-1, 0)
  filter_match_sigma(imgc, img2, matches, vectors, params)
end

function filter_match_sigma(imgc, img2, matches, vectors, params)
  sigma = params["sigma"]
  if length(matches.filters) > 0
    if matches.filters[end]["type"] == 0.5
      undo_filter!(matches)
    end
  end
  println("Sigma filter @ ", sigma)
  filter = (:get_properties, >, sigma, 0.5)
  filter!(matches, filter...)
  update_annotations(imgc, img2, matches, vectors, params)
end

function toggle_flag(imgc, img2, match)
  if is_flagged(match)
    println("Unflag match!")
    unflag!(match)
  else
    println("Flag match!")
    flag!(match)
  end
end

function increase_vectors(imgc, img2, meshset, matches, vectors, params)
  params["vector_scale"] *= 1.1
  update_annotations(imgc, img2, matches, vectors, params)
end

function decrease_vectors(imgc, img2, meshset, matches, vectors, params)
  params["vector_scale"] = max(params["vector_scale"]*0.9, 0.1)
  update_annotations(imgc, img2, matches, vectors, params)
end

function switch_pre_to_post(imgc, img2, meshset, matches, vectors, params)
  params["post_matches"] = !params["post_matches"]
  println(params["post_matches"] ? "Using post matches" : "Using pre matches" )
  k = params["match_index"]
  update_annotations(imgc, img2, matches, vectors, params)
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
  attr = get_properties(match, property_name)
  attr_filtered = get_filtered_properties(match, property_name)
  if length(attr) > 1
    max_bin = min(maximum(attr), 10000)
    min_bin = minimum(attr)
    bins = [linspace(min_bin, max_bin, nbins)]
    if filtered
      attr = attr_filtered
      color="#009900"
    end
    if length(attr) > 1
      p = plt[:hist](attr, bins=bins, color=color)
      grid("on")
      title(property_name)
      return p
    end
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
  inc = max(round(Int64, 2*sr/40), 10)
  dv_bounds = vcat([[i -sr] for i=-sr:inc:sr]...,
                    [[i sr] for i=-sr:inc:sr]...,
                    [[sr i] for i=-sr:inc:sr]...,
                    [[-sr i] for i=-sr:inc:sr]...)
  if length(dv) > 1
    p = plt[:scatter](dv[2,:], -dv[1,:], color=color, alpha=0.5)
    p = plt[:scatter](dv_bounds[:,1], dv_bounds[:,2], color="#0f0f0f", marker="+", alpha=0.2)
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
  view_property_spatial_scatter(match, "norm"; filtered=false, factor=1)
  view_property_spatial_scatter(match, "norm"; filtered=true, factor=1)
  subplot(236)
  view_property_histogram(match, 0.5; filtered=false, nbins=20)
  view_property_histogram(match, 0.5; filtered=true, nbins=20)

  src_index = get_src_index(match)
  dst_index = get_dst_index(match)
  annotate("$src_index-$dst_index",
    xy=[1;1],
    xycoords="figure fraction",
    xytext=[-10,-10],
    textcoords="offset points",
    ha="right",
    va="top",
    fontsize=16)
  annotate(string("Flags:\n", join(keys(match.properties["review"]["flags"]), "\n")),
    xy=[0;1],
    xycoords="figure fraction",
    xytext=[10,-10],
    textcoords="offset points",
    ha="left",
    va="top",
    fontsize=16)
  fig[:canvas][:draw]()
end
