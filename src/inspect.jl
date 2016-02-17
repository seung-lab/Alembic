"""
The only function called by tracers to inspect montage points
"""
function inspect_montages(meshset_ind, match_ind)
  indrange = get_index_range((1,1,-2,-2), (8,173,-2,-2))
  meshset = load(indrange[meshset_ind])
  println("\n", meshset_ind, ": ", indrange[meshset_ind], " @ ", match_ind, " / ", length(meshset.matches))
  imgc, img2, matches, vectors, params = inspect_matches(meshset, match_ind, "seam");
  enable_inspection(imgc, img2, meshset, matches, vectors, params, "montage", (meshset_ind, match_ind))
end

"""
The only function called by tracers to inspect prealignment points
"""
function inspect_prealignments(meshset_ind)
  match_ind = 1
  index_pairs = collect(get_sequential_index_pairs((1,1,-2,-2), (2,149,-2,-2)))
  indexA, indexB = index_pairs[meshset_ind]
  meshset = load(indexB, indexA)
  println("\n", meshset_ind, ": ", (indexB, indexA), " @ ", match_ind, " / ", length(meshset.matches))
  imgc, img2, matches, vectors, params = inspect_matches(meshset, match_ind, "thumb")
  enable_inspection(imgc, img2, meshset, matches, vectors, params, "prealignment", (meshset_ind, match_ind))
end

"""
The only function called by tracers to inspect alignment points
"""
function inspect_alignments(meshset, match_ind)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[count_meshes(meshset)].index
  name = string(join(firstindex[1:2], ","),  "-", join(lastindex[1:2], ","),"_aligned.jls")
  match = meshset.matches[match_ind]
  src_dst = string(join(match.src_index[1:2], ","), "-", join(match.dst_index[1:2], ","))
  println("\n", name, ": ", src_dst, " @ ", match_ind, " / ", length(meshset.matches))
  imgc, img2, matches, vectors, params = inspect_matches(meshset, match_ind, "thumb_imfuse");
  enable_inspection(imgc, img2, meshset, matches, vectors, params, "alignment", (1, match_ind));
end



function show_blockmatch(match, ind, params)
  src_patch, src_pt, dst_patch, dst_pt, xc, offset = get_correspondence_patches(match, ind)
  block_r = params["block_r"]
  search_r = params["search_r"]
  N=size(xc, 1)
  M=size(xc, 2)
  if !contains(gethostname(), "seunglab")
    surf([i for i=1:N, j=1:M], [j for i=1:N, j=1:M], xc, cmap=get_cmap("hot"), 
                            rstride=10, cstride=10, linewidth=0, antialiased=false)
  end
  println("offset: ", offset)
  xc_image = xcorr2Image(xc)
  # xc_image = padimage(xc_image, block_r, block_r, block_r, block_r, 1)
  hot = create_hot_colormap()
  xc_color = apply_colormap(xc_image, hot)
  fused_img, _ = imfuse(dst_patch, [0,0], src_patch, offset)

  src = convert(Array{UInt32,2}, src_patch).<< 8
  dst = convert(Array{UInt32,2}, dst_patch).<< 16

  cgrid = canvasgrid(2,2; pad=10)
  opts = Dict(:pixelspacing => [1,1])

  imgc, img2 = view(cgrid[1,1], src; opts...)
  imgc, img2 = view(cgrid[2,1], dst; opts...)
  imgc, img2 = view(cgrid[2,2], fused_img; opts...)
  imgc, img2 = view(cgrid[1,2], xc_color'; opts...)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  return win
end

"""
Extend ImageView to inspect & remove matches
"""
function enable_inspection(imgc::ImageView.ImageCanvas, 
                                img2::ImageView.ImageSlice2d, 
                                meshset, matches, vectors, params, stage, calls)
  println("Enable inspection")
  c = canvas(imgc)
  win = Tk.toplevel(c)
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
  bind(win, ",", path->decrease_distance_filter(imgc, img2, matches, vectors, params))
  bind(win, ".", path->increase_distance_filter(imgc, img2, matches, vectors, params))
  bind(win, "=", path->increase_vectors(imgc, img2, meshset, matches, vectors, params))
  bind(win, "-", path->decrease_vectors(imgc, img2, meshset, matches, vectors, params))
  bind(win, "p", path->switch_pre_to_post(imgc, img2, meshset, matches, vectors, params))
  bind(win, "s", path->save_inspection(meshset))
  bind(win, "<Return>", path->go_to_next_inspection(imgc, img2, meshset, stage, calls))
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
  bind(win, ",", path->path)
  bind(win, ".", path->path)
  bind(win, "=", path->path)
  bind(win, "-", path->path)
  bind(win, "p", path->path)
  bind(win, "s", path->path)
  bind(win, "<Return>", path->path)
  destroy(win)
end

function save_inspection(meshset)
  println("meshset saved")
  # save(meshset)
end

function go_to_next_inspection(imgc, img2, meshset, stage, calls)
  disable_inspection(imgc, img2)

  meshset_ind, match_ind = calls
  match_ind += 1
  if match_ind > length(meshset.matches)
    meshset_ind += 1
    match_ind = 1
  end

  if stage == "alignment"
    inspect_alignments(meshset, match_ind)
  elseif stage == "prealignment"
    inspect_prealignments(meshset_ind)
  elseif stage == "montage"
    inspect_montages(meshset_ind, match_ind)
  else
    inspect_montages(meshset_ind, match_ind)
  end
end

"""
Convention: mask is FALSE if point is to be REMOVED
"""
function update_annotations(imgc, img2, matches, params)
  mask = get_filtered_indices(matches)
  for an in values(imgc.annotations)
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
    bm_win = show_blockmatch(matches, idx, params)
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

function flag_inspection(imgc, img2, matches)
  println("FLAG MATCHES! (not implemented, yet)")
  # flag!(matches)
end

function increase_vectors(imgc, img2, meshset, matches, vectors, params)
  params["vector_scale"] += 1
  k = params["match_index"]
  vectors = make_vectors(meshset, k, params)
  update_annotations(imgc, img2, matches, vectors)
end

function decrease_vectors(imgc, img2, meshset, matches, vectors, params)
  params["vector_scale"] = max(params["vector_scale"]-1, 1)
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

"""
Display matches index k as overlayed images with points
"""
function inspect_matches(meshset, k, prefix="review")
  matches = meshset.matches[k]
  indexA = matches.src_index
  indexB = matches.dst_index

  path = get_review_filename(prefix, indexB, indexA)
  if !isfile(path)
    path = get_review_filename(prefix, indexA, indexB)
  end
  img = h5read(path, "img")
  offset = h5read(path, "offset")
  scale = h5read(path, "scale")

  # Add border to the image for vectors that extend
  # pad = 200
  # img = zeros(UInt32, size(img,1)+pad*2, size(img,2)+pad*2)
  # img[pad:end-pad, pad:end-pad] = img_orig
  # img = vcat(img, ones(UInt32, 400, size(img, 2)))

  params = meshset.properties["params"]["match"]
  params["offset"] = offset
  params["scale"] = scale
  params["match_index"] = k
  params["vector_scale"] = 4
  params["post_matches"] = false
  params["dist"] = 90

  imgc, img2 = view(img, pixelspacing=[1,1])
  vectors = make_vectors(meshset, k, params)
  show_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1))
  update_annotations(imgc, img2, matches, vectors)
  return imgc, img2, matches, vectors, params
end

<<<<<<< HEAD
function make_vectors(meshset, k, params)
  scale = params["scale"]
  offset = params["offset"]
  factor = params["vector_scale"]
  println("Vector scale: ", factor)
  src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences(meshset, k)
  if params["post_matches"]
    src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences_post(meshset, k)
  end
=======
  # src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences(meshset, k)
  src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences_post(meshset, k)
  mask = !indices_to_mask(filtered_inds, length(src_nodes))
>>>>>>> origin/master
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vecs = offset_matches(vectorsA, vectorsB, offset)
  vectors = [hcat(vecs[1]...); hcat(vecs[2]...)]
  return change_vector_lengths([hcat(vecs[1]...); hcat(vecs[2]...)], factor)
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
  end
  return meshset
end

function update_montage_meshsets(waferA, secA, waferB, secB)
  indices = get_index_range((waferA,secA,-2,-2), (waferB,secB,-2,-2))
  logs = compile_review_logs("montage")
  logs = get_most_recent_logs(logs)
  for index in indices
    meshset = load(index)
    meshset = get_meshset_with_edits(meshset, index, logs)
    save(meshset)
    solve!(meshset, method="elastic")
    save(meshset)
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

"""
Write points to remove in a text file
"""
function store_points(path, meshset, k, indices_to_remove, inspection_flagged, username, comment)
  ts = Dates.format(now(), "yymmddHHMMSS")
  src_index = meshset.matches[k].src_index
  dst_index = meshset.matches[k].dst_index
  if length(indices_to_remove) == 0
    indices_to_remove = [0]
  end
  pts_line = [ts, username, src_index, dst_index, k, join(indices_to_remove, ","), comment]'
  if !isfile(path)
    f = open(path, "w")
    close(f)
    pts = pts_line
  else  
    pts = readdlm(path)
    pts = vcat(pts, pts_line)
  end
  pts = pts[sortperm(pts[:, 5]), :]
  println("Saving pts_to_remove:\n", path)
  writedlm(path, pts)
end
