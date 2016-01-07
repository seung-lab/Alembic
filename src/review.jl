"""
Remove matches from meshset corresponding to the indices provided.

Args:

* indices_to_remove: array of indices of match points to remove
* match_index: the index of the matches object in the meshset
* meshset: the meshset containing all the matches

Returns:

* (updated meshset)

  remove_matches_from_meshset!(indices_to_remove, match_index, meshset)
"""
function remove_matches_from_meshset!(meshset, indices_to_remove, match_index=1)
  println("Removing ", length(indices_to_remove), " points from ", match_index)
  println("before: ", length(meshset.matches[match_index].dst_points))
  if length(indices_to_remove) > 0
    no_pts_removed = length(indices_to_remove)
    matches = meshset.matches[match_index]
    flag = trues(matches.n)
    flag[indices_to_remove] = false
    matches.src_points_indices = matches.src_points_indices[flag]
    matches.n -= no_pts_removed
    matches.dst_points = matches.dst_points[flag]
    matches.dst_triangles = matches.dst_triangles[flag]
    matches.dst_weights = matches.dst_weights[flag]
    matches.disp_vectors = matches.disp_vectors[flag]
    meshset.m -= no_pts_removed
    meshset.m_e -= no_pts_removed
  end
  println("after: ", length(meshset.matches[match_index].dst_points))
end

"""
Combine tracer point inspection logs into one table
"""
function compile_tracer_point_inspection_logs(meshset)
  tracers = ["hmcgowan", "bsilverman", "tmacrina_cleanup"]
  logs = []
  for tracer in tracers
    path = get_storage_path(meshset, tracer)
    if isfile(path)
      push!(logs, readdlm(path))
    end
  end
  log = vcat(logs...)
  return log[sortperm(log[:, 1]), :]
end

"""
Return indices of matches that have not been reviewed by any tracer
"""
function list_missing_inspections(meshset)
  log = compile_tracer_point_inspection_logs(meshset)
  match_nums = collect(1:length(meshset.matches))
  inspections = unique(log[:,5])
  match_nums[inspections] = 0
  return match_nums[match_nums .!= 0]
end

function update_meshset_with_edits!(meshset, log)
  for i in 1:size(log,1)
    match_index = log[i,5]
    indices_to_remove = [log[i,6]]
    if typeof(indices_to_remove[1]) != Int64
      indices_to_remove = readdlm(IOBuffer(log[i,6]), ',', Int)
    end
    if indices_to_remove[1] != 0
      remove_matches_from_meshset!(meshset, indices_to_remove, match_index)
    end
  end
end

function crop_meshset!(meshset, a, b)
  meshset.indices = meshset.indices[a:b]
  meshset.meshes = meshset.meshes[a:b]
  meshset.N = length(meshset.meshes)
  meshset.nodes_indices = meshset.nodes_indices[a:b]
  meshset.nodes_indices -= meshset.nodes_indices[1]
  matches_mask = map(x -> x>=(a,a) && x<=(b,b), meshset.matches_pairs)
  matches_pairs = meshset.matches_pairs[matches_mask]
  meshset.matches_pairs = map(x->(x[1]-a+1, x[2]-a+1), matches_pairs)
  meshset.matches = meshset.matches[matches_mask]
  meshset.M = length(meshset.matches)

  meshset.n = 0
  meshset.m = 0
  meshset.m_i = 0
  meshset.m_e = 0
  for mesh in meshset.meshes
    meshset.n += mesh.n
    meshset.m += mesh.m
    meshset.m_i += mesh.m

  end
  for matches in meshset.matches
    meshset.m += matches.n
    meshset.m_e += matches.n
  end
end

function remove_match_dict_from_meshset!(meshset, match_dict)
  for (key, val) in match_dict
    remove_matches_from_meshset!(meshset, val, key)
  end
  @time solve_meshset!(meshset)
  save(meshset)
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

function find_pt_idx(vectors, pt)
  i = 1
  while i <= size(vectors, 2)
    if vectors[1:2,i] == pt
      return i
    end
    i += 1
  end
  return 0
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

function updatexylabel(xypos, imgc, img, x, y)
    w = width(imgc.c)
    xu,yu = ImageView.device_to_user(Tk.getgc(imgc.c), x, y)
    # Image-coordinates
    xi,yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
    if ImageView.isinside(imgc.canvasbb, x, y)
      if 1 <= xi <= size(img,2) && 1 <= yi <= size(img,1)
        val = img[yi,xi]
        str = "$yi, $xi: $val"
        if length(str)*10>w
            ImageView.set_value(xypos, "$yi, $xi")
        else
            ImageView.set_value(xypos, str)
        end
      end
    else
        ImageView.set_value(xypos, "$yi, $xi")
    end
end

function xcorr2Image(xc)
  b = xc' / maximum(xc)
  println("max r-value: ", maximum(xc))
  b[b .> 1] = 1
  b[b .< 0] = 0
  # b = b.*b / maximum(b.*b) * 254 + 1
  b = b * 254 + 1
  # b = (b + 1) ./ 2 * 254 + 1
  b[isnan(b)] = 0
  return round(UInt8, b)
  # return round(UInt8, (b + 1) ./ 2 * 254 + 1)
end

function draw_box(img::Array{UInt32}, bounds)
  srf = create_drawing(img)
  ctx = create_contex(srf)
  draw_rect(ctx, bounds, (1,1,1), 3)
  return get_drawing(srf)
end

function display_blockmatch(src_patch, src_pt, dst_patch, dst_pt, xc, params)
  block_r = params["block_size"]
  search_r = params["search_r"]
  N=size(xc, 1)
  M=size(xc, 2)
  surf([i for i=1:N, j=1:M], [j for i=1:N, j=1:M], xc, cmap=get_cmap("hot"), 
                          rstride=10, cstride=10, linewidth=0, antialiased=false)
  xc_image = xcorr2Image(xc)
  xc_image = padimage(xc_image, block_r, block_r, block_r, block_r, 1)
  hot = create_hot_colormap()
  xc_color = apply_colormap(xc_image, hot)
  xc = padimage(xc', block_r, block_r, block_r, block_r, 1)

  offset = round(Int64, (dst_pt-block_r) - (src_pt-block_r-search_r))
  println(offset)
  # reverse_offset = collect(size(dst_img_orig)) - (collect(size(src_img)) + offset)
  # src_padded = padimage(src_img, reverse(offset)..., reverse(reverse_offset)...)
  fused_img, _ = imfuse(dst_img_orig, [0,0], src_img, offset)
  # view(reinterpret(UFixed8, src_img'), pixelspacing=[1,1])

  bounds = (search_r, search_r, size(src_img)...)
  src = draw_box(convert(Array{UInt32,2}, src_img_full).<< 8, bounds)
  bounds = (offset..., size(src_img)...)
  dst = draw_box(convert(Array{UInt32,2}, dst_img_orig).<< 16, bounds)

  cgrid = canvasgrid(2,2; pad=10)
  opts = Dict(:pixelspacing => [1,1])

  imgc, img2 = view(cgrid[1,1], src'; opts...)
  imgc, img2 = view(cgrid[2,1], dst'; opts...)
  imgc, img2 = view(cgrid[2,2], fused_img; opts...)
  imgc, img2 = view(cgrid[1,2], xc_color'; opts...)
end

function display_blockmatch(params, src_point, dst_point)
  src_index = params["src_index"]
  dst_index = params["dst_index"]
  block_r = params["block_size"]
  search_r = params["search_r"]
  
  src_pt = src_point - params["src_offset"]
  dst_pt = dst_point - params["dst_offset"]
  dst_pt_orig = src_point - params["dst_offset"]

  src_size = params["src_size"]
  dst_size = params["dst_size"]  

  src_bnds = sliceimg(src_pt, block_r, src_size)
  dst_bnds = sliceimg(dst_pt, block_r+search_r, dst_size)
  src_bnds_full = sliceimg(src_pt, block_r+search_r, src_size)
  dst_bnds_orig = sliceimg(dst_pt_orig, block_r+search_r, dst_size)

  src_slice = colon(src_bnds[1:2]...), colon(src_bnds[3:4]...)
  dst_slice = colon(dst_bnds[1:2]...), colon(dst_bnds[3:4]...)
  src_slice_full = colon(src_bnds_full[1:2]...), colon(src_bnds_full[3:4]...)
  dst_slice_orig = colon(dst_bnds_orig[1:2]...), colon(dst_bnds_orig[3:4]...)
  
  src_img = get_h5_slice(get_h5_path(src_index), src_slice)
  dst_img_new = get_h5_slice(get_h5_path(dst_index), dst_slice)
  src_img_full = get_h5_slice(get_h5_path(src_index), src_slice_full)  
  dst_img_orig = get_h5_slice(get_h5_path(dst_index), dst_slice_orig)  

  xc = normxcorr2(src_img, dst_img_orig)
  N=size(xc, 1)
  M=size(xc, 2)
  surf([i for i=1:N, j=1:M], [j for i=1:N, j=1:M], xc, cmap=get_cmap("hot"), 
                          rstride=10, cstride=10, linewidth=0, antialiased=false)
  xc_image = xcorr2Image(xc)
  xc_image = padimage(xc_image, block_r, block_r, block_r, block_r, 1)
  hot = create_hot_colormap()
  xc_color = apply_colormap(xc_image, hot)
  xc = padimage(xc', block_r, block_r, block_r, block_r, 1)

  offset = round(Int64, (dst_point-block_r) - (src_point-block_r-search_r))
  println(offset)
  # reverse_offset = collect(size(dst_img_orig)) - (collect(size(src_img)) + offset)
  # src_padded = padimage(src_img, reverse(offset)..., reverse(reverse_offset)...)
  fused_img, _ = imfuse(dst_img_orig, [0,0], src_img, offset)
  # view(reinterpret(UFixed8, src_img'), pixelspacing=[1,1])

  bounds = (search_r, search_r, size(src_img)...)
  src = draw_box(convert(Array{UInt32,2}, src_img_full).<< 8, bounds)
  bounds = (offset..., size(src_img)...)
  dst = draw_box(convert(Array{UInt32,2}, dst_img_orig).<< 16, bounds)

  cgrid = canvasgrid(2,2; pad=10)
  opts = Dict(:pixelspacing => [1,1])

  imgc, img2 = view(cgrid[1,1], src'; opts...)
  imgc, img2 = view(cgrid[2,1], dst'; opts...)
  imgc, img2 = view(cgrid[2,2], fused_img; opts...)
  imgc, img2 = view(cgrid[1,2], xc_color'; opts...)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  fnotify = ImageView.Frame(win)
  lastrow = 1
  ImageView.grid(fnotify, lastrow+=1, 1, sticky="ew")
  xypos = ImageView.Label(fnotify)
  imgc.handles[:pointerlabel] = xypos
  ImageView.grid(xypos, 1, 1, sticky="ne")
  ImageView.set_visible(win, true)
  c.mouse.motion = (path,x,y)-> updatexylabel(xypos, imgc, xc, x, y)

  return win
end

function edit_matches(imgc, img2, vectors, vectors_t, params)
    e = Condition()

    indices_to_remove = Array{Integer,1}()
    lines = Void
    original_lines = Void
    for annotation in values(imgc.annotations)
      if :lines in fieldnames(annotation.data)
        lines = annotation.data.lines
        original_lines = copy(lines)
      end
    end
    mask = ones(Bool, size(original_lines,2))

    offset = params["thumb_offset"]
    scale = params["scale"]

    c = canvas(imgc)
    win = Tk.toplevel(c)
    bind(c, "<Button-3>", (c, x, y)->inspect_match(parse(Int, x), parse(Int, y)))
    bind(c, "<Control-Button-3>", (c, x, y)->remove_match(parse(Int, x), parse(Int, y)))
    bind(win, "<Delete>", path->path)
    bind(win, "<Escape>", path->end_edit())
    bind(win, "<Destroy>", path->end_edit())

    inspect_window = Void
    inspected_index = Void
    inspected_index_removed = true

    function remove_match(x, y)
        xu,yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
        xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
        # println(xi, ", ", yi)
        limit = (img2.zoombb.xmax - img2.zoombb.xmin) * 0.0125 # 100/8000
        annidx = find_idx_of_nearest_pt(lines[1:2,:], [xi, yi], limit)
        if annidx > 0
            pt = lines[1:2,annidx]
            idx = find_pt_idx(original_lines[1:2,:], pt)
            println(idx, ": ", original_lines[1:2, idx])
            mask[idx] = false
            update_annotations(imgc, img2, original_lines, mask)
            push!(indices_to_remove, idx)
        end
    end

    function inspect_match(x, y)
        xu,yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
        xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
        # println(xi, ", ", yi)
        limit = (img2.zoombb.xmax - img2.zoombb.xmin) * 0.0125 # 100/8000
        annidx = find_idx_of_nearest_pt(lines[1:2,:], [xi, yi], limit)
        if annidx > 0
          annpt = ([lines[2,annidx], lines[1,annidx]] + offset) / scale 
          println(annpt)
          idx = find_idx_of_nearest_pt(vectors_t[1:2,:], annpt, 10)
          if idx > 0
            ptA = vectors[1:2,idx] # - params["src_offset"]
            ptB = vectors[3:4,idx] # - params["dst_offset"]
            println(ptA, ", ", ptB)
            # keep_point = display_blockmatch(params, ptA, ptB)
            # if !keep_point
            #   mask[idx] = false
            #   update_annotations(imgc, img2, original_lines, mask)
            #   push!(indices_to_remove, idx)
            # end
            inspect_window = display_blockmatch(params, ptA, ptB)
            inspected_index_removed = false
            inspected_index = idx
          end
        end
    end

    function remove_inspected_match()
      if !inspected_index_removed
        if inspected_index > 0
          println(inspected_index, ": ", original_lines[1:2, inspected_index])
          mask[inspected_index] = false
          update_annotations(imgc, img2, original_lines, mask)
          push!(indices_to_remove, inspected_index)
          inspected_index_removed = true
          destroy(inspect_window)
          PyPlot.close()
        end
      end
    end

    function end_edit()
        println("End edit\n")
        notify(e)
        bind(c, "<Button-3>", path->path)
        bind(c, "<Control-Button-3>", path->path)
        bind(win, "<Delete>", path->path)
        bind(win, "<Destroy>", path->path)
        bind(win, "<Escape>", path->path)
    end

    println("1) Right click to inspect correspondences\n",
            "2) Ctrl + right click to remove correspondences\n",
            "3) Exit window or press Escape")
    wait(e)

    return indices_to_remove
end

"""
Maintain vector start point, but adjust end point for more prominent visual
"""
function change_vector_lengths(vectors, k)
  v = [vectors[2,:]; 
        vectors[1,:]; 
        (vectors[4,:]-vectors[2,:])*k + vectors[2,:]; 
        (vectors[3,:]-vectors[1,:])*k + vectors[1,:]]
  return v
end

function filter_distance(imgc, img2, lines, disp, dist)
  mask = update_annotations(imgc, img2, lines, disp .< dist)
  return mask, "disp.<$dist"
end

"""
Convention: mask is FALSE if point is to be REMOVED
"""
function update_annotations(imgc, img2, lines, mask)
  for an in values(imgc.annotations)
    if :pts in fieldnames(an.data)
      an.data.pts = lines[1:2, mask]
    elseif :lines in fieldnames(an.data)
      an.data.lines = lines[:, mask]
    end
  end
  println("Removing ", sum(!mask), " matches from GUI")
  ImageView.redraw(imgc)
  return mask
end

function review_matches(meshset, indexA, indexB)
  k = find_matches_index(meshset, indexA, indexB)
  println(indexA, indexB, ": ", k)
  assert(k > 0)
  return review_matches(meshset, k)
end

"""
Display matches index k as overlayed images with points
"""
function review_matches(meshset, k)
  matches = meshset.matches[k]
  indexA = (matches.src_index[1:2]..., -4, -4)
  indexB = (matches.dst_index[1:2]..., -4, -4)

  path = get_outline_filename("thumb_imfuse", indexB, indexA)
  img = h5read(path, "img")
  offset = h5read(path, "offset")
  scale = h5read(path, "scale")

  img = vcat(img, ones(UInt32, 400, size(img, 2)))

  params = meshset.params
  params["scale"] = scale
  params["thumb_offset"] = offset
  params["src_index"] = matches.src_index
  params["dst_index"] = matches.dst_index
  params["src_offset"] = meshset.meshes[find_index(meshset, matches.src_index)].disp
  params["dst_offset"] = meshset.meshes[find_index(meshset, matches.dst_index)].disp
  params["src_size"] = get_image_size(matches.src_index)
  params["dst_size"] = get_image_size(matches.dst_index)

  src_nodes, dst_nodes = get_matched_points(meshset, k)
  vectors = [hcat(src_nodes...); hcat(dst_nodes...)]
  src_nodes, dst_nodes = get_matched_points_t(meshset, k)
  vectors_t = [hcat(src_nodes...); hcat(dst_nodes...)]
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vecs = offset_matches(vectorsA, vectorsB, offset)

  imview = view(img, pixelspacing=[1,1])
  big_vecs = change_vector_lengths([hcat(vecs[1]...); hcat(vecs[2]...)], 10)
  an_pts, an_vectors = show_vectors(imview..., big_vecs, RGB(0,0,1), RGB(1,0,1))
  return imview, vectors, vectors_t, copy(an_vectors.ann.data.lines), params
end

"""
Get indices of the false entries in a binary array
"""
function mask_to_indices(mask)
  return eachindex(mask)[!mask]
end

"""
Store indices of the false entries from a binary mask as matches to remove
"""
function store_mask(path, meshset, k, mask, username, comment)
  indices_to_remove = mask_to_indices(mask)
  store_points(path, meshset, k, indices_to_remove, username, comment)
end

function get_inspection_groundtruth_path()
  # return "/usr/people/tmacrina/seungmount/Omni/alignment/training/1,2-1,16_aligned_EDITED_tmacrina_20151113162553.txt"
  fn = "1,2-1,16_aligned_EDITED_tmacrina_20151113162553.txt"
  return joinpath(inspection_storage_path, fn)
end

function dict_of_inspections(path, meshset)
  d = Dict()
  pts = readdlm(path)

  # remove_matches_from_meshset!(meshset, val, key)
  for i in 1:size(pts,1)
    match_index = pts[i,5]
    if !(match_index in keys(d))
      # d[match_index] = Array{Array{Int64, 1}, 1}()
      accepted = collect(eachindex(meshset.matches[match_index].dst_points))
      rejected = Array{Int64, 1}()
      # push!(d[match_index], accepted)
      # push!(d[match_index], rejected)
      d[match_index] = accepted, rejected
    end
  end

  for i in 1:size(pts,1)
    match_index = pts[i,5]
    indices_to_remove = readdlm(IOBuffer(pts[i,6]), ',', Int)
    flag = trues(length(d[match_index][1]))
    flag[indices_to_remove] = false
    accepted, rejected = d[match_index]
    rejected = vcat(rejected, d[match_index][1][!flag])
    accepted = d[match_index][1][flag]
    d[match_index] = accepted, rejected
  end
  return d
end

function dict_of_inspections(path)
  d = Dict()
  pts = readdlm(path)
  for i in 1:size(pts,1)
    match_index = pts[i,5]
    indices_to_remove = Set(readdlm(IOBuffer(pts[i,6]), ',', Int))
    if match_index in keys(d)
      d[match_index] = union(d[match_index], indices_to_remove)
    else
      d[match_index] = indices_to_remove
    end
  end
  return d
end

function compare_inspections(meshset, pathA, pathB=get_inspection_groundtruth_path())
  dC = Dict()
  dA = dict_of_inspections(pathA, meshset)
  dB = dict_of_inspections(pathB, meshset)
  sections = intersect(Set(keys(dB)), Set(keys(dA)))
  for k in sections
    acceptedA, rejectedA = Set(dA[k][1]), Set(dA[k][2])
    acceptedB, rejectedB = Set(dB[k][1]), Set(dB[k][2])
    # [TP in A, TN in A, FP in A, FN in A] # TN: match properly removed
    dC[k] = [intersect(acceptedA, acceptedB), 
              intersect(rejectedA, rejectedB), 
              intersect(acceptedA, rejectedB), 
              intersect(rejectedA, acceptedB)]
  end
  return dC
end

function print_inspection_report(path, meshset)
  dC = compare_inspections(meshset, path)
  report = ["k" "TP" "TN" "FP" "FN"]
  for k in sort(collect(keys(dC)))
    tp, tn, fn, fp = map(length, dC[k])
    report = vcat(report, [k tp tn fn fp])
  end
  report = vcat(report, sum(report[2:end,:],1))
  path = string(path[1:end-4], "_report.txt")
  println("Saving report:\n", path)
  writedlm(path, report)
  return report
end

function display_discrepant_match(params, vectors, idx)
  ptA = vectors[1:2,idx] # - params["src_offset"]
  ptB = vectors[3:4,idx] # - params["dst_offset"]
  pt_display = ptA*params["scale"]-params["thumb_offset"]
  println(idx, ": ", reverse(pt_display))
  display_blockmatch(params, ptA, ptB)
end

"""
Write points to remove in a text file
"""
function store_points(path, meshset, k, indices_to_remove, username, comment)
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

function get_storage_path(meshset, ts="")
  if ts == ""
    ts = Dates.format(now(), "yyyymmddHHMMSS")
  end
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[end].index
  fn = string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","),
                "_aligned.txt")
  fn = update_filename(fn, ts)
  return joinpath(INSPECTION_DIR, fn)
end

function review_meshset_matches(meshset, ts="", start=1, finish=0)
  if ts == ""
    ts = Dates.format(now(), "yyyymmddHHMMSS")
  end
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[end].index
  storage_path = update_filename(get_name(firstindex, lastindex, "txt"), ts)
  meshset_path = update_filename(get_name(firstindex, lastindex, "jld"), ts)

  if finish == 0
    finish = length(meshset.matches)
  end

  for k in start:finish
    pts_to_remove = review_matches(meshset, k)
    store_points(storage_path, k, pts_to_remove)
  end
end

"""
# Pipeline
imview, vectors, disp, k = prepare_review_matches(meshset, 18);
mask = update_annotations(imview..., vectors, disp .< 70);
store_mask(path, k, mask, "disp.<70");
indices_to_remove = edit_matches(imview..., vectors[:,mask]);
store_points(path, k, indices_to_remove, "manual_after_disp");
"""
function prepare_review_matches(meshset, k)
  # r_values = calculate_r_values(meshset, k)
  disp = calculate_displacements(meshset, k)
  imview, vectors, vectors_t, lines, params = review_matches(meshset, k)
  return imview, vectors, vectors_t, lines, disp, k, params
end

"""
Display matches index k as overlayed images with points & image zoomed to bounds
"""
function prepare_review_matches(meshset, k, bounds)
  imview, v, vt, l, d, k, p = prepare_review_matches(meshset, k)
  b = [bounds[2][1], bounds[2][end], bounds[1][1], bounds[1][end]]
  bb = Graphics.BoundingBox(b*params["scale"]...)
  ImageView.zoombb(imview..., bb)
  return imview, v, vt, l, d, k, p
end

"""
Append 'EDITED' with a timestamp to update filename
"""
function update_filename(fn, ts="")
  if ts == ""
    ts = Dates.format(now(), "yyyymmddHHMMSS")
  end
  return string(fn[1:end-4], "_EDITED_", ts, fn[end-3:end])
end

function calculate_r_values(meshset, k)
  r_values = []
  matches = meshset.matches[k]
  src_index = matches.src_index
  dst_index = matches.dst_index
  src_mesh = meshset.meshes[find_index(meshset, src_index)]
  dst_mesh = meshset.meshes[find_index(meshset, dst_index)]
  src_offset = src_mesh.disp
  dst_offset = dst_mesh.disp

  block_r = meshset.params["block_size"]
  src_img = get_h5_image(get_h5_path(src_index))
  dst_img = get_h5_image(get_h5_path(dst_index))
  src_nodes, dst_nodes = get_matched_points(meshset, k)
  src_pts = [x - src_offset for x in src_nodes]
  dst_pts = [x - dst_offset for x in dst_nodes]
  for (src_pt, dst_pt) in zip(src_pts, dst_pts)
    src_bnds = sliceimg(src_pt, block_r, size(src_img))
    dst_bnds = sliceimg(dst_pt, block_r, size(dst_img))
    src = src_img[colon(src_bnds[1:2]...), colon(src_bnds[3:4]...)]
    dst = dst_img[colon(dst_bnds[1:2]...), colon(dst_bnds[3:4]...)]
    # println(size(src), " ", size(dst))
    push!(r_values, normxcorr2(src, dst)[1])
  end
  return r_values
end

function calculate_displacements(meshset, k)
  src_nodes, dst_nodes = get_matched_points_t(meshset, k)
  return map(norm, src_nodes-dst_nodes)
end

"""
Evaluating filters
"""
function removed_to_truth(n, removed)
  truth = ones(Bool, n)
  for i in removed
    truth[i] = !truth[i]
  end
  return truth
end

function tally_venn(a, b)
  tp = sum(a .* b)
  fp = sum(!a .* b)
  fn = sum(a .* !b)
  tn = sum(!a .* !b)
  return [tp, fp, fn, tn]
end

function assess_filter(criteria, truth, steps)
  counts = []
  for i in steps
    push!(counts, tally_venn(truth, criteria .<= i))
  end
  return hcat(counts...)
end

"""
Display index k match pair meshes are two frames in movie to scroll between
"""
function review_matches_movie(meshset, k, step="post", downsample=2)
  matches = meshset.matches[k]
  src_index = matches.src_index
  dst_index = matches.dst_index
  println("display_matches_movie: ", (src_index, dst_index))
  src_mesh = meshset.meshes[find_index(meshset, src_index)]
  dst_mesh = meshset.meshes[find_index(meshset, dst_index)]

  if step == "post"
    src_nodes, dst_nodes = get_matched_points_t(meshset, k)
    src_index = (src_mesh.index[1], src_mesh.index[2], src_mesh.index[3]-1, src_mesh.index[4]-1)
    dst_index = (dst_mesh.index[1], dst_mesh.index[2], dst_mesh.index[3]-1, dst_mesh.index[4]-1)
    global_bb = get_global_bb(meshset)
    src_offset = [global_bb.i, global_bb.j]
    dst_offset = [global_bb.i, global_bb.j]
  else
    src_nodes, dst_nodes = get_matched_points(meshset, k)
    src_offset = src_mesh.disp
    dst_offset = dst_mesh.disp
  end

  src_nodes = hcat(src_nodes...)[1:2, :] .- src_offset
  dst_nodes = hcat(dst_nodes...)[1:2, :] .- dst_offset

  src_img = get_ufixed8_image(src_index)
  println(size(src_img))
  for i = 1:downsample
    src_img = restrict(src_img)
    src_nodes /= 2
  end
  dst_img = get_ufixed8_image(dst_index)
  println(size(dst_img))
  for i = 1:downsample
    dst_img = restrict(dst_img)
    dst_nodes /= 2
  end

  i_max = min(size(src_img,1), size(dst_img,1))
  j_max = min(size(src_img,2), size(dst_img,2))
  img = Image(cat(3, src_img[1:i_max, 1:j_max], dst_img[1:i_max, 1:j_max]), timedim=3)
  # imgc, img2 = view(make_isotropic(img))
  imgc, img2 = view(img)
  vectors = [src_nodes; dst_nodes]
  # an_pts, an_vectors = draw_vectors(imgc, img2, vectors)
  an_src_pts = draw_points(imgc, img2, src_nodes, RGB(1,0,0))
  an_dst_pts = draw_points(imgc, img2, dst_nodes, RGB(0,1,0))
  # pts_to_remove = edit_matches(imgc, img2, an_dst_pts)
  # return pts_to_remove
end

"""
Excerpt image in radius around given point
"""
function sliceimg(point, radius::Int64, img_size)
  point = ceil(Int64, point)
  imin = max(point[1]-radius, 1)
  jmin = max(point[2]-radius, 1)
  imax = min(point[1]+radius, img_size[1])
  jmax = min(point[2]+radius, img_size[2])
  return imin, imax, jmin, jmax
  # i_range = imin:imax
  # j_range = jmin:jmax
  # return img[i_range, j_range]
end

"""
Excerpt image in radius around given point
"""
function pt_radius_to_bounds(point, radius::Int64)
  point = ceil(Int64, point)
  imin = point[1]-radius
  jmin = point[2]-radius
  imax = point[1]+radius
  jmax = point[2]+radius
  return imin, imax, jmin, jmax
end

"""
Retrieve 1d array of block match pairs from the original images
"""
function get_blockmatch_images(meshset, k, mesh_type, radius)
  matches = meshset.matches[k]
  index = matches.(symbol(mesh_type, "_index"))
  mesh = meshset.meshes[find_index(meshset, index)]
  src_points, dst_points = get_matched_points(meshset, k)
  # src_points_t, dst_points_t = get_matched_points_t(meshset, k)
  scale = meshset.params["scaling_factor"]
  s = [scale 0 0; 0 scale 0; 0 0 1]
  img, _ = imwarp(get_ufixed8_image(mesh), s)
  offset = mesh.disp*scale
  # src_points *= scale
  # dst_points *= scale

  bm_imgs = []
  bm_bounds = []
  # for pt in (mesh_type == "src" ? src_points : dst_points)
  for pt in dst_points
    bounds = sliceimg(pt-offset, radius, size(img))
    push!(bm_imgs, img[range(bounds[1:2]...), range(bounds[3:4]...)])
    push!(bm_bounds, (bounds[1]+offset[1], bounds[3]+offset[2]))
    # push!(bm_imgs, sliceimg(img, pt-offset, radius))
  end
  return bm_imgs, bm_bounds
end

"""
Make one image from two blockmatch images, their difference, & their overlay
"""
function create_filtered_images(src_imgs, dst_imgs)
  gc(); gc();
  images = []
  for (n, (src_img, dst_img)) in enumerate(zip(src_imgs, dst_imgs))
    println(n, "/", length(src_imgs))
    src_img, dst_img = match_padding(src_img, dst_img)
    diff_img = convert(Image{RGB}, src_img-dst_img)
    imgA = grayim(Image(src_img))
    imgB = grayim(Image(dst_img))
    yellow_img = Overlay((imgA,imgB), (RGB(1,0,0), RGB(0,1,0)))
    pairA = convert(Image{RGB}, src_img)
    pairB = convert(Image{RGB}, dst_img)
    left_img = vcat(pairA, pairB)
    right_img = vcat(diff_img, yellow_img)
    img = hcat(left_img, right_img)
    push!(images, img)
  end
  return images
end

"""
Edit blockmatches with paging
"""
function edit_blockmatches(images)
  max_images_to_display = 60
  no_of_images = length(images)
  pages = ceil(Int, no_of_images / max_images_to_display)
  blockmatch_ids = Set()
  for page in 1:pages
    start = (page-1)*max_images_to_display + 1
    finish = start-1 + min(max_images_to_display, length(images)-start-1)
    println(start, finish)
    page_ids = display_blockmatches(images[start:finish], true, start)
    blockmatch_ids = union(blockmatch_ids, page_ids)
  end
  return blockmatch_ids
end

"""
Display blockmatch images in a grid to be clicked on
"""
function display_blockmatches(images, edit_mode_on=false, start_index=1)
  no_of_images = length(images)
  println("Displaying ", no_of_images, " images")
  grid_height = 6
  grid_width = 10
  # aspect_ratio = 1.6
  # grid_height = ceil(Int, sqrt(no_of_images/aspect_ratio))
  # grid_width = ceil(Int, aspect_ratio*grid_height)
  img_canvas_grid = canvasgrid(grid_height, grid_width, pad=1)
  match_index = 0
  blockmatch_ids = Set()
  e = Condition()

  function right_click(k)
    if in(k, blockmatch_ids)
      blockmatch_ids = setdiff(blockmatch_ids, Set(k))
    else
      push!(blockmatch_ids, k)
    end
    println(collect(blockmatch_ids))
  end

  for j = 1:grid_width
    for i = 1:grid_height
      match_index += 1
      n = match_index + start_index - 1
      if match_index <= no_of_images
        img = images[match_index]
        # imgc, img2 = view(img_canvas_grid[i,j], make_isotropic(img))
        imgc, img2 = view(img_canvas_grid[i,j], img)
        img_canvas = canvas(imgc)
        if edit_mode_on
          bind(img_canvas, "<Button-3>", path->right_click(n))
          win = Tk.toplevel(img_canvas)
          bind(win, "<Destroy>", path->notify(e))
        end
      end
    end
  end
  if edit_mode_on
    wait(e)
  end
  return blockmatch_ids
end

"""
Boolean if bounding boxes intersect
"""
function intersects(bbA::BoundingBox, bbB::BoundingBox)
  bb = bbA - bbB
  return !isequal(bb.i, NaN)
end

"""
Given list of bounding boxes, calculate pairs of indices with overlaps
"""
function find_overlaps(boundingboxes)
  bbs = copy(boundingboxes)
  overlap_tuples = []
  i = length(bbs)
  while length(bbs) != 0
    bbi = pop!(bbs)
    for (j, bbj)  in enumerate(bbs)
      if intersects(bbi, bbj)
        push!(overlap_tuples, (i,j))
      end
    end
    i -= 1
  end
  return overlap_tuples
end

"""
Crop image with offset to bounding box
"""
function imcrop(img, offset, bb)
  o = zeros(eltype(img), ceil(bb.h)+1, ceil(bb.w)+1)
  ibb = BoundingBox(offset..., size(img)...)
  d = bb - ibb
  o_start = abs(bb.i-d.i)+1:abs(bb.i-d.i) + d.h
  o_end = abs(bb.j-d.j)+1:abs(bb.j-d.j) + d.w
  im_start = abs(ibb.i-d.i)+1:abs(ibb.i-d.i) + d.h
  im_end = abs(ibb.j-d.j)+1:abs(ibb.j-d.j) + d.w
  o[o_start, o_end] = img[im_start, im_end]
  return o
end

function reshape_seam(img)
  size_threshold = 4000
  d = 10
  if size(img, 1) >= size_threshold && size(img, 2) <= size_threshold
    x = 5
    i = size(img, 2)
    m = ceil(Int64, size(img, 1)/x)
    o = zeros(eltype(img), m, x*(i+d))
    for k in 1:x
      o_slice = 1:m, (k-1)*(i+d)+1:k*i+(k-1)*d
      img_slice = (k-1)*m+1:k*m, 1:i
      if k == x
        img_slice = (k-1)*m+1:size(img,1), 1:i
        s = size(img[img_slice...], 1)
        o_slice = 1:s, (k-1)*(i+d)+1:k*i+(k-1)*d
      end
      o[o_slice...] = img[img_slice...]
    end
    img = o

  elseif size(img, 2) >= size_threshold && size(img, 1) <= size_threshold
    x = 3
    i = size(img, 1)
    m = ceil(Int64, size(img, 2)/x)
    o = zeros(eltype(img), x*(i+d), m)

    for k in 1:x
      o_slice = (k-1)*(i+d)+1:k*i+(k-1)*d, 1:m
      img_slice = 1:i, (k-1)*m+1:k*m
      if k == x
        img_slice = 1:i, (k-1)*m+1:size(img,2)
        s = size(img[img_slice...], 2)
        o_slice = (k-1)*(i+d)+1:k*i+(k-1)*d, 1:s
      end
      o[o_slice...] = img[img_slice...]
    end
    img = o
  end
  return img
end

function find_matches_index(meshset, indexA, indexB)
  matches_index = 0
  for (k, matches) in enumerate(meshset.matches)
    if matches.src_index == indexA && matches.dst_index == indexB
      return k
    end
  end
  return matches_index
end

function adjust_points_list_by_offset(points, offset)
  for k in 1:length(points)
    points[k] -= offset
  end
  return points
end


function transform_point(pt, tform)
  return ([pt..., 1]'*tform)[1:2]
end

function scale_matches(pts, scale)
  return [x*scale for x in pts]
end

function transform_matches(pts, tform)
  return [transform_point(x, tform) for x in pts]
end

function offset_matches(src_pts, dst_pts, offset)
  src_pts = [x - offset for x in src_pts]
  dst_pts = [x - offset for x in dst_pts]
  return src_pts, dst_pts
end

function get_matches(meshset, indexA, indexB)
  match_no = find_matches_index(meshset, indexA, indexB)
  if match_no > 0
    return get_matched_points_t(meshset, match_no)
  end
  return [], []
end

function make_thumbnail(img, vectors, colors, match_nums, factor=1.0)
  tf = 1.0
  vf = 1.0
  pf = 1.0
  surface = create_drawing(img')
  ctx = create_contex(surface)
  for (k, (vector, color, idx))  in enumerate(zip(vectors, colors, match_nums))
    draw_vectors(ctx, zip(vector...), color*vf, factor)
    draw_points(ctx, vector[1], color*pf)
    draw_indices(ctx, vector[1], [0,10], 22.0, color*tf)
    draw_text(ctx, "#$idx", [50+(k-1)*30,10], [0,0], 28.0, color)
  end
  draw_reference(ctx, size(img)..., factor)
  return get_drawing(surface)
end

function write_thumbnail(img, path, vectors, colors, match_nums, factor=1.0)
  drawing = make_thumbnail(img, vectors, colors, match_nums, factor)
  img = convert_drawing(drawing)
  println("Writing thumbnail\n\t", path)
  FileIO.save(path, reshape_seam(img))
end

"""
`WRITE_SEAMS` - Write out overlays of montaged seams
""" 
function write_seams(meshset, imgs, offsets, indices, fn_label="seam")
    bbs = []
    for (img, offset) in zip(imgs, offsets)
        push!(bbs, BoundingBox(offset..., size(img)...))
    end
    overlap_tuples = find_overlaps(bbs)
    for (k, (i,j)) in enumerate(overlap_tuples)
      println("Writing seam ", k, " / ", length(overlap_tuples))
      img, fuse_offset = imfuse(imgs[i], offsets[i], imgs[j], offsets[j])
      bb = bbs[i] - bbs[j]
      path = get_outline_filename(fn_label, indices[i], indices[j])
      img_cropped = imcrop(img, fuse_offset, bb)
      # imwrite(reshape_seam(img_cropped), path)
      offset = [bb.i, bb.j]
      vectorsA = get_matches(meshset, indices[i], indices[j])
      vectorsB = get_matches(meshset, indices[j], indices[i])
      vectorsA = offset_matches(vectorsA..., offset)
      vectorsB = offset_matches(vectorsB..., offset)
      matchij = find_matches_index(meshset, indices[i], indices[j])
      matchji = find_matches_index(meshset, indices[j], indices[i])
      vectors = (vectorsA, vectorsB)
      match_nums = (matchij, matchji)
      colors = ([0,0,0], [1,1,1])
      # factor = 100
      write_thumbnail(img_cropped, path, vectors, colors, match_nums)
    end
end

"""
Load aligned meshset mesh nodes
"""
function load_nodes(meshset, section_range=1:10, downsample=1)
  global_bb = GLOBAL_BB
  global_nodes = []
  for mesh in meshset.meshes[section_range]
    nodes = hcat(mesh.nodes_t...) .- collect((global_bb.i, global_bb.j))
    for i = 1:downsample
      nodes /= 2
    end
    push!(global_nodes, [nodes])
  end
  return global_nodes
end

"""
Provide prompt to user to enter int array elements one at a time
"""
function enter_array()
  println("Type integers and press return. Press return twice to end.")
  output = []
  is_int = true
  while is_int
    a = readline(STDIN)
    if typeof(parse(a)) == Int64
      push!(output, parse(a))
    else
      is_int = false
    end
  end
  return output
end

"""
Display outline with matches plotted
"""
function plot_matches_outline(meshset, match_no, factor=5)
  scale = 0.20
  matches = meshset.matches[match_no]
  src_index = matches.src_index
  dst_index = matches.dst_index
  src_mesh = meshset.meshes[find_index(meshset, src_index)]
  dst_mesh = meshset.meshes[find_index(meshset, dst_index)]

  src_bb = find_mesh_bb(points_to_Nx3_matrix(src_mesh.nodes_t*scale))
  dst_bb = find_mesh_bb(points_to_Nx3_matrix(dst_mesh.nodes_t*scale))
  bb = src_bb+dst_bb
  img, offset = create_image(bb)
  offset = collect((offset..., 0))

  src_pts, dst_pts = get_matched_points_t(meshset, match_no)
  # src_pts = src_mesh.nodes
  # dst_pts = src_mesh.nodes_t
  src_pts = points_to_Nx3_matrix(src_pts*scale) .- offset'
  dst_pts = points_to_Nx3_matrix(dst_pts*scale) .- offset'

  set_canvas_size(img[1], bb.w/2, bb.h/2)
  src_bb = BoundingBox(src_bb.i-bb.i, src_bb.j-bb.j, src_bb.h, src_bb.w)
  dst_bb = BoundingBox(dst_bb.i-bb.i, dst_bb.j-bb.j, dst_bb.h, dst_bb.w)
  draw_box(img..., src_bb, RGB(0,1,0))
  draw_box(img..., dst_bb, RGB(1,0,0))

  vectors = [src_pts[:,1:2]'; dst_pts[:,1:2]']
  draw_vectors(img..., vectors, RGB(0,0,1), RGB(1,0,1), factor)
  draw_indices(img..., vectors[1:2,:], 14.0, [-8,-20])
  draw_reference_vector(img..., factor*scale, bb)
  return img
end

function indices2string(indexA, indexB)
  return string(join(indexA[1:2], ","), "-", join(indexB[1:2], ","))
end

function get_outline_filename(prefix, src_index, dst_index=(0,0,0,0))
  dir = ALIGNED_DIR
  indstring = indices2string(src_index, dst_index)
  if dst_index[1] == 0
    indstring = join(src_index[1:2], ",")
  end
  fn = string(prefix, "_", indstring, ".h5")
  if is_premontaged(src_index)
    dir = MONTAGED_DIR
    ind = string(join(src_index, ","), "-", join(dst_index, ","))
    fn = string(prefix, "_", ind, ".jpg")
  elseif is_montaged(src_index)
    dir = PREALIGNED_DIR
    fn = string(prefix, "_", indstring, ".tif")
  elseif is_prealigned(src_index)
    dir = PREALIGNED_DIR
    fn = string(prefix, "_", indstring, ".tif")
  end
  # return joinpath(dir, "review/151106_1,2-1,16_cleaned", fn)
  return joinpath(dir, "review", fn)
end

function write_meshset_match_outlines(meshset, factor=5)
  for k in 1:length(meshset.matches)
    imgc, img2 = plot_matches_outline(meshset, k, factor)
    matches = meshset.matches[k]
    path = get_outline_filename("outline", matches.src_index, matches.dst_index)
    path = string(path[1:end-4], ".png")
    write_canvas(imgc, path)
    close_image(imgc)
  end
end

function write_dir_outlines(dir, start=1, finish=0, factor=5)
  if finish == 0
    finish = length(sort_dir(dir, "jld"))
  end
  for fn in sort_dir(dir, "jld")[start:finish]
    meshset = JLD.load(joinpath(dir, fn))["MeshSet"]
    write_meshset_match_outlines(meshset, factor)
  end
end

"""
Display outline for montaged meshsets with matches plotted
"""
function plot_matches_outline(meshset, factor=5)
  scale = 0.05
  bb_list = [find_mesh_bb(points_to_Nx3_matrix(m.nodes*scale)) for m in meshset.meshes]
  bb = sum(bb_list)
  offset = [bb.i, bb.j, 0]
  src_pts = []
  dst_pts = []
  for (k, matches) in enumerate(meshset.matches)
    src, dst = get_matched_points_t(meshset, k)
    src = points_to_Nx3_matrix(src*scale) .- offset'
    dst = points_to_Nx3_matrix(dst*scale) .- offset'
    push!(src_pts, src)
    push!(dst_pts, dst)
  end
  src_pts = vcat(src_pts...)
  dst_pts = vcat(dst_pts...)

  img, offset = create_image(bb)
  set_canvas_size(img[1], bb.w/2, bb.h/2)
  # draw_box(img..., src_bb, RGB(0,1,0))
  # draw_box(img..., dst_bb, RGB(1,0,0))

  vectors = [src_pts[:,1:2]'; dst_pts[:,1:2]']
  draw_vectors(img..., vectors, RGB(0,0,1), RGB(1,0,1), factor)
  # draw_indices(img..., vectors[1:2,:], 14.0, [-8,-20])
  return img[1], img[2]
end

function review_montages(username, i=1)
  path = joinpath(MONTAGED_DIR, "review")
  assert(isdir(path))
  review_path = joinpath(path, string("montage_review_", username, ".txt"))
  break_review = false
  d = readdir(path)
  fn = d[i]
  println(i, " / ", length(d) , ": ", fn)
  img = Images.imread(joinpath(path, fn))
  imgc, img2 = view(img, pixelspacing=[1,1])

  e = Condition()
  c = canvas(imgc)
  win = Tk.toplevel(c)
  bind(win, "g", path->good_exit())
  bind(win, "b", path->bad_exit())
  bind(win, "<KP_0>", path->good_exit())
  bind(win, "<KP_4>", path->bad_exit())
  bind(win, "<KP_7>", path->go_back())
  bind(win, "<KP_9>", path->go_next())
  # bind(win, "<space>", path->toggle_rating())
  bind(win, "<Escape>", path->exit_loop())
  # bind(win, "<Enter>", path->end_review())
  bind(win, "<Destroy>", path->exit_loop())

  isgood = true
  j = i

  function good_exit()
    println("good")
    log_montage_review(username, fn, true)
    go_next()
  end

  function bad_exit()
    println("bad")
    log_montage_review(username, fn, false)
    go_next()
  end

  function go_next()
    j = i+1
    end_review()
  end

  function go_back()
    j = i-1
    end_review()
  end

  # function toggle_rating()
  #   isgood = !isgood
  #   isgood ? println("good") : println("bad")
  #   log_montage_review(username, fn, isgood)
  # end

  function exit_loop()
    break_review = true
    # log_montage_review(username, fn, isgood)
    end_review()
  end

  function end_review()
    bind(win, "g", path->path)
    bind(win, "b", path->path)
    bind(win, "<KP_0>", path->path)
    bind(win, "<KP_4>", path->path)
    bind(win, "<KP_7>", path->path)
    bind(win, "<KP_9>", path->path)
    # bind(win, "<space>", path->path)
    bind(win, "<Escape>", path->path)
    # bind(win, "<Enter>", path->path)
    bind(win, "<Destroy>", path->path)
    destroy(win)
    notify(e)
  end

  wait(e)

  if !break_review
    return review_montages(username, j)
  end
  return i
end

function get_montage_review_path(username)
  return joinpath(INSPECTION_DIR, string("montage_review_", username, ".txt"))
end

"""
Write points to remove in a text file
"""
function log_montage_review(username, fn, isgood)
  path = get_montage_review_path(username)
  ts = Dates.format(now(), "yymmddHHMMSS")
  row = [ts, username, fn, Int(isgood)]'
  if !isfile(path)
    f = open(path, "w")
    close(f)
    table = row
  else  
    table = readdlm(path)
    i = findfirst(table[:,3], fn)
    if i != 0
      table = vcat(table[1:i-1, :], row, table[i+1:end,:])
    else
      table = vcat(table, row)
    end
  end
  table = table[sortperm(table[:, 3]), :]
  # println("Saving montage_review:\n", path)
  writedlm(path, table)
end

"""
Combine stack error log files of all tracers to create one log matrix
"""
function compile_tracer_montage_review_logs()
  tracers = ["hmcgowan", "bsilverman", "merlinm", "kpw3"]
  logs = []
  for tracer in tracers
    path = get_montage_review_path(tracer)
    if isfile(path)
      push!(logs, readdlm(path))
    end
  end
  return vcat(logs...)
end

function show_montage_review_progress()
  path = joinpath(MONTAGED_DIR, "review")
  montages = readdir(path)
  x = 1:length(montages)
  y = zeros(Int64, length(montages))
  logs = compile_tracer_montage_review_logs()
  tracers = unique(logs[:, 2])
  for tracer in tracers
    log = logs[logs[:,2] .== tracer, :]
    log_time = round(Int64, (log[:,1] % 10^6) / 100)
    log_k = map(x -> findfirst(montages, x), log[:,3])
    y[log_k] = log_time
    PyPlot.plot(log_k, log_time, ".")
  end
  PyPlot.plot(x[y.==0], y[y.==0], ".")
  println("Reviewed stack columns: ", length(unique(logs[:,3])), " / ", length(montages))
end
