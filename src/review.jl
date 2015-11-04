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
  println("Removing ", length(indices_to_remove), " points")
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

"""
Provide bindings for GUI to right-click and remove matches from a mesh. End 
manual removal by exiting the image window, or by pressing enter while focused
on the image window.

Args:

* imgc: ImageCanvas object (from image with annotations)
* img2: ImageZoom object (from image with annotations)
* annotation: the annotation object from the image

Returns:

* pts_to_remove: array of indices for points to remove from mesh

  pts_to_remove = edit_matches(imgc, img2, annotation)
"""
function edit_matches(imgc, img2, annotation)
    e = Condition()

    pts_to_remove = Array{Integer,1}()
    pts = copy(annotation.ann.data.pts)

    c = canvas(imgc)
    win = Tk.toplevel(c)
    bind(c, "<Button-3>", (c, x, y)->right_click(parse(Int, x), parse(Int, y)))
    bind(win, "<Destroy>", path->end_edit())

    function right_click(x, y)
        xu,yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
        xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
        # println(xi, ", ", yi)
        limit = (img2.zoombb.xmax - img2.zoombb.xmin) * 0.0125 # 100/8000
        annpts = annotation.ann.data.pts
        annidx = find_idx_of_nearest_pt(annpts, [xi, yi], limit)
        if annidx > 0
            idx = find_idx_of_nearest_pt(pts, [xi, yi], limit)
            println(idx, ": ", [xi, yi])
            annotation.ann.data.pts = hcat(annpts[:,1:annidx-1], 
                                                        annpts[:,annidx+1:end])
            ImageView.redraw(imgc)
            push!(pts_to_remove, idx)
        end
    end

    function end_edit()
        println("End edit\n")
        notify(e)
        bind(c, "<Button-3>", path->path)
        bind(win, "<Destroy>", path->path)
    end

    println("Right click to remove correspondences, then exit window.")
    wait(e)

    return pts_to_remove
end

"""
Display index k match pair meshes are two frames in movie to scroll between
"""
function review_matches(meshset, k)
  matches = meshset.matches[k]
  indexA = (matches.src_index[1:2]..., -4, -4)
  indexB = (matches.dst_index[1:2]..., -4, -4)

  path = get_outline_filename("thumb_imfuse", indexB, indexA)
  img = h5read(path, "img")
  offset = h5read(path, "offset")
  scale = h5read(path, "scale")

  src_nodes, dst_nodes = get_matched_points_t(meshset, k)
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vectors = offset_matches(vectorsA, vectorsB, offset)
  vecs = [hcat(vectors[1]...); hcat(vectors[2]...)];

  imview = view(img, pixelspacing=[1,1])
  an_pts, an_vectors = show_vectors(imview..., vecs, RGB(0,0,1), RGB(1,0,1), 10)
  pts_to_remove = edit_matches(imview..., an_pts)
  return pts_to_remove
end

"""
Write points to remove in a text file
"""
function store_points(path, k, pts_to_remove)
  if !isfile(path)
    f = open(path, "w")
    close(f)
    pts = [k, join(pts_to_remove, ",")]'
  else  
    pts = readdlm(path)
    idx = findfirst(pts[:,1], k)
    if idx != 0
      pts[idx, 2] = join(pts_to_remove, ",")
    else
      pts_line = [k, join(pts_to_remove, ",")]'
      pts = vcat(pts, pts_line)
    end
  end
  pts = pts[sortperm(pts[:, 1],), :]
  println("Saving pts_to_remove:\n\t", path)
  writedlm(path, pts)
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
Append 'EDITED' with a timestamp to update filename
"""
function update_filename(fn, ts="")
  if ts == ""
    ts = Dates.format(now(), "yyyymmddHHMMSS")
  end
  return string(fn[1:end-4], "_EDITED_", ts, fn[end-3:end])
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

# """
# Determine appropriate meshset solution method for meshset
# """
# function resolve!(meshset::MeshSet)
#   index = meshset.meshes[1].index
#   if is_montaged(index)
#   elseif is_prealigned(index)
#     solve_meshset!(meshset)
#   end
# end
#
# """
# Review all matches in the meshsets in given directory via method specified
# """
# function review_matches(dir, method="movie")
#   filenames = sort_dir(dir)
#   for fn in filenames[1:1]
#     println(joinpath(dir, fn))
#     is_file_changed = false
#     meshset = load(joinpath(dir, fn))["MeshSet"]
#     # src_index = meshset.meshes[1].index
#     # dst_index = meshset.meshes[2].index
#     # meshset.meshes[1].index = (src_index[1:2]..., PREALIGNED_INDEX, PREALIGNED_INDEX)
#     # meshset.meshes[2].index = (dst_index[1:2]..., PREALIGNED_INDEX, PREALIGNED_INDEX)
#     new_path = update_filename(fn)
#     log_file = open(joinpath(dir, string(new_path[1:end-4], ".txt")), "w")
#     write(log_file, "Meshset from ", joinpath(dir, fn), "\n")
#     for (k, matches) in enumerate(meshset.matches)
#       println("Inspecting matches at index ", k,)
#       if method=="images"
#         @time src_bm_imgs, _ = get_blockmatch_images(meshset, k, "src", meshset.params["block_size"]-400)
#         @time dst_bm_imgs, _ = get_blockmatch_images(meshset, k, "dst", meshset.params["block_size"]-400)
#         @time filtered_imgs = create_filtered_images(src_bm_imgs, dst_bm_imgs)
#         pts_to_remove = collect(edit_blockmatches(filtered_imgs))
#       else
#         pts_to_remove = review_matches_movie(meshset, k, "pre")
#       end
#       if length(pts_to_remove) > 0
#         is_file_changed = true
#         # save_blockmatch_imgs(meshset, k, pts_to_remove, joinpath(dir, "blockmatches"))
#       end
#       remove_matches_from_meshset!(pts_to_remove, k, meshset)
#       log_line = string("Removed following matches from index ", k, ":\n")
#       log_line = string(log_line, join(pts_to_remove, "\n"))
#       write(log_file, log_line, "\n")
#     end
#     if is_file_changed
#       resolve!(meshset)
#       println("Saving JLD here: ", new_path)
#       @time save(joinpath(dir, new_path), meshset)
#     end
#     close(log_file)
#   end
# end

# """
# Convert matches object for section into filename string
# """
# function matches2filename(meshset, k)
#   src_index = meshset.matches[k].src_index
#   dst_index = meshset.matches[k].dst_index
#   return string(join(dst_index[1:2], ","), "-", join(src_index[1:2], ","))
# end

# """
# Convert cross correlogram to Image for saving
# """
# function xcorr2Image(xc)
#   return grayim((xc .+ 1)./2)
# end

# """
# Save match pair images w/ search & block radii at k index in meshset
# """
# function save_blockmatch_imgs(meshset, k, blockmatch_ids=[], path=joinpath(ALIGNED_DIR, "blockmatches"))
#   dir_path = joinpath(path, matches2filename(meshset, k))
#   if !isdir(dir_path)
#     mkdir(dir_path)
#   end
#   block_radius = meshset.params["block_size"]
#   search_radius = meshset.params["search_r"]
#   scale = meshset.params["scaling_factor"]
#   combined_radius = search_radius+block_radius
#   src_imgs, src_bounds = get_blockmatch_images(meshset, k, "src", combined_radius)
#   dst_imgs, dst_bounds = get_blockmatch_images(meshset, k, "dst", block_radius)
#   src_points, dst_points = get_matched_points(meshset, k)
#   # dst_imgs_adjusted = get_blockmatch_images(meshset, k, "dst", block_radius)
#   println("save_blockmatch_imgs")
#   if blockmatch_ids == []
#     blockmatch_ids = 1:length(src_imgs)
#   end
#   for (idx, (src_img, dst_img, src_point, src_bound)) in enumerate(zip(src_imgs, 
#                                                                 dst_imgs, 
#                                                                 src_points,
#                                                                 src_bounds))
#     if idx in blockmatch_ids
#       println(idx, "/", length(src_imgs))
#       # xc = normxcorr2(src_img, dst_img)
#       xc_peak, xc = get_max_xc_vector(dst_img, src_img)
#       println(size(src_img))
#       println(size(dst_img))
#       println(size(xc))
#       n = @sprintf("%03d", idx)
#       img_mark = "good"
#       dst_path = joinpath(dir_path, string(n , "_1_dst_", k, "_", img_mark, ".png"))
#       imwrite(dst_img, dst_path)
#       src_offset = src_point - collect(src_bound)
#       src_path = joinpath(dir_path, string(n , "_2_src_", k, "_", img_mark, ".png"))
#       imwrite_box(src_img, src_offset, block_radius, src_path)
#       # imwrite(dst_img, joinpath(dir_path, string(n , "_dst_", k, "_", img_mark, ".jpg")))
#       if !isnan(sum(xc))
#         r_max = maximum(xc);
#         rad = round(Int64, (size(xc, 1) - 1)/ 2)  
#         ind = findfirst(r_max .== xc)
#         xc_peak = [rem(ind, size(xc, 1)), cld(ind, size(xc, 1))]
#         xc_path = joinpath(dir_path, string(n , "_3_xc_", k, "_", img_mark, ".png"))
#         imwrite_box(xcorr2Image(xc'), xc_peak, 20, xc_path)
#         # imwrite(xcorr2Image(xc'), joinpath(dir_path, string(n , "_xc_", k, "_", img_mark, ".jpg")))
#       end
#     end
#   end
# end

"""
Excerpt image in radius around given point
"""
function sliceimg(img, point, radius::Int64)
  point = ceil(Int64, point)
  imin = max(point[1]-radius, 1)
  jmin = max(point[2]-radius, 1)
  imax = min(point[1]+radius, size(img,1))
  jmax = min(point[2]+radius, size(img,2))
  return imin, jmin, imax, jmax
  # i_range = imin:imax
  # j_range = jmin:jmax
  # return img[i_range, j_range]
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
    imin, jmin, imax, jmax = sliceimg(img, pt-offset, radius)
    push!(bm_imgs, img[imin:imax, jmin:jmax])
    push!(bm_bounds, (imin+offset[1], jmin+offset[2]))
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
    src_pts, dst_pts = get_matched_points_t(meshset, match_no)
  end
  return src_pts, dst_pts
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
  Images.save(path, reshape_seam(img))
end

"""
`WRITE_SEAMS` - Write out overlays of montaged seams
""" 
function write_seams(meshset, imgs, offsets, indices)
    bbs = []
    for (img, offset) in zip(imgs, offsets)
        push!(bbs, BoundingBox(offset..., size(img)...))
    end
    overlap_tuples = find_overlaps(bbs)
    for (k, (i,j)) in enumerate(overlap_tuples)
      println("Writing seam ", k, " / ", length(overlap_tuples))
      img, fuse_offset = imfuse(imgs[i], offsets[i], imgs[j], offsets[j])
      bb = bbs[i] - bbs[j]
      path = get_outline_filename("seam", indices[i], indices[j])
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
      factor = 100
      write_thumbnail(img_cropped, path, vectors, colors, match_nums, factor)
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
Ripped from ImageView > navigation.jl
"""
function stop_playing!(state::ImageView.NavigationState)
    if state.timer != nothing
        close(state.timer)
        state.timer = nothing
    end
end

"""
Ripped from ImageView > navigation.jl
"""
function updatet(ctrls, state)
  Tk.set_value(ctrls.editt, string(state.t))
  Tk.set_value(ctrls.scalet, state.t)
  enableback = state.t > 1
  Tk.set_enabled(ctrls.stepback, enableback)
  Tk.set_enabled(ctrls.playback, enableback)
  enablefwd = state.t < state.tmax
  Tk.set_enabled(ctrls.stepfwd, enablefwd)
  Tk.set_enabled(ctrls.playfwd, enablefwd)
end

"""
Ripped from ImageView > navigation.jl
"""
function incrementt(inc, ctrls, state, showframe)
    state.t += inc
    updatet(ctrls, state)
    showframe(state)
end

"""
Ripped from ImageView > navigation.jl
"""
function playt(inc, ctrls, state, showframe)
    if !(state.fps > 0)
        error("Frame rate is not positive")
    end
    stop_playing!(state)
    dt = 1/state.fps
    state.timer = Timer(timer -> stept(inc, ctrls, state, showframe), dt, dt)
end

"""
Ripped from ImageView > navigation.jl
"""
function stept(inc, ctrls, state, showframe)
    if 1 <= state.t+inc <= state.tmax
        incrementt(inc, ctrls, state, showframe)
    else
        stop_playing!(state)
    end
end

"""
Create loop of the image
"""
function start_loop(imgc, img2, fps=6)
  state = imgc.navigationstate
  ctrls = imgc.navigationctrls
  showframe = state -> ImageView.reslice(imgc, img2, state)
  set_fps!(state, fps)

  if !(state.fps > 0)
      error("Frame rate is not positive")
  end
  stop_playing!(state)
  dt = 1/state.fps
  state.timer = Timer(timer -> loopt(ctrls, state, showframe), dt, dt)
end

"""
Higher level call for ImageView outputs
"""
function stop_loop(imgc)
  stop_playing!(imgc.navigationstate)
end

"""
Endlessly repeat forward loop (continuous stept)
"""
function loopt(ctrls, state, showframe)
  inc = 1
  if state.t+inc < 1 || state.tmax < state.t+inc
      state.t = 0
  end
  incrementt(inc, ctrls, state, showframe)
end

"""
Building on ImageView > navigation.jl
"""
function set_fps!(state, fps)
  state.fps = fps
end

"""
Cycle through sections of the stack movie, with images staged for easier viewing
"""
function view_stack(meshset, section_range=(1:2), slice=(1:200,1:200), perm=[1,2,3])
  imgs = []
  for mesh in meshset.meshes[section_range]
    index = (mesh.index[1:2]..., mesh.index[3]-1, mesh.index[4]-1)
    img = get_h5_slice(get_h5_path(index), slice)
    push!(imgs, img)
  end

  println(slice)
  # img_stack = cat(3, imgs..., reverse(imgs[2:end-1])...)  # loop it
  img_stack = cat(3, imgs...)
  # img_stack = permutedims(cat(3, imgs...), perm)
  # img_movie = Image(img_stack, timedim=3)
  # if perm != [1,2,3]
    imgc, img2 = view(Image(permutedims(img_stack, [3,2,1]), timedim=3))
    start_loop(imgc, img2, 10)
    imgc, img2 = view(Image(permutedims(img_stack, [1,3,2]), timedim=3))
    start_loop(imgc, img2, 10)
    imgc, img2 = view(Image(permutedims(img_stack, [1,2,3]), timedim=3), pixelspacing=[1,1])
    start_loop(imgc, img2, 10)
  # end
  # start_loop(imgc, img2, 10)

  e = Condition()
  c = canvas(imgc)
  win = Tk.toplevel(c)

  function exit_movie()
    stop_loop(imgc)
    notify(e)
  end
  bind(win, "<Destroy>", path->notify(e))
  bind(win, "<Escape>", path->exit_movie())

  wait(e)
  destroy(win)
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