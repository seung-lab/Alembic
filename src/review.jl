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
    bind(win, "<Return>", path->end_edit())
    bind(win, "<KP_Enter>", path->end_edit())
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
        println("End edit")
        notify(e)
        bind(c, "<Button-3>", path->path)
        bind(win, "<Return>", path->path)
        bind(win, "<KP_Enter>", path->path)
        bind(win, "<Destroy>", path->path)
    end

    println("Right click to remove correspondences, then press enter.")
    wait(e)

    return pts_to_remove
end

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
function remove_matches_from_meshset!(indices_to_remove, match_index, meshset)
  matches = meshset.matches[match_index]
  flag = trues(matches.n)
  flag[indices_to_remove] = false
  matches.src_pointIndices = matches.src_pointIndices[flag]
  matches.dst_points = matches.dst_points[flag]
  matches.dst_triangles = matches.dst_triangles[flag]
  matches.dst_weights = matches.dst_weights[flag]
  matches.dispVectors = matches.dispVectors[flag]
  no_pts_removed = length(indices_to_remove)
  matches.n -= no_pts_removed
  meshset.m -= no_pts_removed
  meshset.m_e -= no_pts_removed
end

"""
Count number of displacement vector lengths k-sigma from the group mean

Args:

* vectors: 4xN array of points (start and end points of displacement vectors)
* k: number of sigmas to set the threshold for counting

Returns:

* integer for number of displacement vectors above the k-sigma threshold

  num = count_outliers(vectors, k)
"""
function count_outliers(vectors, k)
  d = sum((vectors[1:2,:] - vectors[3:4,:]).^2, 1).^(1/2)
  d_mean = mean(d)
  d_std = mean((d.-d_mean).^2).^(1/2)
  return sum((d.-d_mean)./d_std .> k)
end

"""
Loop through meshset & manually intervene at bad seams (3-sigma displacements)

Args:

* meshset: MeshSet with Matches objects

Returns:
* (updated meshset)

  fix_bad_seams(meshset)
"""
function review_matches_side_by_side(meshset, k, downsample=4)
  matches = meshset.matches[k]
  src_index = matches.src_index
  dst_index = matches.dst_index
  println("fix_bad_matches: ", (src_index, dst_index))
  src_mesh = meshset.meshes[findIndex(meshset, src_index)]
  dst_mesh = meshset.meshes[findIndex(meshset, dst_index)]
  src_nodes, dst_nodes = get_matched_points(meshset, k)
  src_offset = src_mesh.disp
  dst_offset = dst_mesh.disp
  src_nodes = hcat(src_nodes...)
  dst_nodes = hcat(dst_nodes...)  

  src_img = getFloatImage(src_mesh)
  for i = 1:downsample/2
    src_img = restrict(src_img)
    src_offset /= 2
    src_nodes /= 2
  end
  dst_img = getFloatImage(dst_mesh)
  for i = 1:downsample/2
    dst_img = restrict(dst_img)
    dst_offset /= 2
    dst_nodes /= 2
  end
  vectors = vcat(dst_nodes .- dst_offset, src_nodes .- dst_offset)

  c = canvasgrid(1,2)
  imgc, img2 = view(c[1,1], make_isotropic(src_img))
  draw_points(imgc, img2, src_nodes .- src_offset, RGB(0,0,1))
  imgc, img2 = view(c[1,2], make_isotropic(dst_img))
  imgc, img2, points, lines = draw_vectors(imgc, img2, vectors, RGB(1,0,0))
  pts_to_remove = edit_matches(imgc, img2, points)
  println(pts_to_remove)
  if length(pts_to_remove) > 0
    remove_matches_from_meshset!(pts_to_remove, k, meshset)
  end
  c = 0; src_img = 0; dst_img = 0; imgc = 0; img2 = 0;
  gc()
end

"""
Run meshset back through elastic solver and save
"""
function resolve_meshset(meshset)
  match_coeff = 10
  eta_gradient = 0.01
  eta_newton = 1.0
  show_plot = false
  grad_threshold = 1/1000
  newton_tolerance = 10.0^-7
  @time solveMeshSet!(meshset, match_coeff, eta_gradient, 
                        eta_newton, grad_threshold, newton_tolerance)
end

"""
Review all matches in a meshset
"""
function review_matches_in_matchset(meshset)
  for (idx, mesh) in enumerate(meshset.meshes)
    review_matches_side_by_side(meshset, idx)
  end
  resolve_meshset(meshset)
end

"""
Run through directory and review all meshes
"""
function fix_meshsets_in_dir(dir)
  jld_filenames = filter(x -> x[end-2:end] == "jld", readdir(dir))
  for fn in jld_filenames[1:1]
    println("Reviewing ", fn[1:end-4])
    meshset = load(joinpath(dir, fn))["MeshSet"]
    # review_matches_in_matchset(meshset)
    # @time MeshModule.MeshSet2JLD(joinpath(dir, fn[1:end-4], "_EDITED", now(), ".jld"), meshset)
  end
end

"""
Return first JLD file in provided directory
"""
function load_sample_meshset(dir)
  jld_filenames = filter(x -> x[end-2:end] == "jld", readdir(dir))
  fn = jld_filenames[1]
  return load(joinpath(dir, fn))["MeshSet"]
end

function sliceimg(img, point, radius)
  point = ceil(Int64, point)
  i_range = point[1]-radius:point[1]+radius
  j_range = point[2]-radius:point[2]+radius
  return sub(img, i_range, j_range)
end

"""
Retrieve 1d array of block match pairs from the original images
"""
function get_blockmatch_images(meshset, k, radius=block_size)
  matches = meshset.matches[k]
  src_index = matches.src_index
  dst_index = matches.dst_index
  src_mesh = meshset.meshes[findIndex(meshset, src_index)]
  dst_mesh = meshset.meshes[findIndex(meshset, dst_index)]
  src_points, dst_points = get_matched_points(meshset, k)
  # src_points_t, dst_points_t = get_matched_points_t(meshset, k)
  src_offset = src_mesh.disp
  dst_offset = dst_mesh.disp
  src_nodes = hcat(src_points...)
  dst_nodes = hcat(dst_points...)
  src_img = getFloatImage(src_mesh)
  dst_img = getFloatImage(dst_mesh)

  src_bm_imgs = []
  dst_bm_imgs = []

  # block_size in Params
  # search_r in Params
  for (src_pt, dst_pt) in zip(src_points, dst_points)
    push!(src_bm_imgs, sliceimg(src_img, src_pt-src_offset, radius))
    push!(dst_bm_imgs, sliceimg(dst_img, dst_pt-dst_offset, radius))
  end
  return src_bm_imgs, dst_bm_imgs
end

"""
Allow selecting of blockmatch images to remove correspondence
"""
function select_blockmatch(imgc, img2, k)
    c = canvas(imgc)
    bind(c, "<Button-3>", path->return k)
end

function right_click(k)
  println(k)
  # push!(blockmatch_ids, k)
end

"""
Display blockmatch images in a grid to be clicked on
"""
function display_blockmatches(src_imgs, dst_imgs, style="pair", aspect_ratio=1.4)
  a = ceil(Int, sqrt(length(src_imgs)/aspect_ratio))
  b = ceil(Int, aspect_ratio*a)
  c = canvasgrid(a, b, pad=5)
  n = 0
  blockmatch_ids = []
  e = Condition()

  for j = 1:b
    for i = 1:a
      n += 1
      if n <= length(src_imgs)
        if style=="diff"
          img = src_imgs[n] - dst_img[n]
        elseif style=="yellow"
          A = grayim(Image(src_imgs[n]))
          B = grayim(Image(dst_imgs[n]))
          img = Overlay((A,B), (RGB(1,0,0), RGB(0,1,0)))
        else
          img = vcat(src_imgs[n], ones(2, size(src_imgs[n], 2)) .* 255, dst_imgs[n])
        end
        imgc, img2 = view(c[i,j], make_isotropic(img))
        img_canvas = canvas(imgc)
        bind(img_canvas, "<Button-3>", img_canvas->right_click(n))
        win = Tk.toplevel(img_canvas)
        bind(win, "<Destroy>", path->notify(e))
      end
    end
  end
  wait(e)
  return blockmatch_ids 
end