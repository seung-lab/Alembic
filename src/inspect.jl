function show_blockmatch(match, ind, params)
  src_patch, src_pt, dst_patch, dst_pt, xc, offset = get_correspondence_patches(match, ind)
  block_r = params["block_r"]
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
Remove points displayed on ImageView with right-click
"""
function edit_matches(imgc, img2, matches, vectors, params)
    e = Condition()

    break_review = false
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

    c = canvas(imgc)
    win = Tk.toplevel(c)
    bind(c, "<Button-3>", (c, x, y)->inspect_match(parse(Int, x), parse(Int, y)))
    bind(c, "<Control-Button-3>", (c, x, y)->remove_match(parse(Int, x), parse(Int, y)))
    bind(win, "<Delete>", path->path)
    bind(win, "<Escape>", path->exit_loop())
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
          # annpt = ([lines[2,annidx], lines[1,annidx]] + offset) / scale 
          annpt = [lines[2,annidx], lines[1,annidx]]
          println(annidx, ": ", annpt)
          idx = find_idx_of_nearest_pt(vectors[1:2,:], annpt, 10)
          if idx > 0
            ptA = vectors[1:2,idx] # - params["src_offset"]
            ptB = vectors[3:4,idx] # - params["dst_offset"]
            println(idx, ": ", ptA, ", ", ptB)
            # keep_point = display_blockmatch(params, ptA, ptB)
            # if !keep_point
            #   mask[idx] = false
            #   update_annotations(imgc, img2, original_lines, mask)
            #   push!(indices_to_remove, idx)
            # end
            inspect_window = show_blockmatch(matches, idx, params)
            # inspect_window = display_blockmatch(params, ptA, ptB)
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

  function exit_loop()
    break_review = true
    end_edit()
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

    # println("1) Right click to inspect correspondences\n",
    #         "2) Ctrl + right click to remove correspondences\n",
    #         "3) Exit window or press Escape")
    wait(e)

    return indices_to_remove, break_review
end

"""
Display matches index k as overlayed images with points
"""
function inspect_matches(meshset, k)
  matches = meshset.matches[k]
  indexA = matches.src_index
  indexB = matches.dst_index

  path = get_outline_filename("seam", indexB, indexA)
  if !isfile(path)
    path = get_outline_filename("seam", indexA, indexB)
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
  params["scale"] = scale
  params["thumb_offset"] = offset
  params["src_index"] = matches.src_index
  params["dst_index"] = matches.dst_index
  params["src_offset"] = get_offset(meshset.meshes[find_mesh_index(meshset, indexA)])
  params["dst_offset"] = get_offset(meshset.meshes[find_mesh_index(meshset, indexB)])
  params["src_size"] = get_image_sizes(matches.src_index)
  params["dst_size"] = get_image_sizes(matches.dst_index)

  # src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences(meshset, k)
  # vectors = [hcat(src_nodes...); hcat(dst_nodes...)]
  src_nodes, dst_nodes, filtered_inds = get_globalized_correspondences_post(meshset, k)
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vecs = offset_matches(vectorsA, vectorsB, offset)
  vectors = [hcat(vecs[1]...); hcat(vecs[2]...)]

  imview = view(img, pixelspacing=[1,1])
  big_vecs = change_vector_lengths([hcat(vecs[1]...); hcat(vecs[2]...)], 20)
  an_pts, an_vectors = show_vectors(imview..., big_vecs, RGB(0,0,1), RGB(1,0,1))
  return imview, matches, vectors, copy(an_vectors.ann.data.lines), params
end

function get_montage_review_path(username)
  return joinpath(INSPECTION_DIR, string("montage_inspection_", username, ".txt"))
end

function inspect_montages(username, meshset_ind, match_ind)
  path = get_montage_review_path(username)
  indrange = get_index_range((1,1,-2,-2), (8,173,-2,-2))
  meshset = load(indrange[meshset_ind])
  println(meshset_ind, ": ", indrange[meshset_ind], " @ ", match_ind, " / ", length(meshset.matches))
  imview, matches, vectors, lines, params = inspect_matches(meshset, match_ind);
  indices_to_remove, break_review = edit_matches(imview..., matches, vectors, params);
  if !break_review
    store_points(path, meshset, k, indices_to_remove, username, "manual")
    match_ind += 1
    if match_ind > length(meshset.matches)
      meshset_ind += 1
      match_ind = 1
    end
    inspect_montages(username, meshset_ind, match_ind)
  end
end