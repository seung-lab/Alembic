"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp_mesh(mesh::Mesh)
  img = get_image(mesh)
  src_nodes, dst_nodes = get_globalized_nodes_h(mesh);
  src_nodes = src_nodes'
  dst_nodes = dst_nodes'
  offset = get_offset(mesh);
  node_dict = incidence_to_dict(mesh.edges')
  triangles = dict_to_triangles(node_dict)
  return @time ImageRegistration.meshwarp(img, src_nodes, dst_nodes, triangles, offset), mesh.index
end

"""
Multiple dispatch so Dodam doesn't have to type sooo much
"""
function render_montaged(index::Index; render_full=true, render_review=true)
  render_montaged(index, index; render_full=render_full, render_review=render_review)
end

"""
Cycle through indices and render montages
"""
function render_montaged(firstindex::Index, lastindex::Index; 
                                          render_full=true, render_review=true)
  firstindex = montaged(firstindex)
  lastindex = montaged(lastindex)
  for index in get_index_range(firstindex, lastindex)
    meshset = load(index)
    render_montaged(meshset; render_full=render_full, render_review=render_review)
  end 
end

function render_montaged(meshset::MeshSet; render_full=false, render_review=true, flagged_only=true)
  assert(is_premontaged(meshset.meshes[1].index))
  index = montaged(meshset.meshes[1].index)
  if is_flagged(meshset) 
    println("The meshset has a flag. Continuing anyway....")
  end

  try
    new_fn = get_filename(index)
    println("Rendering ", new_fn)
    warps = pmap(meshwarp_mesh, meshset.meshes);
    imgs = [x[1][1] for x in warps];
    offsets = [x[1][2] for x in warps];
    indices = [x[2] for x in warps];
    # review images
    if render_review
      write_seams(meshset, imgs, offsets, indices, flagged_only)
    end
    if render_full
      println(typeof(imgs))
      img, offset = merge_images(imgs, offsets)
      println("Writing ", new_fn)
      f = h5open(get_path(index), "w")
      chunksize = min(1000, min(size(img)...))
      @time f["img", "chunk", (chunksize,chunksize)] = img
      close(f)
      update_offset(index, [0,0], size(img))
    end
  catch e
    println(e)
    log_error(index; comment=e)
  end

end

"""
Calculate prealignment transforms from first section through section_num
"""
function calculate_cumulative_tform(index::Index, startindex=ROI_FIRST)
  index = montaged(index)
  startindex = montaged(startindex)
  cumulative_tform = eye(3)
  if index != startindex
    index_pairs = get_sequential_index_pairs(startindex, index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexB, indexA)
      # reset cumulative tform if the mesh is fixed
      if is_fixed(meshset.meshes[2])
        println(meshset.meshes[2].index, ": fixed")
        cumulative_tform = eye(3)
      end
      # tform = affine_approximate(meshset)
      offset = get_offset(indexB)
      translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
      tform = regularized_solve(meshset, lambda=0.9)
      cumulative_tform = cumulative_tform*translation*tform
    end
  end
  return cumulative_tform
end

function render_prealigned(index::Index; render_full=false, render_review=true)
  src_index = montaged(index)
  dst_index = get_preceding(src_index)
  meshset = load(src_index, dst_index)
  offset = get_offset(src_index)
  translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
  tform = translation*regularized_solve(meshset, lambda=0.9)
  render_prealigned(src_index, dst_index, tform; 
                          render_full=render_full, render_review=render_review)
end

function render_prealigned(src_index::Index, dst_index::Index, tform; 
                                          render_full=false, render_review=true)
  println("Loading images for rendering... 1/2")
  src_img = get_image(montaged(src_index))
  println("Loading images for rendering... 2/2")
  dst_img = get_image(montaged(dst_index))
  render_prealigned(src_index, dst_index, src_img, dst_img, tform; 
                          render_full=render_full, render_review=render_review)
end

function render_prealigned(src_index::Index, dst_index::Index, src_img, dst_img, 
                                  tform; render_full=false, render_review=true)
  scale = 0.05
  s = [scale 0 0; 0 scale 0; 0 0 1]

  if render_full
    src_index = prealigned(src_index)
    println("Warping images for rendering... 1/1")
    src_warped, src_offset = imwarp(src_img, tform, [0,0])
    update_offset(src_index, src_offset, size(src_warped))
    println("Writing full image:", get_filename(src_index))
    f = h5open(get_path(src_index), "w")
    chunksize = min(1000, min(size(src_warped)...))
    @time f["img", "chunk", (chunksize,chunksize)] = src_warped
    close(f)    
  end

  if render_review
    println("Warping images for rendering... 1/2")
    src_thumb, src_thumb_offset = imwarp(src_img, tform*s, [0,0])
    println("Warping images for rendering... 2/2")
    dst_thumb, dst_thumb_offset = imwarp(dst_img, s, [0,0])
    path = get_review_path(src_index, dst_index)
    O, O_bb = imfuse(dst_thumb, dst_thumb_offset, src_thumb, src_thumb_offset)
    f = h5open(path, "w")
    chunksize = min(1000, min(size(O)...))
    @time f["img", "chunk", (chunksize,chunksize)] = O
    f["offset"] = O_bb
    f["scale"] = scale
    f["tform"] = tform
    close(f)
    println("Writing review image:\n\t", path)
  end
end

# function render_prealigned(firstindex::Index, lastindex::Index; render_full=true, render_review=false)
#   firstindex = montaged(firstindex)
#   lastindex = montaged(lastindex)
#   dst_img = nothing
  
#   for index in get_index_range(firstindex, lastindex)
#     src_index = index
#     dst_index = get_preceding(src_index)
#     meshset = load(src_index, dst_index) 
#     println("Loading src_image for rendering...")
#     src_img = get_image(src_index)
#     if dst_img == nothing
#       println("Loading dst_image for rendering...")
#       dst_img = get_image(dst_index)
#     end

#     cumulative_tform = calculate_cumulative_tform(dst_index)
#     offset = get_offset(src_index)
#     translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
#     tform = translation*regularized_solve(meshset, lambda=0.9)

#     render_prealigned(src_index::Index, dst_index::Index, src_img, dst_img, 
#                                   tform; render_full=render_full, render_review=render_review)

#     println("Swapping src_image to dst_image")
#     dst_img = copy(src_img)
#   end
# end

"""
Render aligned images
"""
function render_aligned_review(firstindex::Index, lastindex::Index, start=1, finish=0)
  firstindex = prealigned(firstindex)
  lastindex = prealigned(lastindex)
  meshset = load(firstindex, lastindex)
  render_aligned_review(meshset, start, finish)
end

function render_aligned_review(meshset, start=1, finish=length(meshset.matches))
  scale = 0.10
  s = [scale 0 0; 0 scale 0; 0 0 1]
  images = Dict()
  BB = GLOBAL_BB
  
  # Check images dict for thumbnail, otherwise render it - just moving prealigned
  function retrieve_image(mesh)
    index = prealigned(mesh.index)
    if !(index in keys(images))
      println("Making review for ", mesh.index)
      # @time (img, offset), _ = meshwarp_mesh(mesh)
      # GLOBAL_BB = BoundingBox(-4000,-4000,38000,38000)
      img = get_image(mesh)
      offset = get_offset(mesh)
      # @time img = rescopeimage(img, offset, BB)
      img, offset = imwarp(img, s, offset)
      images[index] = img, offset
    end
    return images[index]
  end

  for (k, match) in enumerate(meshset.matches[start:finish])
    src_index = match.src_index
    dst_index = match.dst_index

    src_mesh = meshset.meshes[find_mesh_index(meshset, src_index)]
    dst_mesh = meshset.meshes[find_mesh_index(meshset, dst_index)]

    src_img, src_offset = retrieve_image(src_mesh)
    dst_img, dst_offset = retrieve_image(dst_mesh)
    # offset = [BB.i, BB.j] * scale
    O, O_bb = imfuse(src_img, src_offset, dst_img, dst_offset)

    indexA = prealigned(src_index)
    indexB = prealigned(dst_index)

    path = get_review_path(indexB, indexA)
    println("Writing thumbnail:\n\t", path)
    f = h5open(path, "w")
    @time f["img", "chunk", (1000,1000)] = O
    f["offset"] = O_bb # same as offset
    f["scale"] = scale
    close(f)
  end
end

"""
Render aligned images
"""
function render_aligned(waferA, secA, waferB, secB, start=1, finish=0)
  indexA = prealigned(waferA, secA)
  indexB = prealigned(waferB, secB)
  meshset = load(indexA, indexB)
  render_aligned(meshset, start, finish)
end

function render_aligned(meshset, start=1, finish=0)
  scale = 0.10
  s = [scale 0 0; 0 scale 0; 0 0 1]

  if start <= 0
    start = 1
  end
  if finish <= 0
    finish = length(meshset.meshes)
  end
  images = Dict()

  for (k, mesh) in enumerate(meshset.meshes[start:finish])
    index = aligned(mesh.index)
    println("Warping ", mesh.index)
    @time (img, offset), _ = meshwarp_mesh(mesh)
    println("Writing ", get_name(index))
    f = h5open(get_path(index), "w")
    @time f["img", "chunk", (1000,1000)] = img
    close(f)
    # Log image offsets
    update_offset(index, offset, size(img))
    images[index] = imwarp(img, s) 
    # Rescope the image & save
    write_finished(index, img, offset, GLOBAL_BB)
  end

  for (k, match) in enumerate(meshset.matches)
    src_index = aligned(match.src_index)
    dst_index = aligned(match.dst_index)

    if start <= find_mesh_index(meshset, src_index) <= finish &&
          start <= find_mesh_index(meshset, dst_index) <= finish

      src_img, src_offset = retrieve_image(src_mesh)
      dst_img, dst_offset = retrieve_image(dst_mesh)
      O, O_bb = imfuse(src_img, src_offset, dst_img, dst_offset)

      path = get_review_path(dst_index, src_index)
      println("Writing thumbnail:\n\t", path)
      f = h5open(path, "w")
      @time f["img", "chunk", (1000,1000)] = O
      f["offset"] = O_bb # same as offset
      f["scale"] = scale
      close(f)
    end
  end
end

function render_finished(waferA, secA, waferB, secB)
  indexA = aligned(waferA, secA)
  indexB = aligned(waferB, secB)
  for index in get_index_range(indexA, indexB)
    img = get_image(index)
    offset = get_offset(index)
    write_finished(index, img, offset)
  end
end

function write_finished(index, img, offset, BB=GLOBAL_BB)
  println("Rescoping ", get_name(index))
  @time img = rescopeimage(img, offset, BB)
  index = finished(index)
  println("Writing ", get_name(index))
  f = h5open(get_path(index), "w")
  @time f["img", "chunk", (1000,1000)] = img
  f["offset"] = offset
  f["bb"] = [BB.i, BB.j, BB.w, BB.h]
  close(f)
end


"""
Prealignment where offsets are global
"""
function render_prealigned(firstindex::Index, lastindex::Index; render_full=true, render_review=false)
  firstindex = montaged(firstindex)
  lastindex = montaged(lastindex)
  dir = PREALIGNED_DIR
  scale = 0.05
  s = [scale 0 0; 0 scale 0; 0 0 1]
  fixed = Dict()

  cumulative_tform = calculate_cumulative_tform(firstindex)
  # cumulative_tform = eye(3)

  # return Dictionary of staged image to remove redundancy in loading
  function stage_image(mesh, cumulative_tform, tform)
    stage = Dict()
    stage["index"] = montaged(mesh.index)
    img = get_image(mesh)
    println("tform:\n", tform)
    if cumulative_tform*tform == eye(3)
      stage["img"], stage["offset"] = img, [0,0]
    else
      println("Warping ", get_index(mesh))
      @time stage["img"], stage["offset"] = imwarp(img, cumulative_tform*tform, [0,0])
    end
    println("Creating thumbnail for ", get_index(mesh))
    stage["thumb_fixed"], stage["thumb_offset_fixed"] = imwarp(img, s, [0,0])
    stage["thumb_moving"], stage["thumb_offset_moving"] = imwarp(img, tform*s, [0,0])
    stage["scale"] = scale
    return stage
  end

  function save_image(stage)
    new_fn = string(join(stage["index"][1:2], ","), "_prealigned.h5")
    update_offset(prealigned(stage["index"]), stage["offset"], size(stage["img"]))
    println("Writing image:\n\t", new_fn)
    # @time imwrite(stage["img"], joinpath(dir, fn))
    f = h5open(joinpath(dir, new_fn), "w")
    chunksize = min(1000, min(size(stage["img"])...))
    @time f["img", "chunk", (chunksize,chunksize)] = stage["img"]
    close(f)
  end

  index_pairs = get_sequential_index_pairs(firstindex, lastindex)
  for (k, (indexA, indexB)) in enumerate(index_pairs)
    println("\nRendering ", indexA, " & ", indexB)
    meshset = load(indexB, indexA)
    if k==1
      fixed = stage_image(meshset.meshes[2], cumulative_tform, eye(3))
    end
    offset = get_offset(indexB)
    translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
    tform = translation*regularized_solve(meshset, lambda=0.9)
    moving = stage_image(meshset.meshes[1], cumulative_tform, tform)
    
    # save full scale image
    if render_full
      save_image(moving)
    end

    if render_review
      # save thumbnail of fused images
      path = get_review_path(moving["index"], fixed["index"])
      O, O_bb = imfuse(fixed["thumb_fixed"], fixed["thumb_offset_fixed"], 
                            moving["thumb_moving"], moving["thumb_offset_moving"])
      f = h5open(path, "w")
      chunksize = min(1000, min(size(O)...))
      @time f["img", "chunk", (chunksize,chunksize)] = O
      f["offset"] = O_bb
      f["scale"] = scale
      f["tform"] = tform
      close(f)
      println("Writing thumb:\n\t", path)
    end

    # propagate for the next section
    fixed = moving
    cumulative_tform = cumulative_tform*tform
  end
end