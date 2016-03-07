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
function render_montaged(wafer_no, section_no; render_full=true, render_review=true)
  render_montaged(wafer_no, section_no, wafer_no, section_no, render_full=render_full, render_review=render_review)
end

function render_montaged_review(fn)
  meshset = load(joinpath(MONTAGED_DIR, fn))
  warps = pmap(meshwarp_mesh, meshset.meshes);
  imgs = [x[1][1] for x in warps];
  offsets = [x[1][2] for x in warps];
  indices = [x[2] for x in warps];
  # review images
  write_seams(meshset, imgs, offsets, indices, fn[1:end-4])
end

function render_montaged(meshset)
  try
    warps = pmap(meshwarp_mesh, meshset.meshes);
    imgs = [x[1][1] for x in warps];
    offsets = [x[1][2] for x in warps];
    indices = [x[2] for x in warps];
    # review images
    write_seams(meshset, imgs, offsets, indices)
  catch
    idx = (meshset.meshes[1].index[1:2]..., -2, -2)
    log_render_error(MONTAGED_DIR, idx, comment="")
  end
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
      path = get_review_filename(indices[i], indices[j])
      try 
        img, fuse_offset = imfuse(imgs[i], offsets[i], imgs[j], offsets[j])
        bb = bbs[i] - bbs[j]
        img_cropped = imcrop(img, fuse_offset, bb)
        f = h5open(path, "w")
        @time f["img", "chunk", (50,50)] = img_cropped
        f["offset"] = [bb.i, bb.j]
        f["scale"] = 1.0
        close(f)
      catch e
        idx = (indices[i], indices[j])
        log_render_error(MONTAGED_DIR, idx, e)
      end
    end
end

"""
Cycle through JLD files in montaged directory and render montage
"""
function render_montaged(waferA, secA, waferB, secB; render_full=true, render_review=true)
  indexA = (waferA, secA, -2, -2)
  indexB = (waferB, secB, -2, -2)
  for index in get_index_range(indexA, indexB)
    idx = (index[1:2]..., 1, 1)
    meshset = load(idx, idx)
    try
      new_fn = string(idx[1], ",", idx[2], "_montaged.h5")
      println("Rendering ", new_fn)
      warps = pmap(meshwarp_mesh, meshset.meshes);
      imgs = [x[1][1] for x in warps];
      offsets = [x[1][2] for x in warps];
      indices = [x[2] for x in warps];
      # review images
      if render_review
        write_seams(meshset, imgs, offsets, indices)
      end
      if render_full
        println(typeof(imgs))
        img, offset = merge_images(imgs, offsets)
        println("Writing ", new_fn)
        f = h5open(joinpath(MONTAGED_DIR, new_fn), "w")
        @time f["img", "chunk", (1000,1000)] = img
        close(f)
        update_offset((index[1:2]...,-2,-2), [0,0], size(img))
      end
    catch e
      idx = (meshset.meshes[1].index[1:2]..., -2, -2)
      log_render_error(MONTAGED_DIR, idx, e)
    end
  end 
end

"""
Calculate prealignment transforms from first section through section_num
"""
function calculate_cumulative_tform(index, dir=PREALIGNED_DIR)
  cumulative_tform = eye(3)
  if index != (1,1,-2,-2)
    index_pairs = get_sequential_index_pairs((2,149,-2,-2), index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexB, indexA)
      # reset cumulative tform is the mesh is fixed
      if haskey(meshset.meshes[2].properties, "fixed")
        if meshset.meshes[2].properties["fixed"]
          println(meshset.meshes[2].index, ": fixed")
          cumulative_tform = eye(3)
        end
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

"""
Copy a section from one process step to the next
"""
function copy_section_through(index)
  println("copy_section_through INCOMPLETE")
  return
end

"""
Prealignment where offsets are global
"""
function render_prealigned(waferA, secA, waferB, secB; render_full=true, render_review=true)
  indexA = (waferA, secA, -2, -2)
  indexB = (waferB, secB, -2, -2)
  dir = PREALIGNED_DIR
  scale = 0.05
  s = [scale 0 0; 0 scale 0; 0 0 1]
  fixed = Dict()

  cumulative_tform = calculate_cumulative_tform(indexA)
  # cumulative_tform = eye(3)

  # return Dictionary of staged image to remove redundancy in loading
  function stage_image(mesh, cumulative_tform, tform)
    stage = Dict()
    stage["index"] = (mesh.index[1:2]..., -2, -2)
    img = get_image(mesh)
    println("tform:\n", tform)
    println("Warping ", get_index(mesh))
    @time stage["img"], stage["offset"] = imwarp(img, cumulative_tform*tform, [0,0])
    println("Warping thumbnail for ", get_index(mesh))
    stage["thumb_fixed"], stage["thumb_offset_fixed"] = imwarp(img, s, [0,0])
    stage["thumb_moving"], stage["thumb_offset_moving"] = imwarp(img, tform*s, [0,0])
    stage["scale"] = scale
    return stage
  end

  function save_image(stage, dir)
    new_fn = string(join(stage["index"][1:2], ","), "_prealigned.h5")
    update_offset(prealigned(stage["index"]), stage["offset"], size(stage["img"]))
    println("Writing image:\n\t", new_fn)
    # @time imwrite(stage["img"], joinpath(dir, fn))
    f = h5open(joinpath(dir, new_fn), "w")
    chunksize = min(1000, min(size(stage["img"])...))
    @time f["img", "chunk", (chunksize,chunksize)] = stage["img"]
    close(f)
  end

  println("Cumulative tform:\n", cumulative_tform)
  index_pairs = get_sequential_index_pairs(indexA, indexB)
  for (k, (indexA, indexB)) in enumerate(index_pairs)
    println("\nPrealigning ", indexA, " & ", indexB)
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
      save_image(moving, dir)
    end

    if render_review
      # save thumbnail of fused images
      path = get_review_filename(moving["index"], fixed["index"])
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

"""
Render aligned images
"""
function render_aligned_review(waferA, secA, waferB, secB, start=1, finish=0)
  indexA = (waferA, secA, -3, -3)
  indexB = (waferB, secB, -3, -3)
  meshset = load(indexA, indexB)
  render_aligned_review(meshset, start, finish)
end

function render_aligned_review(meshset, start=1, finish=0)
  scale = 0.10
  s = [scale 0 0; 0 scale 0; 0 0 1]

  if start <= 0
    start = 1
  end
  if finish <= 0
    finish = length(meshset.matches)
  end
  images = Dict()
  BB = BoundingBox(-4000,-4000,46000,46000)
  
  # Check images dict for thumbnail, otherwise render it - just moving prealigned
  function retrieve_image(mesh)
    index = (mesh.index[1:2]..., -4, -4)
    if !(index in keys(images))
      println("Making review for ", mesh.index)
      # @time (img, offset), _ = meshwarp_mesh(mesh)
      # GLOBAL_BB = BoundingBox(-4000,-4000,38000,38000)
      img = get_image(mesh)
      offset = get_offset(mesh)
      @time img = rescopeimage(img, offset, BB)
      img, _ = imwarp(img, s)
      images[index] = img
    end
    return images[index]
  end

  for (k, match) in enumerate(meshset.matches[start:finish])
    src_index = match.src_index
    dst_index = match.dst_index

    src_mesh = meshset.meshes[find_mesh_index(meshset, src_index)]
    dst_mesh = meshset.meshes[find_mesh_index(meshset, dst_index)]

    src_img = retrieve_image(src_mesh)
    dst_img = retrieve_image(dst_mesh)
    offset = [BB.i, BB.j] * scale
    O, O_bb = imfuse(src_img, offset, dst_img, offset)

    indexA = (src_index[1:2]..., -4, -4)
    indexB = (dst_index[1:2]..., -4, -4)

    path = get_review_filename(indexB, indexA)
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
  indexA = (waferA, secA, -3, -3)
  indexB = (waferB, secB, -3, -3)
  meshset = load(indexA, indexB)
  render_aligned(meshset, start, finish)
end

function render_aligned(meshset, start=1, finish=0)
  dir = ALIGNED_DIR
  scale = 0.10
  s = [scale 0 0; 0 scale 0; 0 0 1]

  if start <= 0
    start = 1
  end
  if finish <= 0
    finish = length(meshset.meshes)
  end
  images = Dict()
  BB = BoundingBox(-4000,-4000,42000,42000)

  for (k, mesh) in enumerate(meshset.meshes[start:finish])
    index = (mesh.index[1:2]..., -4, -4)
    println("Warping ", mesh.index)
    @time (img, offset), _ = meshwarp_mesh(mesh)
    @time img = rescopeimage(img, offset, BB)
    new_fn = string(join(mesh.index[1:2], ","), "_aligned.h5")
    println("Writing ", new_fn)
    f = h5open(joinpath(dir, new_fn), "w")
    @time f["img", "chunk", (1000,1000)] = img
    close(f)
    # Log image offsets
    update_offset(index, offset, size(img))

    img, _ = imwarp(img, s)
    images[index] = img    
  end

  for (k, match) in enumerate(meshset.matches)
    src_index = match.src_index
    dst_index = match.dst_index

    if start <= find_mesh_index(meshset, src_index) <= finish &&
          start <= find_mesh_index(meshset, dst_index) <= finish

      src_img = images[src_index]
      dst_img = images[dst_index]

      offset = [BB.i, BB.j] * scale
      O, O_bb = imfuse(src_img, offset, dst_img, offset)

      indexA = (src_index[1:2]..., -4, -4)
      indexB = (dst_index[1:2]..., -4, -4)

      path = get_review_filename(indexB, indexA)
      println("Writing thumbnail:\n\t", path)
      f = h5open(path, "w")
      @time f["img", "chunk", (1000,1000)] = O
      f["offset"] = O_bb # same as offset
      f["scale"] = scale
      close(f)
    end
  end
end

"""
Write any render errors to a log file
"""
function log_render_error(dir, idx, comment="")
  ts = parse(Dates.format(now(), "yymmddHHMMSS"))
  path = joinpath(dir, "render_error_log.txt")
  new_row = [ts, idx, comment]'
  if !isfile(path)
    f = open(path, "w")
    close(f)
    log = new_row
  else  
    log = readdlm(path)
    log = vcat(log, new_row)
  end
  log = log[sortperm(log[:, 1]), :]
  println("Logging render error:\n", path)
  writedlm(path, log)
end
