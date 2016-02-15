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
function render_montaged(wafer_no, section_no, render_full=false)
  render_montaged(wafer_no, section_no, wafer_no, section_no, render_full)
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
Cycle through JLD files in montaged directory and render montage
"""
function render_montaged(waferA, secA, waferB, secB, render_full=false)
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
      write_seams(meshset, imgs, offsets, indices)
      # write_seams_with_points(meshset, imgs, offsets, indices)
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
    index_pairs = get_sequential_index_pairs((1,167,-2,-2), index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexB, indexA)
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
function render_prealigned(waferA, secA, waferB, secB, render_full=false)
  indexA = (waferA, secA, -2, -2)
  indexB = (waferB, secB, -2, -2)
  dir = PREALIGNED_DIR
  scale = 0.05
  s = [scale 0 0; 0 scale 0; 0 0 1]
  fixed = Dict()

  cumulative_tform = calculate_cumulative_tform(indexA)
  # cumulative_tform = eye(3)
  log_path = joinpath(dir, "prealigned_offsets.txt")

  # return Dictionary of staged image to remove redundancy in loading
  function stage_image(mesh, cumulative_tform, tform)
    stage = Dict()
    stage["index"] = (mesh.index[1:2]..., -3, -3)
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

  function save_image(stage, dir, log_path)
    new_fn = string(join(stage["index"][1:2], ","), "_prealigned.h5")
    update_offset(stage["index"], stage["offset"], size(stage["img"]))
    println("Writing image:\n\t", new_fn)
    # @time imwrite(stage["img"], joinpath(dir, fn))
    f = h5open(joinpath(dir, new_fn), "w")
    @time f["img", "chunk", (1000,1000)] = stage["img"]
    close(f)
  end

  println("Cumulative tform:\n", cumulative_tform)
  index_pairs = get_sequential_index_pairs(indexA, indexB)
  for (k, (indexA, indexB)) in enumerate(index_pairs)
    println("\nPrealigning ", indexA, " & ", indexB)
    meshset = load(indexB, indexA)
    if k==1
      fixed = stage_image(meshset.meshes[2], cumulative_tform, eye(3))
      if is_first_section(indexA)
        # save_image(fixed, dir, log_path)
      end
    end
    offset = get_offset(indexB)
    translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
    tform = translation*regularized_solve(meshset, lambda=0.9)
    moving = stage_image(meshset.meshes[1], cumulative_tform, tform)
    
    # save full scale image
    if render_full
      save_image(moving, dir, log_path)
    end

    # save thumbnail of fused images
    path = get_review_filename("thumb", moving["index"], fixed["index"])
    O, O_bb = imfuse(fixed["thumb_fixed"], fixed["thumb_offset_fixed"], 
                          moving["thumb_moving"], moving["thumb_offset_moving"])
    f = h5open(path, "w")
    @time f["img", "chunk", (1000,1000)] = O
    f["offset"] = O_bb
    f["scale"] = scale
    f["tform"] = tform
    close(f)

    # propagate for the next section
    fixed = moving
    cumulative_tform = cumulative_tform*tform
  end
end

"""
Cycle through JLD files in aligned directory and render alignment
"""
function render_aligned(waferA, secA, waferB, secB, render_full=false, start=1, finish=0)
  indexA = (waferA, secA, -3, -3)
  indexB = (waferB, secB, -3, -3)
  meshset = load(indexA, indexB)
  render_aligned(meshset, render_full, start, finish)
end

function render_aligned(meshset, render_full=false, start=1, finish=0)
  dir = ALIGNED_DIR
  scale = 0.10
  s = [scale 0 0; 0 scale 0; 0 0 1]

  # Log file for image offsets
  log_path = joinpath(dir, "aligned_offsets.txt")
  if start <= 0
    start = 1
  end
  if finish <= 0
    finish = length(meshset.meshes)
  end
  images = Dict()
  
  # Check images dict for thumbnail, otherwise render, save, & resize it
  function retrieve_image(mesh)
    index = (mesh.index[1:2]..., -4, -4)
    if !(index in keys(images))
      println("Warping ", mesh.index)
      if render_full
        @time (img, offset), _ = meshwarp_mesh(mesh)
        @time img = rescopeimage(img, offset, GLOBAL_BB)
        new_fn = string(join(mesh.index[1:2], ","), "_aligned.h5")
        println("Writing ", new_fn)
        f = h5open(joinpath(dir, new_fn), "w")
        @time f["img", "chunk", (1000,1000)] = img
        close(f)
        # Log image offsets
        update_offset(index, offset, size(img))
      end

      # 
      img = get_image(mesh)
      offset = get_offset(mesh)
      @time img = rescopeimage(img, offset, GLOBAL_BB)
      img, _ = imwarp(img, s)
      # path = get_review_filename("thumb", index)
      # println("Writing thumbnail:\n\t", path)
      # f = h5open(path, "w")
      # @time f["img", "chunk", (1000,1000)] = img
      # f["offset"] = [GLOBAL_BB.i, GLOBAL_BB.j] * scale
      # f["scale"] = scale
      # close(f)
      images[index] = img
    end
    return images[index]
  end

  indices = 1:length(meshset.matches)

  # for (k, matches) in zip(reverse(indices), reverse(meshset.matches))
  for (k, matches) in enumerate(meshset.matches[start:finish])
    src_index = matches.src_index
    dst_index = matches.dst_index
    # if start <= src_index[2] <= finish && start <= dst_index[2] <= finish
    src_mesh = meshset.meshes[find_mesh_index(meshset, src_index)]
    dst_mesh = meshset.meshes[find_mesh_index(meshset, dst_index)]

    src_img = retrieve_image(src_mesh)
    dst_img = retrieve_image(dst_mesh)
    offset = [GLOBAL_BB.i, GLOBAL_BB.j] * scale
    O, O_bb = imfuse(src_img, offset, dst_img, offset)

    indexA = (src_index[1:2]..., -4, -4)
    indexB = (dst_index[1:2]..., -4, -4)

    path = get_review_filename("thumb_imfuse", indexB, indexA)
    println("Writing thumbnail:\n\t", path)
    f = h5open(path, "w")
    @time f["img", "chunk", (1000,1000)] = O
    f["offset"] = O_bb # same as offset
    f["scale"] = scale
    close(f)

    # end
  end
end

"""
Write any render errors to a log file
"""
function log_render_error(dir, idx, comment="")
  ts = Dates.format(now(), "yymmddHHMMSS")
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