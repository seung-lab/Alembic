"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp_mesh(mesh::Mesh)
  img = get_image(mesh)
  src_nodes, dst_nodes = get_globalized_nodes(mesh);
  src_nodes = src_nodes'
  dst_nodes = dst_nodes'
  offset = get_offset(mesh) - get_topleft_offset(mesh)
  node_dict = incidence_to_dict(mesh.edges')
  triangles = dict_to_triangles(node_dict)
  return @time ImageRegistration.meshwarp(img, src_nodes, dst_nodes, triangles, offset), mesh.index
end

"""
Load image from hdf5, then apply meshwarp
"""
function meshwarp_h5(mesh::Mesh)
  @time img = get_h5_image(mesh)
  src_nodes = hcat(mesh.nodes...)'
  dst_nodes = hcat(mesh.nodes_t...)'
  offset = mesh.disp
  node_dict = incidence_to_dict(mesh.edges')
  triangles = dict_to_triangles(node_dict)
  return @time meshwarp(img, src_nodes, dst_nodes, triangles, offset)
end  

"""
Multiple dispatch so Dodam doesn't have to type sooo much
"""
function render_montaged(wafer_no, section_no, render_full=false)
  render_montaged(wafer_no, section_no, wafer_no, section_no, render_full)
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
    new_fn = string(idx[1], ",", idx[2], "_montaged.h5")
    println("Rendering ", new_fn)
    warps = pmap(meshwarp_mesh, meshset.meshes);
    imgs = [x[1][1] for x in warps];
    offsets = [x[1][2] for x in warps];
    indices = [x[2] for x in warps];
    # review images
    write_seams(meshset, imgs, offsets, indices)
    if render_full
      println(typeof(imgs))
      img, offset = merge_images(imgs, offsets)
      println("Writing ", new_fn)
      f = h5open(joinpath(MONTAGED_DIR, new_fn), "w")
      @time f["img", "chunk", (1000,1000)] = img
      close(f)
      update_offsets((index[1:2]...,-2,-2), [0,0], size(img))
    end
  end
end

"""
Calculate prealignment transforms from first section through section_num
"""
function calculate_cumulative_tform(index, dir=PREALIGNED_DIR)
  cumulative_tform = eye(3)
  if index != (1,1,-2,-2)
    index_pairs = get_sequential_index_pairs((1,1,-2,-2), index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexA, indexB)
      # tform = affine_approximate(meshset)
      offset = load_offset(indexB)
      translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
      tform = regularized_solve(meshset, lambda=0.9)
      cumulative_tform = cumulative_tform*translation*tform
    end
  end
  return cumulative_tform
end

"""
Scale & transform moving image for thumbnail using first matches type in meshset
"""
function write_prealignment_thumbnail(moving_img, fixed_img, meshset, scale=0.05)
  moving = Dict()
  fixed = Dict()
  s = [scale 0 0; 0 scale 0; 0 0 1]
  moving_offset = collect(meshset.meshes[2].disp)
  tform = regularized_solve(meshset, lambda=0.9)
  moving["nodes"], fixed["nodes"] = get_matched_points(meshset, 1)
  moving["nodes"] = transform_matches(moving["nodes"], tform)
  fixed["index"] = (meshset.meshes[1].index[1:2]..., -3, -3)
  moving["index"] = (meshset.meshes[2].index[1:2]..., -3, -3)
  fixed["thumb_fixed"], fixed["thumb_offset_fixed"] = imwarp(fixed_img, s)
  moving["thumb_moving"], moving["thumb_offset_moving"] = imwarp(moving_img, tform*s, moving_offset)
  fixed["scale"] = scale
  moving["scale"] = scale
  write_thumbnail_from_dict(fixed, moving)
end

"""
Fuse and save thumbnail images
"""
function write_thumbnail_from_dict(A, B)
  path = get_outline_filename("thumb", B["index"], A["index"])
  # path = string(path[1:end-4], ".png")
  O, O_bb = imfuse(A["thumb_fixed"], A["thumb_offset_fixed"], 
                            B["thumb_moving"], B["thumb_offset_moving"])
  vectorsA = scale_matches(A["nodes"], A["scale"])
  vectorsB = scale_matches(B["nodes"], B["scale"])
  vectors = (offset_matches(vectorsA, vectorsB, O_bb),)
  match_nums = (1,)
  colors = ([1,1,1],)
  factor = 20
  write_thumbnail(O, path, vectors, colors, match_nums, factor)
end

"""
Return Dictionary of staged image to remove redundancy in loading
"""
function stage_image(mesh, cumulative_tform, tform, scale=0.05)
  s = [scale 0 0; 0 scale 0; 0 0 1]
  stage = Dict()
  stage["index"] = (mesh.index[1:2]..., -3, -3)
  img = get_uint8_image(mesh)
  println("tform:\n", tform)
  println("Warping ", mesh.name)
  @time stage["img"], stage["offset"] = imwarp(img, cumulative_tform*tform, [0,0])
  println("Warping thumbnail for ", mesh.name)
  stage["thumb_fixed"], stage["thumb_offset_fixed"] = imwarp(img, s, [0,0])
  stage["thumb_moving"], stage["thumb_offset_moving"] = imwarp(img, tform*s, [0,0])
  stage["scale"] = scale
  return stage
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
function render_prealigned(waferA, secA, waferB, secB)
  indexA = (waferA, secA, -2, -2)
  indexB = (waferB, secB, -2, -2)
  dir = PREALIGNED_DIR
  fixed = Dict()

  cumulative_tform = calculate_cumulative_tform(indexA)
  log_path = joinpath(dir, "prealigned_offsets.txt")

  function save_image(stage, dir, log_path)
    new_fn = string(join(stage["index"][1:2], ","), "_prealigned.h5")
    update_offsets(stage["index"], stage["offset"], size(stage["img"]))
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
    meshset = load(indexA, indexB)
    if k==1
      fixed = stage_image(meshset.meshes[1], cumulative_tform, eye(3))
      if is_first_section(indexA)
        # save_image(fixed, dir, log_path)
      end
    end
    offset = load_offset(indexB)
    translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
    tform = regularized_solve(meshset, lambda=0.9)
    moving = stage_image(meshset.meshes[2], cumulative_tform, translation*tform)
    save_image(moving, dir, log_path)
    moving["nodes"], fixed["nodes"] = get_matched_points(meshset, 1)
    moving["nodes"] = transform_matches(moving["nodes"], tform)
    write_thumbnail_from_dict(fixed, moving)
    fixed = moving
    cumulative_tform = cumulative_tform*translation*tform
  end
end

function update_ext(fn, ext)
  i = searchindex(fn, ".")
  s = fn
  if i > 0
    s = string(fn[1:i], ext)
  end
  return s
end

"""
Cycle through JLD files in aligned directory and render alignment
"""
function render_aligned(waferA, secA, waferB, secB, start=1, finish=0)
  indexA = (waferA, secA, -3, -3)
  indexB = (waferB, secB, -3, -3)
  dir = ALIGNED_DIR
  scale = 0.10
  s = [scale 0 0; 0 scale 0; 0 0 1]

  # Log file for image offsets
  log_path = joinpath(dir, "aligned_offsets.txt")
  meshset = load_aligned(indexA, indexB)
  if start == 0
    start = 1
  end
  if finish == 0
    finish = length(meshset.meshes)
  end
  images = Dict()
  
  # Check images dict for thumbnail, otherwise render, save, & resize it
  function retrieve_image(mesh)
    index = (mesh.index[1:2]..., -4, -4)
    if !(index in keys(images))
      println("Warping ", mesh.name)
      @time img, offset = meshwarp_h5(mesh)
      @time img = rescopeimage(img, offset, GLOBAL_BB)
      new_fn = string(join(mesh.index[1:2], ","), "_aligned.h5")
      println("Writing ", new_fn)
      f = h5open(joinpath(dir, "round2", new_fn), "w")
      @time f["img", "chunk", (1000,1000)] = img
      close(f)
      # @time imwrite(img, joinpath(dir, new_fn))
      img, _ = imwarp(img, s)
      path = get_outline_filename("thumb", index)
      println("Writing thumbnail:\n\t", path)
      f = h5open(path, "w")
      @time f["img", "chunk", (1000,1000)] = img
      f["offset"] = [GLOBAL_BB.i, GLOBAL_BB.j] * scale
      f["scale"] = scale
      close(f)
      # path = update_ext(path, "tif")
      # @time Images.save(path, img)
      # Log image offsets
      update_offsets(index, offset, size(img))
      # end
      images[index] = img
    end
    return images[index]
  end

  indices = 1:length(meshset.matches)

  # for (k, matches) in zip(reverse(indices), reverse(meshset.matches))
  for (k, matches) in enumerate(meshset.matches)
    src_index = matches.src_index
    dst_index = matches.dst_index
    if start <= src_index[2] <= finish && start <= dst_index[2] <= finish
      src_mesh = meshset.meshes[find_index(meshset, src_index)]
      dst_mesh = meshset.meshes[find_index(meshset, dst_index)]

      src_nodes, dst_nodes = get_matched_points_t(meshset, k)
      vectorsA = scale_matches(src_nodes, scale)
      vectorsB = scale_matches(dst_nodes, scale)

      src_img = retrieve_image(src_mesh)
      dst_img = retrieve_image(dst_mesh)
      offset = [GLOBAL_BB.i, GLOBAL_BB.j] * scale
      O, O_bb = imfuse(src_img, offset, dst_img, offset)

      indexA = (src_index[1:2]..., -4, -4)
      indexB = (dst_index[1:2]..., -4, -4)

      path = get_outline_filename("thumb_imfuse", indexB, indexA)
      println("Writing thumbnail:\n\t", path)
      f = h5open(path, "w")
      @time f["img", "chunk", (1000,1000)] = O
      f["offset"] = O_bb # same as offset
      f["scale"] = scale
      close(f)

      # path = update_ext(get_outline_filename("thumb_pts", indexB, indexA), "tif")
      # vectors = (offset_matches(vectorsA, vectorsB, O_bb),)
      # match_nums = (find_matches_index(meshset, src_index, dst_index),)
      # colors = ([1,1,1],)
      # factor = 1
      # write_thumbnail(O, path, vectors, colors, match_nums, factor)
      # println("Writing thumbnail:\n", path)
      # f = h5open(path, "w")
      # @time f["img", "chunk", (1000,1000)] = drawing
      # f["offset"] = [GLOBAL_BB.i, GLOBAL_BB.j] * scale
      # f["scale"] = scale
      # close(f)
    end
  end
end
