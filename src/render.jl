"""
Multiple dispatch for imwarp on Mesh object
"""
function imwarp(meshset::MeshSet)
  # tform = recompute_affine(meshset)
  tform = affine_approximate(meshset)
  img = get_uint8_image(meshset.meshes[2])
  @time img, offset = imwarp(img, tform)
  return img, offset
end

"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp(mesh::Mesh)
  @time img = get_uint8_image(mesh)
  src_nodes = hcat(mesh.nodes...)'
  dst_nodes = hcat(mesh.nodes_t...)'
  offset = mesh.disp
  node_dict = incidence2dict(mesh.edges)
  triangles = dict2triangles(node_dict)
  return @time meshwarp(img, src_nodes, dst_nodes, triangles, offset), mesh.index
end

"""
Load image from hdf5, then apply meshwarp
"""
function meshwarp_h5(mesh::Mesh)
  @time img = get_h5_image(mesh)
  src_nodes = hcat(mesh.nodes...)'
  dst_nodes = hcat(mesh.nodes_t...)'
  offset = mesh.disp
  node_dict = incidence2dict(mesh.edges)
  triangles = dict2triangles(node_dict)
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
    new_fn = string(idx[1], ",", idx[2], "_montaged.tif")
    println("Rendering ", new_fn)
    warps = pmap(meshwarp, meshset.meshes);
    imgs = [x[1][1] for x in warps];
    offsets = [x[1][2] for x in warps];
    indices = [x[2] for x in warps];
    # review images
    # write_seams_with_matches(meshset, imgs, offsets, indices)
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
function calculate_global_tform(index, dir=PREALIGNED_DIR)
  global_tform = eye(3)
  if index != (1,1,-2,-2)
    index_pairs = get_sequential_index_pairs((1,1,-2,-2), index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexA, indexB)
      # tform = affine_approximate(meshset)
      offset = load_offset(indexB)
      translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
      tform = regularized_solve(meshset, lambda=0.9)
      global_tform = global_tform*translation*tform
    end
  end
  return global_tform
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
  fixed["nodes"], moving["nodes"] = get_matched_points_t(meshset, 1)
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
  path = get_outline_filename(B["index"], A["index"], "thumb")
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
Write thumbnail image, with whatever drawings included
"""
function write_imageview(path, imgc, img2)
  println("Writing ", path)
  write_canvas(imgc, path)
  close_image(imgc)
end


"""
Return Dictionary of staged image to remove redundancy in loading
"""
function stage_image(mesh, global_tform, tform, scale=0.05)
  s = [scale 0 0; 0 scale 0; 0 0 1]
  stage = Dict()
  stage["index"] = (mesh.index[1:2]..., -3, -3)
  img = get_uint8_image(mesh)
  offset = load_offset(stage["index"])
  println("tform: ", tform)
  println("montaged_offset: ", offset)
  println("Warping ", mesh.name)
  translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
  @time stage["img"], stage["offset"] = imwarp(img, global_tform*translation*tform, [0,0])
  println("Warping thumbnail for ", mesh.name)
  stage["thumb_fixed"], stage["thumb_offset_fixed"] = imwarp(img, s, [0,0])
  stage["thumb_moving"], stage["thumb_offset_moving"] = imwarp(img, translation*tform*s, [0,0])
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
function render_prealigned(indexA, indexB)
  dir = PREALIGNED_DIR
  fixed = Dict()

  global_tform = calculate_global_tform(indexA)
  log_path = joinpath(dir, "prealigned_offsets.txt")

  function save_image(stage, dir, log_path)
    new_fn = string(join(stage["index"][1:2], ","), "_prealigned.h5")
    update_offsets(stage["index"], stage["offset"], size(stage["img"]))
    println("Writing ", new_fn)
    # @time imwrite(stage["img"], joinpath(dir, fn))
    f = h5open(joinpath(dir, new_fn), "w")
    @time f["img", "chunk", (1000,1000)] = stage["img"]
    close(f)
  end

  index_pairs = get_sequential_index_pairs(indexA, indexB)
  for (k, (indexA, indexB)) in enumerate(index_pairs)
    println("\nPrealigning ", indexA, " & ", indexB)
    meshset = load(indexA, indexB)
    if k==1
      fixed = stage_image(meshset.meshes[1], global_tform, eye(3))
      if is_first_section(indexA)
        # save_image(fixed, dir, log_path)
      end
    end
    offset = load_offset(indexB)
    translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
    tform = regularized_solve(meshset, lambda=0.9)
    # tform = affine_approximate(meshset)
    moving = stage_image(meshset.meshes[2], global_tform, tform)
    save_image(moving, dir, log_path)
    moving["nodes"], fixed["nodes"] = get_matched_points(meshset, 1)
    moving["nodes"] = transform_matches(moving["nodes"], tform)
   # write_thumbnail_from_dict(fixed, moving)
    fixed = moving
    global_tform = global_tform*translation*tform
  end
end

"""
Cycle through JLD files in aligned directory and render alignment
"""
function render_aligned(indexA, indexB, start=1, finish=0)
  dir = ALIGNED_DIR
  scale = 0.02
  s = [scale 0 0; 0 scale 0; 0 0 1]

  # Log file for image offsets
  log_path = joinpath(dir, "aligned_offsets.txt")
  meshset = load(indexA, indexB)
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
    if index in keys(images)
      img = images[index]
    else
      println("Warping ", mesh.name)
      @time img, offset = meshwarp_h5(mesh)
      @time img = rescopeimage(img, offset, GLOBAL_BB)
      println("Writing ", mesh.name)
      new_fn = string(join(mesh.index[1:2], ","), "_aligned.h5")
      f = h5open(joinpath(dir, new_fn), "w")
      @time f["img", "chunk", (1000,1000)] = img
      close(f)
      # @time imwrite(img, joinpath(dir, new_fn))
      img, _ = imwarp(img, s)

      # Log image offsets
      update_offsets(index, offset, size(img))
      # end
      images[index] = img
    end
    return img
  end

  # map(warp_pad_write, meshset.meshes)
  for (k, matches) in enumerate(meshset.matches)
    src_index = matches.src_index
    dst_index = matches.dst_index
    if start <= src_index[2] <= finish && start <= dst_index[2] <= finish
      src_mesh = meshset.meshes[find_index(meshset, src_index)]
      dst_mesh = meshset.meshes[find_index(meshset, dst_index)]

      src_nodes, dst_nodes = get_matched_points_t(meshset, k)
      src_index = (src_index[1:2]..., src_index[3]-1, src_index[4]-1)
      dst_index = (dst_index[1:2]..., dst_index[3]-1, dst_index[4]-1)
      src_offset = [GLOBAL_BB.i, GLOBAL_BB.j]
      dst_offset = [GLOBAL_BB.i, GLOBAL_BB.j]

      src_img = retrieve_image(src_mesh)
      dst_img = retrieve_image(dst_mesh)

      src_offset *= scale
      dst_offset *= scale

      O, O_bb = imfuse(src_img, src_offset, dst_img, dst_offset)

      # src_nodes = hcat(src_nodes...)[1:2, :]*scale .- src_offset
      # dst_nodes = hcat(dst_nodes...)[1:2, :]*scale .- dst_offset
      # vectors = [src_nodes; dst_nodes]
      path = get_thumbnail_path(dst_index, src_index)
      write_imageview(path, view_isotropic(O))
    end
  end
end
