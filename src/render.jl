# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150807
#
# Functions to load meshes, warp images via piecewise affine transforms, and
# display images with meshes.

"""
Multiple dispatch for imwarp on Mesh object
"""
function imwarp(meshset::MeshSet)
  # tform = recompute_affine(meshset)
  tform = affine_approximate(meshset)
  img = get_ufixed8_image(meshset.meshes[2])
  @time img, offset = imwarp(img, tform)
  return img, offset
end

"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp(mesh::Mesh)
  @time img = get_ufixed8_image(mesh)
  src_nodes = hcat(mesh.nodes...)'
  dst_nodes = hcat(mesh.nodes_t...)'
  offset = mesh.disp
  node_dict = incidence2dict(mesh.edges)
  triangles = dict2triangles(node_dict)
  return @time meshwarp(img, src_nodes, dst_nodes, triangles, offset)
end

function get_global_bb(meshset)
  bbs = []
  println("Calculating global bounding box")
  for mesh in meshset.meshes
      nodes = hcat(mesh.nodes_t...)'
      push!(bbs, snap_bb(find_mesh_bb(nodes)))
  end
  global_bb = sum(bbs)
  global_bb.h += 1
  global_bb.w += 1
  println(global_bb)
  return global_bb
end    

"""
Cycle through JLD files in montaged directory and render montage
"""
function render_montaged(section_range::Array{Int64})
  log_path = joinpath(MONTAGED_DIR, "montaged_offsets.txt")
  filenames = sort_dir(MONTAGED_DIR)[section_range]
  for filename in filenames
    println("Rendering ", filename[1:end-4])
    meshset = JLD.load(joinpath(MONTAGED_DIR, filename))["MeshSet"]
    warps = pmap(meshwarp, meshset.meshes)
    img, offset = merge_images([[x[i] for x in warps] for i=1:2]...)
    img = grayim(img)
    img["spatialorder"] = ["y","x"]
    println("Writing ", filename[1:end-4])
    new_filename = string(filename[1:end-4], ".tif")
    @time imwrite(img, joinpath(MONTAGED_DIR, new_filename))
    update_offset_log!(log_path, new_filename, [0,0], size(img))

    imfuse_section(meshset)
  end
end

"""
Write thumbnail image with vectors and match indices overlayed
"""
function write_thumbnail(path, img, vectors, factor, fontsize=24.0, offset=[-20,-20])
  imgc, img2 = view(img, pixelspacing=[1,1])
  a, b = draw_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1), factor)
  c = draw_indices(imgc, img2, vectors[1:2,:], fontsize, offset)
  println("Writing ", path)
  Cairo.write_to_png(imgc.c.back, path)
  destroy(toplevel(imgc))
end

"""
Calculate prealignment transforms from first section through section_num
"""
function calculate_global_tform(index, dir=PREALIGNED_DIR)
  global_tform = eye(3)
  if index != (1,1,-2,-2)
    index_pairs = create_sequential_index_pairs((1,1,-2,-2), index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexA, indexB)
      # tform = affine_approximate(meshset)
      offset = load_offset(indexB)
      translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
      tform = translation*regularized_solve(meshset, lambda=0.9)
      global_tform = global_tform*tform
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
  moving_nodes, fixed_nodes = get_matched_points(meshset, 1)
  fixed["index"] = meshset.meshes[1].index
  moving["index"] = meshset.meshes[2].index
  fixed["nodes"] = points_to_Nx3_matrix(fixed_nodes)
  moving["nodes"] = points_to_Nx3_matrix(moving_nodes)*tform
  fixed["thumb_fixed"], fixed["thumb_offset_fixed"] = imwarp(fixed_img, s)
  moving["thumb_moving"], moving["thumb_offset_moving"] = imwarp(moving_img, tform*s, moving_offset)
  fixed["scale"] = scale
  moving["scale"] = scale
  save_prealignment_thumbnails(fixed, moving)
end

"""
Return Dictionary of staged image to remove redundancy in loading
"""
function stage_image(mesh, global_tform, tform, scale=0.05)
  s = [scale 0 0; 0 scale 0; 0 0 1]
  stage = Dict()
  stage["index"] = mesh.index
  img = get_ufixed8_image(mesh)
  offset = load_offset(mesh.index)
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
Fuse and save thumbnail images
"""
function save_prealignment_thumbnails(A, B)
  dir = PREALIGNED_DIR
  fn = string(join(B["index"][1:2], ","), "_prealigned_thumbnail.png")
  println("Saving thumbnail ", fn)
  path = joinpath(dir, "review", fn)
  O, O_bb = imfuse(A["thumb_fixed"], A["thumb_offset_fixed"], 
                            B["thumb_moving"], B["thumb_offset_moving"])
  moving_nodes = A["nodes"][:,1:2]'*A["scale"]
  fixed_nodes = B["nodes"][:,1:2]'*B["scale"]
  moving_nodes .-= O_bb
  fixed_nodes .-= O_bb
  vectors = [moving_nodes; fixed_nodes]
  write_thumbnail(path, O, vectors, 1.0)
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
    fn = string(join(stage["index"][1:2], ","), "_prealigned.tif")
    update_offset_log!(log_path, fn, stage["offset"], size(stage["img"]))
    println("Writing ", fn)
    @time imwrite(stage["img"], joinpath(dir, fn))
  end

  index_pairs = create_sequential_index_pairs(indexA, indexB)
  for (k, (indexA, indexB)) in enumerate(index_pairs)
    println("\nPrealigning ", indexA, " & ", indexB)
    meshset = load(indexA, indexB)
    if k==1
      fixed = stage_image(meshset.meshes[1], global_tform, eye(3))
      if is_first_section(indexA)
        save_image(fixed, dir, log_path)
      end
    end
    offset = load_offset(indexB)
    translation = [1 0 0; 0 1 0; offset[1] offset[2] 1]
    tform = regularized_solve(meshset, lambda=0.9)
    # tform = affine_approximate(meshset)
    moving = stage_image(meshset.meshes[2], global_tform, tform)
    save_image(moving, dir, log_path)
    moving_nodes, fixed_nodes = get_matched_points(meshset, 1)
    fixed["nodes"] = points_to_Nx3_matrix(fixed_nodes)
    moving["nodes"] = points_to_Nx3_matrix(moving_nodes)*tform
    save_prealignment_thumbnails(fixed, moving)
    fixed = 0
    fixed = moving
    moving = 0
    global_tform = global_tform*translation*tform
  end
end

"""
Cycle through JLD files in aligned directory and render alignment
"""
function render_aligned(file_index, start=1, finish=0)
  dir = ALIGNED_DIR
  scale = 0.02
  s = [scale 0 0; 0 scale 0; 0 0 1]

  # Log file for image offsets
  log_path = joinpath(dir, "aligned_offsets.txt")

  filename = sort_dir(dir)[file_index]
  println("Rendering meshes in ", filename)
  meshset = JLD.load(joinpath(dir, filename))["MeshSet"]
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
      path = get_path(index)
      # if isfile(path)
      #   img, _ = imwarp(get_ufixed8_image(index), s)
      # else
      println("Warping ", mesh.name)
      @time img, offset = meshwarp(mesh)
      @time img = rescopeimage(img, offset, GLOBAL_BB)
      println("Writing ", mesh.name)
      new_fn = string(join(mesh.index[1:2], ","), "_aligned.h5")
      f = h5open(joinpath(dir, new_fn), "w")
      @time f["img", "chunk", (1000,1000)] = ufixed8_to_uint8(img)
      close(f)
      # @time imwrite(img, joinpath(dir, new_fn))
      img, _ = imwarp(img, s)

      # Log image offsets
      update_offset_log!(log_path, new_fn, offset, size(img))
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

      src_nodes = hcat(src_nodes...)[1:2, :]*scale .- src_offset
      dst_nodes = hcat(dst_nodes...)[1:2, :]*scale .- dst_offset
      vectors = [src_nodes; dst_nodes]
      thumbnail_fn = string(join(dst_index[1:2], ","), "-", join(src_index[1:2], ","), "_aligned_thumbnail.png")
      path = joinpath(dir, "review", thumbnail_fn)
      write_thumbnail(path, O, vectors, 0.5, 6.0, [-10,-10])
    end
  end
end

function write_alignment_blockmatches(section_range::Array{Int64})
  filenames = sort_dir(ALIGNED_DIR)[section_range]
  for filename in filenames
    println("Rendering meshes in ", filename)
    meshset = JLD.load(joinpath(ALIGNED_DIR, filename))["MeshSet"]
    save_blockmatch_imgs(meshset, k, [], joinpath(ALIGNED_DIR, "blockmatches"))
  end
end  
