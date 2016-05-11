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
  return @time meshwarp(img, src_nodes, dst_nodes, triangles, offset), get_index(mesh)
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

  # try
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
  # catch e
  #   println(e)
  #   log_error(index; comment=e)
  # end

end

"""
Calculate prealignment transforms from first section through section_num

Notes on transform composition:
* Matrix operations happen from right to left, so repeat the orders of tform
  calculations.
    1. previous transforms: cumulative_tform
    2. monoblock_match: translation
    3. montage to montage matching: tform

    tform * translation * cumulative_tform

  * Be aware that aligned images are already positioned in global space, but
    they are not *rescoped* to it (that's the finished images directory). So the
    aligned image offset needs to be accounted for as an additional translation.
    If the image is fixed, it's assumed to be an aligned image, so we pull its
    offset, and calculate the additional translation.
"""
function prepare_prealignment(index::Index, startindex=montaged(ROI_FIRST))
  src_index = montaged(index)
  dst_index = get_preceding(src_index)

  cumulative_tform = eye(3)
  tform = eye(3)
  for index in get_index_range(startindex, src_index)[2:end]
    cumulative_tform = tform*cumulative_tform
    src_index = index
    dst_index = get_preceding(src_index)
    meshset = load(src_index, dst_index)
    dst_index = get_index(meshset.meshes[2])
    src_offset = get_offset(src_index)
    translation = make_translation_matrix(src_offset)
    if is_fixed(meshset.meshes[2])
      println("FIXED")
      cumulative_tform = eye(3)
      dst_offset = get_offset(dst_index)
      translation = make_translation_matrix(dst_offset)*translation
    end
    tform = regularized_solve(meshset, lambda=0.9)*translation
  end
  return src_index, dst_index, cumulative_tform, tform
end

function render_prealigned(index::Index; render_full=false, render_review=true, startindex=montaged(ROI_FIRST))
  src_index, dst_index, cumulative_tform, tform = prepare_prealignment(index, startindex)
  render_prealigned(src_index, dst_index, cumulative_tform, tform; 
                          render_full=render_full, render_review=render_review)
end

function render_prealigned(src_index::Index, dst_index::Index, cumulative_tform, 
                                  tform; render_full=false, render_review=true)
  println("Loading images for rendering... 1/2")
  src_img = get_image(src_index)
  println("Loading images for rendering... 2/2")
  dst_img = get_image(dst_index)
  render_prealigned(src_index, dst_index, src_img, dst_img, cumulative_tform, 
                    tform; render_full=render_full, render_review=render_review)
end

function render_prealigned(firstindex::Index, lastindex::Index; 
                                        render_full=true, render_review=false, startindex=montaged(ROI_FIRST), align=false)
  startindex = montaged(startindex)
  dst_img = nothing
  for index in get_index_range(montaged(firstindex), montaged(lastindex))
    src_index, dst_index, cumulative_tform, tform = prepare_prealignment(index, startindex)
    println("Loading src_image for rendering")
    src_img = get_image(src_index)
    if dst_img == nothing
      println("Loading dst_image for rendering")
      dst_img = get_image(dst_index)
    end
    render_prealigned(src_index, dst_index, src_img, dst_img, cumulative_tform, 
                    tform; render_full=render_full, render_review=render_review)
    println("Swapping src_image to dst_image")
    dst_img = copy(src_img)

    println(index, " ", firstindex, " ", align && montaged(index) > montaged(firstindex))
    if align && montaged(index) > montaged(firstindex)
      println("Aligning meshes between $dst_index, $src_index")
      reload_registry(prealigned(src_index))
      ms = MeshSet(prealigned(dst_index), prealigned(src_index); solve=false, fix_first=(dst_index==startindex))
      render_aligned_review(ms)
    end
  end
end

function render_prealigned(src_index::Index, dst_index::Index, src_img, dst_img, 
                cumulative_tform, tform; render_full=false, render_review=true)
  scale = 0.05
  s = make_scale_matrix(scale)

  if render_review
    dst_offset = [0,0]
    if is_aligned(dst_index)
      dst_offset = get_offset(dst_index)
      println("dst image is aligned, so translate:\t$dst_offset")
    end
    println("Warping prealigned review image... 1/2")
    src_thumb, src_thumb_offset = imwarp(src_img, tform*s, [0,0])
    println("Warping prealigned review image... 2/2")
    dst_thumb, dst_thumb_offset = imwarp(dst_img, s, dst_offset)
    path = get_review_path(src_index, dst_index)
    write_review_image(path, src_thumb, src_thumb_offset, dst_thumb, dst_thumb_offset, scale, tform)

    println("Warping aligned review image... 1/2")
    src_thumb, src_thumb_offset = imwarp(src_img, tform*cumulative_tform*s, [0,0])
    println("Warping aligned review image... 2/2")
    dst_thumb, dst_thumb_offset = imwarp(dst_img, cumulative_tform*s, [dst_offset])
    aligned_path = get_review_path(prealigned(src_index), prealigned(dst_index))
    write_review_image(aligned_path, src_thumb, src_thumb_offset, dst_thumb, dst_thumb_offset, scale, tform*cumulative_tform)
  end

  if render_full
    src_index = prealigned(src_index)
    println("Warping full image... 1/1")
    @time src_warped, src_offset = imwarp(src_img, tform*cumulative_tform, [0,0])
    update_offset(src_index, src_offset, size(src_warped))
    path = get_path(src_index)
    println("Writing full image:\n ", path)
    f = h5open(path, "w")
    chunksize = min(1000, min(size(src_warped)...))
    @time f["img", "chunk", (chunksize,chunksize)] = src_warped
    close(f)    
  end
end

function write_review_image(path, src_img, src_offset, dst_img, dst_offset, scale, tform)
  O, O_bb = imfuse(dst_img, dst_offset, src_img, src_offset) # dst - red, src - green
  println("Writing review image:\n ", path)
  f = h5open(path, "w")
  chunksize = min(1000, min(size(O)...))
  @time f["img", "chunk", (chunksize,chunksize)] = O
  f["offset"] = O_bb
  f["scale"] = scale
  f["tform"] = tform
  close(f)
end

"""
Check images dict for thumbnail, otherwise render it - just moving prealigned
"""
function retrieve_image(images, index; tform=eye(3))
  if !(index in keys(images))
    println("Making review for ", index)
    img = get_image(index)
    offset = get_offset(index)
    img, offset = imwarp(img, tform, offset)
    images[index] = img, offset
  end
  return images[index]
end

"""
Render aligned images
"""
function render_aligned_review(firstindex::Index, lastindex::Index, start=1, finish=0)
  firstindex, lastindex = prealigned(firstindex), prealigned(lastindex)
  meshset = load(firstindex, lastindex)
  render_aligned_review(meshset, start, finish)
end

function render_aligned_review(meshset, start=1, finish=length(meshset.matches); images=Dict())
  scale = 0.05
  s = make_scale_matrix(scale)

  for (k, match) in enumerate(meshset.matches[start:finish])
    src_index = get_src_index(match)
    dst_index = get_dst_index(match)

    src_img, src_offset = retrieve_image(images, src_index; tform=s)
    dst_img, dst_offset = retrieve_image(images, dst_index; tform=s)
    path = get_review_path(prealigned(src_index), prealigned(dst_index))
    write_review_image(path, src_img, src_offset, dst_img, dst_offset, scale, s)
  end
end

"""
Render aligned images
"""
function render_aligned(firstindex::Index, lastindex::Index, start=1, finish=0)
  firstindex, lastindex = prealigned(firstindex), prealigned(lastindex)
  meshset = load(firstindex, lastindex)
  render_aligned(meshset, start, finish)
end

@fastmath @inbounds function render_aligned(meshset::MeshSet, start=1, finish=length(meshset.meshes))
  scale = 0.05
  s = make_scale_matrix(scale)
  images = Dict()
  for mesh_ind in start:finish
    if mesh_ind != finish
      fetch = prefetch(get_index(meshset.meshes[mesh_ind + 1]));
    end
    mesh = meshset.meshes[mesh_ind];
    index = aligned(mesh.index)
    println("Warping ", mesh.index)
    @time (img, offset), _ = meshwarp_mesh(mesh)
    println("Writing ", get_name(index))
    f = h5open(get_path(index), "w")
    @time f["img", "chunk", (1000,1000)] = img
    f["dtype"] = string(typeof(img[1]))
    f["offset"] = offset
    f["size"] = [size(img)...]
    close(f)
    # Log image offsets
    update_offset(index, offset, size(img))
    #images[index] = imwarp(img, s) 
    # Rescope the image & save
    #write_finished(index, img, offset, GLOBAL_BB)
    if mesh_ind != finish
      wait(fetch);
    end
  end
  # render_aligned_review(meshset; images=images)
end

# function render_finished(firstindex::Index, lastindex::Index)
#   for index in get_index_range(aligned(firstindex), aligned(lastindex))
#     img = get_image(index)
#     offset = get_offset(index)
#     write_finished(index, img, offset)
#   end
# end

# @fastmath @inbounds function write_finished(index::Index, img, offset, BB=GLOBAL_BB)
#   println("Rescoping ", get_name(index))
#   @time img = rescopeimage(img, offset, BB)
#   index = finished(index)
#   println("Writing ", get_name(index))
#   f = h5open(get_path(index), "w")
#   @time f["img", "chunk", (1000,1000)] = img
#   f["offset"] = offset
#   f["bb"] = [BB.i, BB.j, BB.w, BB.h]
#   close(f)
# end
