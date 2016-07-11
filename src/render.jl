"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp_mesh(mesh::Mesh)
  img = get_image(mesh)
  src_nodes = hcat(get_nodes(mesh; globalized = true, use_post = false)...)'
  dst_nodes = hcat(get_nodes(mesh; globalized = true, use_post = true)...)'
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
    println(index);
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
    warps = map(meshwarp_mesh, meshset.meshes);
      println(warps[1][1][1][1:10])
    imgs = [x[1][1] for x in warps];
    offsets = [x[1][2] for x in warps];
    indices = [x[2] for x in warps];
    # review images
    if render_review
      write_seams(meshset, imgs, offsets, indices, flagged_only)
    end
    if render_full
      # println(typeof(imgs));
      # println(typeof(imgs[1]));
      # println(size(imgs[1]));
      # println(imgs[1][1:10]);
      # println("test1")
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
    if render_full && montaged(index) == montaged(startindex)
      render_prealigned(src_index, dst_index, src_img, [], eye(3), 
                      eye(3); render_full=true, render_review=false)
    else
      if dst_img == nothing
        println("Loading dst_image for rendering")
        dst_img = get_image(dst_index)
      end
      render_prealigned(src_index, dst_index, src_img, dst_img, cumulative_tform, 
                      tform; render_full=render_full, render_review=render_review)

      println(index, " ", firstindex, " ", align && montaged(index) > montaged(firstindex))
      if align && montaged(index) > montaged(firstindex)
        println("Aligning meshes between $dst_index, $src_index")
        reload_registry(prealigned(src_index))
        ms = MeshSet(prealigned(dst_index), prealigned(src_index); solve=false, fix_first=(dst_index==startindex))
        render_aligned_review(ms)
      end
    end
    println("Swapping src_image to dst_image")
    dst_img = copy(src_img)
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
function render_aligned_review(firstindex::Index, lastindex::Index, start=1, finish=0; scale=0.05)
  firstindex, lastindex = prealigned(firstindex), prealigned(lastindex)
  meshset = load(firstindex, lastindex)
  render_aligned_review(meshset, start, finish, scale=scale)
end

function render_aligned_review(meshset, start=1, finish=length(meshset.matches); images=Dict(), scale=0.05)
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

@fastmath @inbounds function render_aligned(meshset::MeshSet, start=1, finish=length(meshset.meshes))
  sort!(meshset.meshes; by=get_index)
  subsection_imgs = []
  subsection_offsets = []
  for (k, mesh) in enumerate(meshset.meshes)
    if start <= k <= finish
      index = get_index(mesh)
      println("Warping ", index)
      @time (img, offset), _ = meshwarp_mesh(mesh)
      if is_subsection(index)
        println("$index is a subsection")
        push!(subsection_imgs, img)
        push!(subsection_offsets, offset)
      end
      # determine if subsections should be merged 
      is_last_subsection = true
      if k != length(meshset.meshes)
        next_index = get_index(meshset.meshes[k+1])
        if prealigned(index) == prealigned(next_index)
          is_last_subsection = false
          println("Wait to merge...")
        end
      end
      if length(subsection_imgs) > 1 && is_last_subsection
        println("Merge subsections")
        img, offset = merge_images(subsection_imgs, subsection_offsets)
        subsection_imgs = []
        subsection_offsets = []
      end
      # render if subsections have been merged or is not a split section
      if is_last_subsection || !is_subsection(index)
        println("Writing ", get_name(aligned(index)))
        f = h5open(get_path(aligned(index)), "w")
        chunksize = min(1000, min(size(img)...))
        @time f["img", "chunk", (chunksize, chunksize)] = img
        f["dtype"] = string(typeof(img[1]))
        f["offset"] = offset
        f["size"] = [size(img)...]
        close(f)
        # Log image offsets
        update_offset(aligned(index), offset, size(img))
      end
    end
  end
end

function split_prealigned(index::Index)
  mask_path = get_mask_path(index)
  if isfile(mask_path)
    println("Splitting $index with mask")
    mask = load_mask(mask_path)
    img = get_image(index)
    subimgs = segment_by_mask(img, mask)
    offset = get_offset(index)
    for (i, subimg) in subimgs
      n = length(subimgs)
      subindex = subsection(index, i)
      path = get_path(subindex)
      println("Saving subsection $subindex: $i / $n")
      f = h5open(path, "w")
      chunksize = min(1000, min(size(subimg)...))
      @time f["img", "chunk", (chunksize,chunksize)] = subimg
      f["dtype"] = string(typeof(subimg[1]))
      f["offset"] = offset
      f["size"] = [size(subimg)...]
      close(f)
      # Log image offsets
      update_offset(subindex, offset, size(subimg))
    end
  end
end

function split_prealigned(firstindex::Index, lastindex::Index)
  for index in get_index_range(firstindex, lastindex)
    mask_path = get_mask_path(index)
    if isfile(mask_path)
      split_prealigned(index)
    end
  end
end
