"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp_mesh(image, mesh::Mesh)
  scale = get_scale(:render)/get_scale(:match)
  src_nodes = get_nodes(mesh; globalized=true, use_post=false, scale=scale)
  dst_nodes = get_nodes(mesh; globalized=true, use_post=true, scale=scale)
  offset = get_offset(mesh)*scale;
  edges = get_edges(mesh)
  return @time meshwarp(image, src_nodes, dst_nodes, 
                    incidence_to_triangles(edges), offset), get_index(mesh)
end

"""
`UNSAFE_MASK_IMAGE!` - Mask src image when mask != mask_id & write to dst image
  
    unsafe_mask_image!(src, mask, mask_id, dst)

* `src`: Array
* `mask`: Array, identical size and shape to src (unchecked)
* `mask_id`: Number (should be of at least one element in mask)
* `dst`: Array, identical size and shape to src (unchecked)
* `val`: Value to be fill in the dst array when src does not equal mask_id

This method is unsafe, because there is no check that src, mask, & dst are all
the same size and type. This method will catastrophically fail if they are not.
"""
function unsafe_mask_image!(src, mask, mask_id, dst, val=0)
  @simd for i in 1:length(src)
    @inbounds dst[i] = src[i]
  end
  for i in 1:length(src)
    @inbounds if mask[i] != mask_id
      dst[i] = val
    end
  end
end

"""
`MERGE_IMAGES` - Place images in global reference image

    merged_image, 2D_slice = merge_images(images, offsets)

* `images`: 1D array, images (2D arrays)
* `offsets`: 1D array, 2-element array positions of corresponding image 
  in 'images` in global space

""" 
function merge_images{T}(images::Array{Array{T,2},1}, offsets)
    # T = typeof(images[1][1])
    bbs = []
    for (image, offset) in zip(images, offsets)
        push!(bbs, ImageRegistration.BoundingBox(offset..., size(image)...))
    end
    global_ref = sum(bbs)
    merged_image = zeros(T, global_ref.h, global_ref.w)
    no_images = length(images)
    for (idx, (image, bb)) in enumerate(zip(images, bbs))
        println("Merging image # ", idx , " / ", no_images)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        merged_image[i:i+h, j:j+w] = max.(merged_image[i:i+h, j:j+w], image)
        images[idx] = typeof(image)(0,0)
    end
    return merged_image, bb_to_slice(global_ref)
end

"""
Rescope image from one bounding box to another
"""
function rescope{T}(image::Array{T}, src_slice, dst_slice)
    src_bb = ImageRegistration.slice_to_bb(src_slice)
    dst_bb = ImageRegistration.slice_to_bb(dst_slice)
    src_offset = ImageRegistration.get_offset(src_bb)
    dst_offset = ImageRegistration.get_offset(dst_bb)
    dst = zeros(T, dst_bb.h, dst_bb.w)
    if intersects(src_bb, dst_bb)
        src_roi = translate_bb(dst_bb-src_bb, -src_offset+[1,1])
        dst_roi = translate_bb(dst_bb-src_bb, -dst_offset+[1,1])
        dst[bb_to_slice(dst_roi)...] = image[bb_to_slice(src_roi)...]
    end
    return dst
end

function limit_range(pts::Vector, offset::Int64, max_val::Int64)
  pts = pts-offset
  pts[pts .<= 0] = 1
  pts[pts .> max_val] = max_val
  return pts
end

function intersect_poly_bbox(pts::Matrix, offset::Vector, image_size::Vector)
  for j in 1:size(pts,2)
    pts[:,j] = limit_range(pts[:,j], offset[j], image_size[j])
  end
  return pts
end

"""
`RENDER` - Scale model to resolution & transform images as piecewise affine

    render!(ms)
    render!(ms, z_range)

* `ms`: MeshSet to render out
* `z_range`: list of z indices to be rendered (must have meshes in ms)

Render will compile all the subsections into one complete section & write it to
the dst_image directory specified in params at the render mip level.
"""
function render(ms::MeshSet, z_range=unique(collect_z(ms)))
  for z in z_range
    meshes = get_subsections(ms, z)
    subsection_images = Array{Array{UInt8,2},1}()
    subsection_offsets = []
    src_image = get_image(z, "src_image", mip=get_mip(:render))
    src_offset = get_offset("src_image", mip=get_mip(:render))
    src_size = get_image_size("src_image", mip=get_mip(:render)) 
    if use_roi_mask()
      src_roi = get_image(z, "roi", mip=get_mip(:render), input_mip=get_mip(:roi));
      mask_value = PARAMS[:roi][:mask_value]
      unsafe_mask_image!(src_image, src_roi, mask_value, src_image)
    end
    if use_defect_mask()
      src_image_sub = deepcopy(src_image)
      src_mask = get_image(z, "mask", mip=get_mip(:render), input_mip=get_mip(:match))
      for mesh in meshes
        index = get_index(mesh)
        println("Applying defect mask to $index")
        mask_id = get_subsection(index)
        unsafe_mask_image!(src_image, src_mask, mask_id, src_image_sub)
        println("Defect mask applied")
        println("Warping ", index)
        @time (dst_image, offset), _ = meshwarp_mesh(src_image_sub, mesh)
        push!(subsection_images, dst_image)
        push!(subsection_offsets, offset)
      end
      image, slice = merge_images(subsection_images, subsection_offsets)
      slice = tuple(slice..., z:z)
      println("Saving")
      save_image("dst_image", image, slice, mip=get_mip(:render))
    else
      mesh = get_mesh(ms, z)
      println("Warping ", z)
      @time (dst_image, offset), _ = meshwarp_mesh(src_image, mesh)
      @time save_image("dst_image", dst_image, offset, z, mip=get_mip(:render))
    end
  end
end
