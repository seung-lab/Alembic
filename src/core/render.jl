"""
Multiple dispatch for meshwarp on Mesh object
"""
function meshwarp_mesh(mesh::Mesh)
  img = get_image(mesh, mip=get_mip(:render))
  scale = get_scale(:match)/get_scale(:render)
  src_nodes = get_nodes(mesh; globalized=true, use_post=false, scale=scale)
  dst_nodes = get_nodes(mesh; globalized=true, use_post=true, scale=scale)
  offset = get_offset(mesh);
  #=print("incidence_to_dict: ")
  @time node_dict = incidence_to_dict(mesh.edges') #'
  print("dict_to_triangles: ")
  @time triangles = dict_to_triangles(node_dict)=#
  return @time meshwarp(img, src_nodes, dst_nodes, incidence_to_triangles(mesh.edges), offset), get_index(mesh)
end

"""
`MERGE_IMAGES` - Place images in global reference image

    merged_img, 2D_slice = merge_images(imgs, offsets)

* `imgs`: 1D array, images (2D arrays)
* `offsets`: 1D array, 2-element array positions of corresponding image 
  in 'imgs` in global space

""" 
function merge_images{T}(imgs::Array{Array{T,2},1}, offsets)
    # T = typeof(imgs[1][1])
    bbs = []
    for (img, offset) in zip(imgs, offsets)
        push!(bbs, ImageRegistration.BoundingBox(offset..., size(img)...))
    end
    global_ref = sum(bbs)
    merged_img = zeros(T, global_ref.h, global_ref.w)
    no_imgs = length(imgs)
    for (idx, (img, bb)) in enumerate(zip(imgs, bbs))
        println("Merging image # ", idx , " / ", no_imgs)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        merged_img[i:i+h, j:j+w] = max.(merged_img[i:i+h, j:j+w], img)
        imgs[idx] = typeof(img)(0,0)
    end
    return merged_img, bb_to_slice(global_ref)
end

"""
Rescope image from one bounding box to another
"""
function rescope{T}(img::Array{T}, src_slice, dst_slice)
    src_bb = ImageRegistration.slice_to_bb(src_slice)
    dst_bb = ImageRegistration.slice_to_bb(dst_slice)
    src_offset = ImageRegistration.get_offset(src_bb)
    dst_offset = ImageRegistration.get_offset(dst_bb)
    dst = zeros(T, dst_bb.h, dst_bb.w)
    if intersects(src_bb, dst_bb)
        src_roi = translate_bb(dst_bb-src_bb, -src_offset+[1,1])
        dst_roi = translate_bb(dst_bb-src_bb, -dst_offset+[1,1])
        dst[bb_to_slice(dst_roi)...] = img[bb_to_slice(src_roi)...]
    end
    return dst
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
@fastmath @inbounds function render!(ms::MeshSet, z_range=unique(collect_z(ms)))
  for z in z_range
    meshes = get_subsections(ms, z)
    subsection_imgs = Array{Array{UInt8,2},1}()
    subsection_offsets = []
    for mesh in meshes
      index = get_index(mesh)
      println("Warping ", index)
      @time (img, offset), _ = meshwarp_mesh(mesh)
      push!(subsection_imgs, img)
      push!(subsection_offsets, offset)
    end
    img, slice = merge_images(subsection_imgs, subsection_offsets)
    slice = tuple(slice..., z:z)
    save_image(z, "dst_image", img, slice)
  end
end
