function get_bb(index::Index)
  return ImageRegistration.BoundingBox(get_offset(index)..., get_image_size(index)...)
end

"""
Is point contained within the bounding box (border included)?
"""
function point_is_contained(bb::ImageRegistration.BoundingBox, pt::Point)
  return (bb.i <= pt[1] <= (bb.i+bb.h)) && (bb.j <= pt[2] <= (bb.j+bb.w))
end

"""
At least one point of line is contained within bounding box (border included)?
"""
function line_is_contained(bb::ImageRegistration.BoundingBox, line)
  return point_is_contained(bb, line[1:2]) || point_is_contained(bb, line[3:4])
end

"""
Given list of bounding boxes, calculate pairs of indices with overlaps, symmetric
"""
function find_overlaps(boundingboxes)
  bbs = copy(boundingboxes)
  overlap_tuples = []
  i = length(bbs)
  while length(bbs) != 0
    bbi = pop!(bbs)
    for (j, bbj)  in enumerate(bbs)
      if intersects(bbi, bbj)
        push!(overlap_tuples, (i,j))
        push!(overlap_tuples, (j,i))
      end
    end
    i -= 1
  end
  return overlap_tuples
end

function find_boundingboxes(meshset)
  nodes = [get_globalized_nodes_h(mesh)[2]' for mesh in meshset.meshes]
  return map(ImageRegistration.find_mesh_bb, nodes)
end

"""
Crop image with offset to bounding box
"""
function imcrop(img, offset, bb)
  o = zeros(eltype(img), ceil(bb.h)+1, ceil(bb.w)+1)
  ibb = ImageRegistration.BoundingBox(offset..., size(img)...)
  d = bb - ibb
  o_start = abs(bb.i-d.i)+1:abs(bb.i-d.i) + d.h
  o_end = abs(bb.j-d.j)+1:abs(bb.j-d.j) + d.w
  im_start = abs(ibb.i-d.i)+1:abs(ibb.i-d.i) + d.h
  im_end = abs(ibb.j-d.j)+1:abs(ibb.j-d.j) + d.w
  o[o_start, o_end] = img[im_start, im_end]
  return o
end

"""
`WRITE_SEAMS` - Write out overlays of montaged seams
""" 
function write_seams(meshset, imgs, offsets, indices, flagged_only=true)
  bbs = []
  for (img, offset) in zip(imgs, offsets)
      push!(bbs, ImageRegistration.BoundingBox(offset..., size(img)...))
  end
  overlap_tuples = find_overlaps(bbs) # could include tag for asymmetric list
  total_seams = flagged_only ? count_flags(meshset) : length(overlap_tuples)
  for (k, (i,j)) in enumerate(overlap_tuples)
    src_index, dst_index = indices[i], indices[j]
    ind = find_match_index(meshset, src_index, dst_index)
    if ind > 0
      if !flagged_only || is_flagged(meshset.matches[ind])
        println("Writing match #", k, " of ", total_seams, " seams")
        path = get_path("review", (src_index, dst_index))
        img, fuse_offset = imfuse(imgs[i], offsets[i], imgs[j], offsets[j])
        bb = bbs[i] - bbs[j]
        img_cropped = imcrop(img, fuse_offset, bb)
        f = h5open(path, "w")
        chunksize = min(50, min(size(img_cropped)...))
        @time f["img", "chunk", (chunksize,chunksize)] = img_cropped
        f["offset"] = [bb.i, bb.j]
        f["scale"] = 1.0
        close(f)
      end
    end
  end
end

"""
Create CairoSurface of bounding boxes
"""
function draw_bbs(bbs, indices)
  padding = [100, 100]
  fontsize = 36
  index = indices[1]
  global_bb = snap_bb(sum(bbs))
  bbs = map(translate_bb, bbs, repeated(-ImageRegistration.get_offset(global_bb)+padding))
  sz = ImageRegistration.get_size(global_bb) + 2*padding
  drw = create_drawing(ones(UInt32, sz...))
  ctx = get_context(drw)
  rects = map(get_rect, bbs)
  colors = ([1,0,1], [0,1,1])
  txt = join(index[1:2], ",")
  draw_text(ctx, txt, [50,50], [-10,-10], fontsize, [1,1,1])
  for (k, (idx, rect)) in enumerate(zip(indices, rects))
    draw_rect(ctx, rect, colors[k%2+1])
    ctr = [rect[1]+rect[3]/2, rect[2]+rect[4]/2]
    txt = join(idx[3:4], ",")
    draw_text(ctx, txt, ctr, [0,-10], fontsize, colors[k%2+1])
  end
  return drw
end

function draw_polys(polys, indices, roi=nothing)
  padding = [100, 100]
  fontsize = 36
  index = indices[1]
  x = vcat([vertices[:,1] for vertices in polys]...)
  y = vcat([vertices[:,2] for vertices in polys]...)
  min_x = minimum(x)
  min_y = minimum(y)
  max_x = maximum(x)
  max_y = maximum(y)
  sz = [max_x-min_x+1, max_y-min_y+1] + 2*padding
  sz = round(Int64, sz)
  polys = [[vertices[:,1]+padding[1]-min_x vertices[:,2]+padding[2]-min_y] for vertices in polys]
  drw = create_drawing(ones(UInt32, sz...))
  ctx = get_context(drw)
  colors = ([1,0,1], [0,1,1])
  txt = join(index[1:2], ",")
  draw_text(ctx, txt, [50,50], [-10,-10], fontsize, [1,1,1])
  for (k, (idx, poly)) in enumerate(zip(indices, polys))
    draw_poly(ctx, poly, colors[k%2+1])
    min_poly_x = minimum(poly[:,1])
    max_poly_x = maximum(poly[:,1])
    min_poly_y = minimum(poly[:,2])
    max_poly_y = maximum(poly[:,2])
    ctr = [(max_poly_y-min_poly_y)/2 + min_poly_y, (max_poly_x-min_poly_x)/2 + min_poly_x]
    txt = join(idx[3:4], ",")
    draw_text(ctx, txt, ctr, [0,-10], fontsize, colors[k%2+1])
  end
  if roi != nothing
    roi = [roi[:,1]+padding[1]-min_x roi[:,2]+padding[2]-min_y] 
    draw_poly(ctx, roi, [1,0,0])
  end    
  return drw
end

function view_polys(polys, indices, roi=nothing)
  drw = draw_polys(polys, indices, roi)
  img = convert_drawing(get_drawing(drw))'
  return ImageView.view(img, pixelspacing=[1,1])
end

function view_bbs(bbs::Array{ImageRegistration.BoundingBox, 1}, indices)
  drw = draw_bbs(bbs, indices)
  img = convert_drawing(get_drawing(drw))'
  return ImageView.view(img, pixelspacing=[1,1])
end

"""
Create CairoSurface of tile outlines based on premontage registry
"""
function draw_premontage_review(index::Index; scale=0.05)
  indices = get_index_range(premontaged(index), premontaged(index))
  bbs = map(scale_bb, map(get_bb, indices), repeated(scale))
  return draw_bbs(bbs, indices)
end

function view_premontage_review(index::Index; scale=0.05)
  drw = draw_premontage_review(index, scale=scale)
  img = convert_drawing(get_drawing(drw))'
  return ImageView.view(img, pixelspacing=[1,1])
end

function save_premontage_review(index::Index; scale=0.05)
  drw = draw_premontage_review(index, scale=scale)
  idx = (index[1:2]..., 0, 0)
  fn = get_path("outline", premontaged(index))
  println("Saving premontage review: $fn")
  Cairo.write_to_png(drw, fn)
end

function view_prematch_review(index::Index)
  path = get_path("review", (index, get_preceding(index))); 
  img_orig = h5read(path, "img"); 
  ImageView.view(img_orig);
end

function split_premontage_review(index::Index)
  f = get_path("outline", premontaged(index))
  im = data(FileIO.load(fn))
  h = round(Int64, size(im)[2]/2)
  fn = split(f, ".png")[1]
  fn1 = joinpath("split", string(fn, "1", ".png"))
  fn2 = joinpath("split", string(fn, "2", ".png"))
  FileIO.save(fn1, Images.Image(im[:,1:h+200])')
  FileIO.save(fn2, Images.Image(im[:,h-200:end])')
end
