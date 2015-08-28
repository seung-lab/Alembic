using Base.Test
import Base: +, ==

type BoundingBox{T}
  i::T
  j::T
  h::T    # height
  w::T    # width
end

BoundingBox() = BoundingBox(0,0,0,0)
BoundingBox(a, b, c, d) = BoundingBox(promote(a,b,c,d)...)

"""
Add bounding boxes by finding BoundingBox that encompasses both
"""
function +(bbA::BoundingBox, bbB::BoundingBox)
  i = min(bbA.i, bbB.i)
  j = min(bbA.j, bbB.j)
  h = max(bbA.h+bbA.i, bbB.h+bbB.i) - i
  w = max(bbA.w+bbA.j, bbB.w+bbB.j) - j
  return BoundingBox(i,j,h,w)
end

"""
Test if BoundingBoxes have the same origin and dimensions
"""
function ==(bbA::BoundingBox, bbB::BoundingBox)
  return bbA.i == bbB.i && bbA.j == bbB.j && bbA.w == bbB.w && bbA.h == bbB.h
end

"""
Convert bounding box object to a polygon point list (counter-clockwise)
"""
function bb2pts(r)
  return [r.i r.j;
          r.i+r.h r.j;
          r.i+r.h r.j+r.w;
          r.i r.j+r.w;
          r.i r.j];
end

"""
Return xmin, ymin, xmax, ymax of a bounding box (opposing corner definition)
"""
function minsandmax(bb)
  return bb.i, bb.j, bb.i+bb.h, bb.j+bb.w
end

"""
`FIND_MESH_BB` - Find bounding box around mesh

    BoundingBox(xlow, ylow, height, width) = find_mesh_bb(nodes)

* `nodes`: 2xN matrix of mesh nodes
* `BoundingBox`: smallest integer-valued rectangle containing all mesh nodes

""" 
function find_mesh_bb(nodes)
    xlow = floor(Int64,minimum(nodes[:,1]))
    ylow = floor(Int64,minimum(nodes[:,2]))
    xhigh = ceil(Int64,maximum(nodes[:,1]))
    yhigh = ceil(Int64,maximum(nodes[:,2]))
    return BoundingBox(xlow, ylow, xhigh-xlow, yhigh-ylow)
end

"""
Snap bounding box to integer values (smallest rectangle containing original)
Returns BoundingBox of integers
"""
function snap_bb(bb)
  r = bb2pts(bb)
  i = floor(Int,bb.i)
  j = floor(Int,bb.j)
  h = ceil(Int,r[2,1]) - i
  w = ceil(Int,r[3,2]) - j
  return BoundingBox{Int}(i, j, h, w)
end

"""
Apply affine transform to points in a bounding box & find new bounding box
"""
function tform_bb(bb, tform)
  tform_pts = [bb2pts(bb) ones(size(bb2pts(bb),1),1)] * tform
  i = minimum(tform_pts[:,1])
  j = minimum(tform_pts[:,2])
  h = maximum(tform_pts[:,1]) - i
  w = maximum(tform_pts[:,2]) - j
  return BoundingBox(i, j, h, w)
end

"""
Convert Tuple for image size into a BoundingBox at (0,0)
"""
function sz2bb(sz)
  return BoundingBox(0, 0, sz[1]-1, sz[2]-1);  # should origin be at 1,1, 0,0, or 0.5,0.5?
end
