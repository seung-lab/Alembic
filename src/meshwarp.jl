# Adapted from PiecewiseAffineTransforms.jl
# https://github.com/dfdx/PiecewiseAffineTransforms.jl

"""
`MESHWARP` - Apply piecewise affine transform to image using bilinear interpolation

    warped_img, [bb.i, bb.j] = meshwarp(img, src, dst, trigs, offset, interp)

* `img`: 2D array, image (todo: extend to Image type)
* `src`: Nx2 array, positions of mesh nodes in `img` (global space)
* `dst`: Nx2 array, positions of mesh nodes in `warped_img` (global space)
* `trigs`: Nx3 array, each row contains indices of nodes of a triangle
* `offset`: 2-element array, position of img[1,1] in global space (optional - default is [0,0])
* `interp`: bool determining whether to use bilinear interpolation or not
    (optional - default is true)
* `warped_img`: with pixel values the same type as original image (for Int type, pixel values are rounded)
* `warped_offset`: 2-element array, position of warped_img[1,1] in global space 

The bounding box of the warped image is defined as the smallest
integer-valued rectangle that contains the `dst` mesh.

This means that `warped_offset` is constrained to be integer-valued,
though `offset` is allowed to have floating point values.  The integer
constraint removes the need for further interpolation in any
subsequent fusing of multiple warped tiles.

See definitions in IMWARP documentation for further help.

""" 
function meshwarp{T}(img::Array{T},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                    trigs::Matrix{Int64}, offset=[0,0], interp=true)
    bb = snap_bb(find_mesh_bb(dst))
    warped_img = similar(img, bb.h+1, bb.w+1)
    warped_offset = [bb.i, bb.j]

    for t=1:size(trigs, 1)    
        tr = squeeze(trigs[t, :], 1)
        # coordinates of the source triangle
        X = src[tr, 1]
        Y = src[tr, 2]
        # coordinates of the destination triangle
        U = dst[tr, 1]
        V = dst[tr, 2]
        # Create matrix to transform from destination triangle to source image
        src_tri = [X Y ones(3,1)]
        dst_tri = [U V ones(3,1)]
        M = dst_tri \ src_tri # dst_tri * M = src_tri
        # Create list of coordinates in warped image that represent this
        # triangle (coordinates in the global space).
        us, vs = poly2source(U, V)
        
        # For every pixel in target triangle we find corresponding pixel in 
        # source and copy its value
        if interp
            for n=1:length(vs)
                # Convert warped coordinate to pixel space
                i, j = us[n]-warped_offset[1]+1, vs[n]-warped_offset[2]+1
                # Use warped coordinate in global space for transform
                u, v = us[n], vs[n]
                # x, y = M * [u, v, 1]
                x, y = M[1,1]*u + M[2,1]*v + M[3,1], M[1,2]*u + M[2,2]*v + M[3,2]
                # Convert original image coordinate to pixel space
                x, y = x-offset[1]+1, y-offset[2]+1
                fx, fy = floor(Int64, x), floor(Int64, y)
                wx, wy = x-fx, y-fy
                if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
                    # Expansion of p = [1-wy wy] * img[fy:fy+1, fx:fx+1] * [1-wx; wx]
                    p = round(((1-wy)*img[fx,fy] + wy*img[fx,fy+1]) * (1-wx) + ((1-wy)*img[fx+1,fy] + wy*img[fx+1,fy+1]) * wx)
                    warped_img[i, j] = p
                end
            end
        else
            for n=1:length(vs)
                i, j = us[n]-warped_offset[1]+1, vs[n]-warped_offset[2]+1
                u, v = us[n], vs[n]
                # x, y = M * [u, v, 1]
                x = M[1,1]*u + M[2,1]*v + M[3,1]
                y = M[1,2]*u + M[2,2]*v + M[3,2]
                x, y = round(Int64, x-offset[1]+1), round(Int64, y-offset[2]+1)
                if 1 <= x && x <= size(img, 1) && 1 <= y && y <= size(img, 2)
                    warped_img[i, j] = img[x, y]
                end
            end
        end
        
    end
    return warped_img, [bb.i, bb.j]
end

"""
`POLY2SOURCE` - Run fillpoly and findn over the basic bounding box of a given triangle
""" 
function poly2source(pts_i, pts_j)
    # Find bb of vertices (vertices in global space)
    top, bottom = floor(Int64,minimum(pts_i)), ceil(Int64,maximum(pts_i))
    left, right = floor(Int64,minimum(pts_j)), ceil(Int64,maximum(pts_j))
    # Create image based on number of pixels in bb that will identify triangle
    mask = zeros(Bool, bottom-top+1, right-left+1)
    # Convert vertices into pixel space and fill the mask to identify triangle
    fillpoly!(mask, pts_j-left+1, pts_i-top+1, true)
    # Create list of pixel coordinates that are contained by the triangle
    us, vs = findn(mask)
    # Convert that list of pixel coordinates back into global space
    us += top-1
    vs += left-1
    return us, vs
end

"""
`FILLPOLY!` - Fill pixels contained inside a polygon

    M = fillpoly!(M, px, py, value)

* `M`: 2D array, image
* `px`: 1D array, x-components of polygon vertices
* `py`: 1D array, y-components of polygon vertices
* `value`: fill value

Features:

* This seems intended to work for nonconvex polygons.
* The vertices should be listed sequentially in px and py. Clockwise/counterclockwise ordering doesn't seem to matter.
* It seems that the polygon is closed automatically if it isn't already, i.e., the last vertex in px, py is not equal to the first.  Does the code work if the input polygon is already closed?
* The code treats the "edge case," where an edge of the polygon lies exactly on a grid line.

Bugs:

* The original code (https://github.com/dfdx/PiecewiseAffineTransforms.jl/blob/master/src/polyline.jl) used implicit conversion to Int64 (presumably rounding).  Tommy/Shang replaced this by ceil.  This might produce inconsistent results, with the same pixel belonging to more than one triangle.
* The "corner case" where a grid line intersects a single vertex of the polygon does not appear to be properly treated.  The corner case is nongeneric if px and py are floats.  But the corner case could be common if px and py are ints, which seems encouraged by the parametric typing.
""" 
function fillpoly!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T)
    @assert length(py) == length(px)    
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    # Scan poly from left to right
    for x=left:right     # loop over grid lines
        ys = Set{Int64}()
        m = length(px)
        for n=1:length(px)  # loop over edges (m,1),(1,2),(2,3),...,(m-1,m)
            # grid line intersects edge in one of two ways
            if (px[n] <= x && x <= px[m]) || (px[m] <= x && x <= px[n])
                if px[n] == px[m]  # intersection is entire edge
                    push!(ys, ceil(Int64, py[n]))
                    push!(ys, ceil(Int64, py[m]))
                else # intersection is point
                      y = py[n] + (x-px[n]) * (py[m]-py[n])/(px[m]-px[n])
                    push!(ys, ceil(Int64, y))
                end            
            end
            m = n
        end
        # generically, two intersections for a convex polygon
        # generically, even number of intersections for a nonconvex polygon
        ys = sort([y for y in ys])  # sort the intersection points
        # if odd number of intersection points, add duplicate end point
        if length(ys) % 2 == 1
            push!(ys, ys[end])
        end
        # Place value in matrix at all the y's between min and max for given x
        for n=1:2:length(ys)           
            M[ys[n]:ys[n+1], x] = value
        end
    end
    return M
end

function demo_meshwarp()
# Demo the updated meshwarp function that runs faster than original package
    img = get_ufixed8_image(joinpath(BUCKET, "test_images", "turtle.jpg"))
    src_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    320.0 290.0]'
    dst_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    400.0 460.0]'
    incidence = [1 1 1 0 0 0 0 0;
                -1 0 0 1 1 0 0 0;
                0 0 0 -1 0 1 1 0;
                0 -1 0 0 0 0 -1 1;
                0 0 -1 0 -1 -1 0 -1]
    triangles = [1 2 5;
                1 4 5;
                2 3 5;
                3 4 5];
    node_dict = incidence2dict(incidence)
    draw_mesh(img, src_nodes, node_dict)
    println(size(img))

    warp = meshwarp(img, src_nodes, dst_nodes, triangles)
    draw_mesh(warp, dst_nodes, node_dict)
    println(size(warp))
end
