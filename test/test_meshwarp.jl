using Base.Test

function warp_pts(affine, pts)
    pts = hcat(pts, ones(size(pts,1)))
    tpts = pts * affine
    return tpts[:,1:2]
end

# test_meshwarp()
img = reshape(float(collect(1:121).%2), 11, 11) # 11x11 checkerboard
triangles = [1 2 3]

tform = [1 0 0;
        0 1 0;
        0 0 1]
src = [0.0 0.0;
        0.0 10.0;
        10.0 0.0]   
dst = warp_pts(tform, src)
offset = [0,0]
img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
@test warped_offset == [0,0] 

tform = [1 0 0;
        0 1 0;
        0 0 1]
src = [0.0 0.0;
        0.0 10.0;
        10.0 0.0]   
dst = warp_pts(tform, src)
offset = [5,10]
img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
@test_approx_eq_eps img_warped zeros(11,11) 1e-10
@test warped_offset == [0,0] 

tform = [1 0 0;
        0 1 0;
        0 0 1]
src = [0.0 0.0;
        0.0 20.0;
        20.0 0.0]   
dst = warp_pts(tform, src)
offset = [5,10]
img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
@test img_warped[6,11] == 1.0 
@test warped_offset == [0,0] 

src = [2.0 2.0;
        10.0 2.0;
        6.0 10.0]
dst = [4.0 2.0;
        8.0 2.0;
        6.0 10.0]
offset = [0, 0]
img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
@test warped_offset == [4,2]

# test_poly2mask()
vertices=[1.5 2.5; 2.5 1.5; 3.5 4.5]
img, off = poly2mask(vertices[:,1],vertices[:,2])
""" off should be [1,1] and img should be:
false  false  false  false  false
false   true   true  false  false
false  false   true   true  false
false  false  false  false  false
"""
vertices=[1.5 2.4; 2.5 1.5; 3.5 4.5]
img, off = poly2mask(vertices[:,1],vertices[:,2])
""" off should be [1,1] and img should be:
false  false  false  false  false
false   true  false  false  false
false  false   true  false  false
false  false  false  false  false,
"""
vertices=[1 2; 2.5 1.5; 3.5 4.5]
img, off = poly2mask(vertices[:,1],vertices[:,2])
""" corner case: off should be [1,1] and img should be:
false   true  false  false  false
false   true   true  false  false
false  false   true   true  false
false  false  false  false  false
"""
# p=plot(x=vertices[:,1],y=vertices[:,2])
# display(p)
# img,off