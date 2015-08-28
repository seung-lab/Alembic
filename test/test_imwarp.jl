using Base.Test

# test imwarp()
img = reshape(float(collect(1:121).%2), 11, 11) # 11x11 checkerboard
warp = img

tform = [1 0 0;
        0 1 0;
        0 0 1]
offset = [0, 0]
img_warped, warped_offset = imwarp(img, tform, offset)
@test_approx_eq warp img_warped
@test warped_offset == [0, 0]

tform = [1 0 0;
        0 1 0;
        0 0 1]
offset = [10, 20]
img_warped, warped_offset = imwarp(img, tform, offset)
@test_approx_eq warp img_warped
@test warped_offset == [10, 20]

tform = [1 0 0;
        0 1 0;
        10 10 1]
offset = [0, 0]
img_warped, warped_offset = imwarp(img, tform, offset)
@test_approx_eq warp img_warped
@test warped_offset == [10, 10]

tform = [cos(0.5) -sin(0.5) 0;
        sin(0.5) cos(0.5) 0;
        0 0 1]
offset = [0.0, 0.0]
img_warped, warped_offset = imwarp(img, tform, offset)
@test warped_offset == [0, -5]