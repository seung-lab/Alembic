using Base.Test

# test_padimage()
o = ones(5,2)
po = padimage(o, 0, 0, 0, 0)
@test size(po, 1) == 5
@test size(po, 2) == 2

o = ones(5,2)
po = padimage(o, 1, 2, 3, 4)
@test size(po, 1) == 11
@test size(po, 2) == 6
@test po[1,1] == 0

o = ones(5,2)
po = padimage(o, 0, 1, 1, 0)
@test size(po, 1) == 6
@test size(po, 2) == 3
@test po[6,3] == 0    

o = convert(Array{Int64, 2}, ones(5,2))
po = padimage(o, 1, 2, 3, 4)
@test size(po, 1) == 11
@test size(po, 2) == 6

# test_rescopeimage()
img = ones(1,1)
offset = [0,0]
bb = BoundingBox(0,0,0,0)
@test rescopeimage(img, offset, bb) == img

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(1,1,0,0)
@test rescopeimage(img, offset, bb) == img

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(0,0,1,1)
@test rescopeimage(img, offset, bb) == [0 0; 0 1]

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(-1,-1,2,2)
@test rescopeimage(img, offset, bb) == [0 0 0; 0 0 0; 0 0 1]

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(0,0,2,2)
@test rescopeimage(img, offset, bb) == [0 0 0; 0 1 0; 0 0 0]

img = ones(1,1)
offset = [0,0]
bb = BoundingBox(-1,-1,1,1)
@test rescopeimage(img, offset, bb) == [0 0; 0 1]

img = rand(5,5)
offset = [0,0]
bb = BoundingBox(0,0,4,4)
@test rescopeimage(img, offset, bb) == img

bb = BoundingBox(0,0,1,1)
@test rescopeimage(img, offset, bb) == img[1:2, 1:2]

bb = BoundingBox(0,0,4,5)
@test rescopeimage(img, offset, bb)[1:5, 1:5] == img
@test rescopeimage(img, offset, bb)[1:5, 6:6] == zeros(5, 1)

bb = BoundingBox(0,-1,4,4)
@test rescopeimage(img, offset, bb)[1:5, 2:5] == img[1:5, 1:4]
@test rescopeimage(img, offset, bb)[1:5, 1:1] == zeros(5, 1)

bb = BoundingBox(0,-1,4,6)
@test rescopeimage(img, offset, bb)[1:5, 2:6] == img
@test rescopeimage(img, offset, bb)[1:5, 1:1] == zeros(5, 1)
@test rescopeimage(img, offset, bb)[1:5, 7:7] == zeros(5, 1)

bb = BoundingBox(-1,0,6,4)
@test rescopeimage(img, offset, bb)[2:6, 1:5] == img
@test rescopeimage(img, offset, bb)[1:1, 1:5] == zeros(1, 5)
@test rescopeimage(img, offset, bb)[7:7, 1:5] == zeros(1, 5)

bb = BoundingBox(-3,-3,2,2)
@test rescopeimage(img, offset, bb) == zeros(3,3)

img = rand(5,5)
offset = [-1,-1]
bb = BoundingBox(0,0,4,4)
@test rescopeimage(img, offset, bb)[1:4, 1:4] == img[2:5, 2:5]

# test_imfuse()
A = rand(5,5)
B = rand(5,5)
BB_A = [0, 0]
BB_B = [0, 0]
O, BB_O = imfuse(A, BB_A, B, BB_B)
@test O.channels[1] == A
@test O.channels[2] == B
@test BB_O == BB_A
@test BB_O == BB_B

A = rand(5,5)
B = rand(5,5)
BB_A = [0, 0]
BB_B = [2, 3]
O, BB_O = imfuse(A, BB_A, B, BB_B)
@test O.channels[1][1:5, 1:5] == A
@test O.channels[2][3:7, 4:8] == B
@test BB_O == BB_A

# try
# 	a = [1:10][0]
# catch e
# 	log_render_error(pwd(), (1,1,-2,-2), e)
# end
# rm(joinpath(pwd(), "render_error_log.txt"))