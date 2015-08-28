% speed tests
%% imwarp with simple affine
a = uint8(randi(255, 8000, 8000));
tform = affine2d([1 0 0; 0 1 0; 0.5 0.5 1]);
tic; [b, r] = imwarp(a, tform); toc;

%% imwarp with piecewise affine
load('/usr/people/tmacrina/seungmount/research/tommy/Julimaps/test_images/MATLAB_speed_test_vars.mat')
tile = imread('/usr/people/tmacrina/seungmount/research/tommy/Julimaps/EM_images/Tile_r4-c2_S2-W001_sec20.tif');
tic;
tformPiecewiseLinear = images.geotrans.PiecewiseLinearTransformation2D(moving', fixed');
d = imwarp(tile, tformPiecewiseLinear);
toc;