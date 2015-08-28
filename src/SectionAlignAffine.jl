include("convolve.jl")


using Images
include("imwarp.jl")
include("visualize.jl")
#include("render.jl")	# cyclic inclusion


function test()
	getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)));
	sec1 = getimage("./sections/S2-W001_sec24_0.175.tif")
	#sec2 = getimage("./sections/S2-W001_sec24_0.175.tif")
	sec2 = getimage("./sections/S2-W001_sec25_0.175.tif")
	println(size(sec1))
	println(size(sec2))
	trans, points1, points2, res1, res2 = AffineAlignSections(sec1, sec2, 0.175, return_points=true)

	println(trans)
	#out_img, offset = imwarp(sec2, trans)
	#imwrite(sec1, joinpath(".","test_outputs", string("sec1", ".tif")));
	#imwrite(sec2, joinpath(".","test_outputs", string("sec2_", offset[1], "_", offset[2], ".tif")));
	#imwrite(out_img, joinpath(".","test_outputs", string("warped_", offset[1], "_", offset[2], ".tif")));
	p22 = points2 + res2;
	p11 = points1 + res1;
	points1 = [points1[2,:]; points1[1,:]]
	points2 = [points2[2,:]; points2[1,:]]
	draw_vectors(sec2, vcat(points2[1:2,:], p22[2:-1:1,:]))
	#ccp.write_image_from_points(points1[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
	draw_points(sec1, points1)
	draw_points(sec2, points2)
	#draw_points(out_img, points1)
end


function test2()
	getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)));
	sec1 = getimage("./output_images/(1,1)_montage.tif")
	sec2 = getimage("./output_images/(1,2)_montage.tif")
	println(size(sec1))
	println(size(sec2))
	trans, points1, points2, res1, res2 = AffineAlignSections(sec1, sec2, 1, 0.3; return_points=true)

	println(trans)
	downsample = 4
	sec1 = sec1[1:downsample:end, 1:downsample:end];
	sec2 = sec2[1:downsample:end, 1:downsample:end];

	p11 = points1 + res1;
	p22 = points2 + res2;

	p11 = p11[2:-1:1,:]
	p22 = p22[2:-1:1,:]
	points1 = points1[2:-1:1,:]
	points2 = points2[2:-1:1,:]

	p11 = p11/downsample;
	p22 = p22/downsample;
	points1 = points1/downsample;
	points2 = points2/downsample;
	
	trans = AdjustAffineForScaling(trans.', downsample).'
	out_img, offset = imwarp(sec2, inv(trans))
	println(offset)
	fused, fused_offset = imfuse(sec1, [0,0], out_img, offset)
	println(fused_offset)
	fused_offset = fused_offset[2:-1:1]

	half_block_size = 150;
	scalebar = [1; 1; 2*half_block_size/downsample; 2*half_block_size/downsample]
	draw_vectors(make_isotropic(sec1), hcat(vcat(points1, p11), scalebar))
	draw_vectors(make_isotropic(sec2), hcat(vcat(points2, p22), scalebar+100))
	draw_vectors(make_isotropic(fused), hcat(vcat(points1.-fused_offset, p11.-fused_offset), scalebar)) # todo: check offset
    #ccp.write_image_from_points(points1[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
end


"""
Returns the affine transformation A, such that p_in_B = p_in_A * A, where p is point coordinates in row vector.
"""
function AffineAlignSections(img1::Array{}, img2::Array{}, downsample_ratio = 1, accept_xcorr = 0.3; return_points=false)
	# points here are in column vector convention
	points, half_block_size, search_radius = GenerateMatchPoints(img1, img2)
	points1, points2 = GetBlockMatches(img1, img2, points, half_block_size, search_radius, accept_xcorr)
	A = FindAffine(points1, points2)
	residualIn2 = A*points1 - points2;
	rmsIn2 = mean( sum(residualIn2.^2, 1) )^0.5
	points2in1 = inv(A)*points2;
	residualIn1 = points2in1 - points1;
	rmsIn1 = mean( sum(residualIn1.^2, 1) )^0.5

	tolerance_ratio = 0.001
	rmsThres = tolerance_ratio * mean([size(img1)..., size(img2)...])
	rmsTotal = mean([rmsIn1, rmsIn2].^2)^0.5
	if  rmsTotal > rmsThres
		println("WARNING [AffineAlignSections]: high residual. RMS error: ", rmsTotal)
	end

	A = AdjustAffineForScaling(A, downsample_ratio)

	# convert column vector convention to row vector convention
	A = A.'

	if return_points
		return A, points1, points2, residualIn1, residualIn2, rmsIn1, rmsIn2
	else
		return A
	end
end


function GenerateMatchPoints(img1::Array{}, img2::Array{})
	border_ratio = 0.1;
	radius_ratio = 0.05;
	grid_size = 3;
	half_block_size = 150;
	overlap = [min(size(img1), size(img2))...]
	#overlap = size(img1)

	border = round(Int, border_ratio * minimum(overlap));
	search_radius = round(Int, minimum(overlap) * radius_ratio);
	search_radius = search_radius > 200 ? search_radius : 200;
	println("search_radius: ", search_radius)

	border = border > search_radius + half_block_size ? border : search_radius + half_block_size;
	println("border excluded: ", border)

	# four corners for now
	#points = [[x; y] for x = (1+border,size(img1,1)-border), y = (1+border,size(img1,2)-border)]
	#points = [[x; y] for x = (1+border,overlap[1]-border), y = (1+border,overlap[2]-border)]

	step = floor(Int, (overlap - 2*border - 1) ./ (grid_size-1));
	points = [[x; y] for x = 1+border:step[1]:overlap[1]-border, y = 1+border:step[2]:overlap[2]-border]
	points = points[:];
	println("n points: ", size(points))

	return points, half_block_size, search_radius
end

function GetBlockMatches(img1::Array{}, img2::Array{}, points, half_block_size, search_radius, accept_xcorr = 0.3)
	points1 = Array{Int,1}[]
	points2 = Array{Int,1}[]
	for p = points
		offset,r = BlockMatchAtPoint(img1, p, img2, p, half_block_size, search_radius)
		println(p, offset, r)
		if r >= accept_xcorr
			push!(points1, [p; 1])
			push!(points2, [p + offset; 1])
		end
	end
	points1 = hcat(points1[:]...)
	points2 = hcat(points2[:]...)
	return points1, points2
end


function AdjustAffineForScaling(A, scale)
	# Column vector convention. I.e. p_transformed = A * p.
	ss = 1. / scale;
	B = copy(A)
	# equivalent of inv(S) * A * S, where S = [1 0 0; 0 1 0; 0 0 1/scale]
	B[1:end-1, end] *= ss
	return B
end


function FindAffine(points1, points2)
# Inputs:
# 		2D arrays with each column being a point in homogeneous coords.
# Returns:
#    	Affine transformation A, that produces points2 = A * points1

	A = points2 / points1
end

function BlockMatchAtPoint(A, pointInA, B, pointInB, half_block_size, search_r)
# A is the template
# Returns:
#	offset: offset from pointInB to pointInA's actual match in B.
#	r_max:  correlation value at the match point
#	xc:		raw cross-correlation map

	b = half_block_size
	B_radius = half_block_size + search_r;

	A_lower = pointInA - b
	A_upper = pointInA + b
	B_lower = pointInB - B_radius
	B_upper = pointInB + B_radius

	xc = normxcorr2(sub(A, [l:u for (l,u) in zip(A_lower, A_upper)]...),
					sub(B, [l:u for (l,u) in zip(B_lower, B_upper)]...));

	r_max, ind = findmax(xc); 
	i_max, j_max = ind2sub(xc, ind);

	return [i_max - 1 - search_r; j_max - 1 - search_r], r_max, xc;
	
end
