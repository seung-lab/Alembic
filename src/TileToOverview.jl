include("convolve.jl")

using Images
using ImageView
include("visualize.jl")
include("render.jl")

#	tile = getimage("./input_images/Tile_r1-c2_S2-W001_sec1.tif")
#	overview = getimage("./input_images/S2-W001_sec1_overview.tif")

"""
`TileToOverview`

Returns:

* offset of the tile in the overview, in the scale of the overview image

Args:

* tileFile: path to the tile image file.
* overviewFile: path to the overview image file.
* overviewScaling: scaling factor of the overview image compared to the full resolusion tile image.

"""
function TileToOverview(tileFile::String, overviewFile::String, overviewScaling::Real; diagnosis = false)

	getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)));
	tile = getimage(tileFile)
	overview = getimage(overviewFile)

	s = overviewScaling

	tilesize = size(tile)
	range = round(Int, 1:1/s:tilesize[1]), round(Int, 1:1/s:tilesize[2])
	resampled = tile[range...]

	offset, r_max, xc = BlockMatch(resampled, overview)

	if diagnosis
		println(offset, " ", r_max)
		xc -= minimum(xc)
		xc /= maximum(xc)
		view(make_isotropic(xc))
		fused, fused_offset = imfuse(overview, [0,0], resampled, offset)
		view(make_isotropic(fused))
	end

	return offset
end


function BlockMatch(downsampledTile, overview)
# Returns:
#	offset: offset of the tile from the overview.
#	r_max:  correlation value at the match point
#	xc:		raw cross-correlation map

	xc = normxcorr2(downsampledTile, overview);

	r_max, ind = findmax(xc);
	i_max, j_max = ind2sub(xc, ind);

	return [i_max - 1; j_max - 1], r_max, xc;
	
end
