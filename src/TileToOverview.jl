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
function tile_to_overview(tile, overview, overview_scale::Real; 
				diagnosis = false, overlay_array::Array = [], xcorr_overlay::Array = [])

	# Assuming strings are file paths
	if isa(tile, String)
		tile = get_image(tile)
	end
	if isa(overview, String)
		overview = get_image(overview)
	end
	return tile_to_overview(tile, overview, overview_scale; 
			diagnosis = diagnosis, overlay_array = overlay_array, xcorr_overlay = xcorr_overlay)
end

function tile_to_overview(tile_img::Array, overview_img::Array, overview_scale::Real; 
				diagnosis = false, overlay_array::Array = [], xcorr_overlay::Array = [])

	tile = tile_img
	overview = overview_img
	s = overview_scale

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

	# Add onto overlay array
	if length(overlay_array) > 0
		overlay_size = [size(overlay_array)...]
		resampled_tile_size = [size(resampled)...]

		# if size(overlay_array) is bigger than size(overview), center align them.
		offset_in_overlay = [offset...] + floor(Int, (overlay_size - [size(overview)...]) / 2)
		end_in_overlay = offset_in_overlay + resampled_tile_size

		# account for any potential cropping
		cropped_offset = max(zeros(Int, 2), offset_in_overlay)
		cropped_end = min(overlay_size, end_in_overlay)

		range_in_overlay = [a:b for (a,b) in zip(cropped_offset+1, cropped_end)]
		range_in_resampled = [a:b for (a,b) in 
				zip(cropped_offset-offset_in_overlay+1, resampled_tile_size-end_in_overlay+cropped_end)]
		overlay_array[range_in_overlay...] += resampled[range_in_resampled...]
	end

	# Add onto xcorr overlay, masked to tile size centered at peak correlation
	if length(xcorr_overlay) > 0 && size(xcorr_overlay)==size(overview)
		xc_size = [size(xc)...]
		half_tile = floor(Int, [size(resampled)...] / 2)
		mask_begin = [offset...] - half_tile
		mask_end = [offset...] + half_tile

		# account for any potential cropping
		mask_begin = max(zeros(Int, 2), mask_begin)
		mask_end = min(xc_size, mask_end)

		range_in_xc = [a:b for (a,b) in zip(mask_begin+1, mask_end)]
		range_in_overlay = range_in_xc + half_tile
		xcorr_overlay[range_in_overlay...] += xc[range_in_xc...]
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


function tiles_to_overview(tile_img_file_list::Vector{ByteString}, 
                            overview_img,
                            overview_scale::Real;
                            scale_offsets_to_tile_unit = true,
                            tile_img_dir = "",
                            save_fused_img_to = "",
                            save_xcorr_img_to = "",
                            show_review_imgs = false)

                            #params=PARAMS_PREMONTAGE
  if isa(overview_img, String)
  	overview = get_image(overview_img)
  elseif isa(overview_img, Array)
  	overview = overview_img
  else
  	error("Invalid argument overview_img")
  end
  overlay = zeros(Float64, size(overview))
  xcorr = zeros(Float64, size(overview))
  #offsets = Dict{Index, Array{Float64, 1}}()
  offsets = Vector{Tuple{String,Vector}}()  # Can be used to easily construct a Dict{String, Vector} if needed
  for tilefile in tile_img_file_list
  	print(tilefile, "  ")
  	tile = joinpath(tile_img_dir, tilefile)
  	gc()
  	offset = tile_to_overview(tile, overview, overview_scale; overlay_array = overlay, xcorr_overlay = xcorr)
  	if scale_offsets_to_tile_unit
  		push!(offsets, (tilefile, offset / overview_scale))
  	else
    	#offsets[parseName(tilefile)] = offset
    	#offsets[tilefile] = offset
      push!(offsets, (tilefile, offset))
  	end
  	println(offset)
    #break
  end

  overlay = overlay / 255;
  overview = overview / 255;

  if save_fused_img_to != ""
	  #fused, fused_offset = imfuse(overview, [0,0], overlay, [0,0]) # having trouble writing this to file
	  fused = cat(3, overview, min(overlay,1), zeros(size(overlay)))
	  fused = convert(Image, fused)
	  imwrite(fused, save_fused_img_to)
	  if show_review_imgs
		  imgc, imgs = view(fused)
		  ImageView.annotate!(imgc, imgs, ImageView.AnnotationText(0,0, splitdir(save_fused_img_to)[2],
		  		halign="left", valign="top", fontsize=round(minimum(size(overview))/50), color=RGB(0.5,0.5,0.5)))
	  end
  end
  if save_xcorr_img_to != ""
  	  println("xcorr  min: ", minimum(xcorr), "  max: ",maximum(xcorr))
  	  # clamp to [0,1]
  	  xcorr = max(xcorr, 0)
  	  xcorr = min(xcorr, 1)
	  # stretch high value to 1
	  xcorr /= maximum(xcorr)
	  imwrite(xcorr, save_xcorr_img_to)
	  if show_review_imgs
		  imgc, imgs = view(make_isotropic(xcorr))
		  ImageView.annotate!(imgc, imgs, ImageView.AnnotationText(0,0, splitdir(save_xcorr_img_to)[2],
		  		halign="left", valign="top", fontsize=round(minimum(size(overview))/50), color=RGB(0.5,0.5,0.5)))
	  end
  end
  return offsets, overlay
end
