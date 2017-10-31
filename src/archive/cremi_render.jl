"""
Reverts transform that went from index and returns the image at the index
"""
function meshwarp_revert(index::FourTupleIndex, img = get_image(nextstage(index)), interp = false)
  mesh = load("Mesh", index)
  src_nodes = get_nodes(mesh; use_post = false)
  dst_nodes = get_nodes(mesh; use_post = true)
  offset = get_offset(nextstage(index));
  #=print("incidence_to_dict: ")
  @time node_dict = incidence_to_dict(mesh.edges') #'
  print("dict_to_triangles: ")
  @time triangles = dict_to_triangles(node_dict)=#
  @time reverted_img, reverted_offset = meshwarp(img, dst_nodes, src_nodes, incident_to_triangles(mesh.edges), offset, interp)
  original_image_size = get_image_size(index)
  original_offset = get_offset(index)
  i_range = (1:original_image_size[1]) - reverted_offset[1] + original_offset[1] 
  j_range = (1:original_image_size[2]) - reverted_offset[2] + original_offset[2] 
  @inbounds return reverted_img[i_range, j_range]
end

"""
One-off code for CREMI
"""
function segmentation_revert(fn, savepath)
	segs = h5read(fn, "main")
	#cutout_range_i = 401:2900
	#cutout_range_j = 201:2700
	cutout_range_i = 651:2650
	cutout_range_j = 451:2450

	segs_reverted = zeros(UInt32, 1250, 1250, 125)
	for i in 1:125
	  	sec_num = i + 37
		img = zeros(UInt32, get_image_size(aligned(1,sec_num))...)
		offset = get_offset(aligned(1,sec_num))
		img[cutout_range_i - offset[1] + 1, cutout_range_j - offset[2] + 1] = segs[:,:,i]
		img_reverted = meshwarp_revert(prealigned(1,sec_num), img)
		segs_reverted[:,:,i] = img_reverted[912:2161, 912:2161]
	end
		f = h5open(savepath, "w"); @time f["main"] = segs_reverted; 
		close(f)
	return segs_reverted

end