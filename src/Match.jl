type Match
  src_index          # source mesh index
  dst_index          # destination mesh index

  src_points::Points	    # source point coordinate in image
  dst_points::Points        # destination point coordinate in image
  properties::Array{Dict{Any, Any}} # in the same array order, contains parameters
end

function count_correspondences(match::Match) return size(src_points, 1);	end

function get_ranges(pt, src_mesh, dst_mesh, params)
	get_ranges(pt, src_mesh, dst_mesh, params["block_size"], params["search_r"]);
end

function get_ranges(pt, src_mesh, dst_mesh, block_size::Int64, search_r::Int64)
	# convert to local coordinates in both src / dst images, and then round up to an integer
	src_pt = ceil(Int64, pt);
	dst_pt = pt + get_offsets(src_mesh) - get_offsets(dst_mesh);
	dst_pt = ceil(Int64, dst_pt);

	block_range = -block_size:block_size;
	search_range = -(block_size+search_r):(block_size+search_r);
	search_range_padded = -(block_size+search_r):(block_size+search_r+1);

	src_range = src_pt[1] + block_range, src_pt[2] + block_range;
	dst_range = dst_pt[1] + search_range, dst_pt[2] + search_range;
	dst_range_padded = dst_pt[1] + search_range_padded, dst_pt[2] + search_range_padded;
	
	return src_range, dst_range, dst_range_padded;
end


# retrieve the ranges, only works with UInt8
function get_patch(img, range)
	range_in_img = intersect(range[1], 1:size(img)[1]), intersect(range[2], 1:size(img)[2]);

	# if the image is fully contained, then return the intersection - otherwise find location and pad by mean
	if length(range_in_img[1]) == length(range[1]) && length(range_in_img[2]) == length(range[2])
		return img[range_in_img...];
	else
		indices_within_range = findin(range[1], range_in_img[1]), findin(range[2], range_in_img[2])
		intersect_img = img[range_in_img...];
		avg = mean(intersect_img);
		if isnan(avg) return fill(UInt8(0), length(range[1]), length(range[2])); end
		avg = round(UInt8, avg);
		padded_img = fill(avg, length(range[1]), length(range[2]));
		padded_img[indices_within_range...] = intersect_img;
		return padded_img;
	end
end

function get_match(pt, src_mesh, src_image, dst_mesh, dst_image, params)
	src_range, dst_range, dst_range_padded = get_ranges(pt, src_mesh, dst_mesh, params);
	xc = normxcorr2(get_patch(src_image, src_range), get_patch(dst_image, dst_range));
	r_max = maximum(xc)
	if isnan(r_max) return nothing end;
  	ind = findfirst(r_max .== xc)
	i_max, j_max = rem(ind, size(xc, 1)), cld(ind, size(xc, 1));
  	if i_max == 0 
    		i_max = size(xc, 1)
  	end
  	rad_i = round(Int64, (size(xc, 1) - 1)/ 2)  
  	rad_j = round(Int64, (size(xc, 2) - 1)/ 2)  
	di, dj = i_max - 1 - rad_i, j_max - 1 - rad_j;
	properties = Dict{Any, Any}();
	properties["di"] = di;
	properties["dj"] = dj;
	properties["r_val"] = r_max;

	return vcat(pt + get_offsets(src_mesh) - get_offsets(dst_mesh) + [i_max - 1 - rad_i; j_max - 1 - rad_j], properties);

end

#=
"""
Write out blockmatches during the blockmatching process
"""
function write_blockmatches(A, B, xc, idx, partial_fn)        
  #=imwrite(grayim((A[idx]/255)'), string(partial_fn, "_src.jpg"))
  imwrite(grayim((B[idx]/255)'), string(partial_fn, "_dst.jpg"))
  if (!isnan(sum(xc[idx])))   
    imwrite(grayim(xc[idx]'), string(partial_fn, "_xc.jpg"))
  end=#
  imwrite(grayim((A[idx]/255)'), string(partial_fn, "_src.tif"))
  imwrite(grayim((B[idx]/255)'), string(partial_fn, "_dst.tif"))
  if (!isnan(sum(xc[idx])))   
    imwrite(grayim(xc[idx]'), string(partial_fn, "_xc.tif"))
  end
end
=#
function Match(src_mesh::Mesh, dst_mesh::Mesh, params=get_params(src_mesh))
	if src_mesh == dst_mesh
		return nothing
	end

	println("Matching $(src_mesh.index) -> $(dst_mesh.index):#########\n");

	src_image_local = get_image(src_mesh);
	dst_image_local = get_image(dst_mesh);

	src_image = SharedArray(eltype(src_image_local), size(src_image_local));
	dst_image = SharedArray(eltype(dst_image_local), size(dst_image_local));
	src_image[:, :] = src_image_local[:, :];
	dst_image[:, :] = dst_image_local[:, :];
  	src_index = src_mesh.index;
	dst_index = dst_mesh.index;

	dst_allpoints = pmap(get_match, src_mesh.src_nodes, repeated(src_mesh), repeated(src_image), repeated(dst_mesh), repeated(dst_image), repeated(params));
	matched_inds = find(i -> i != nothing, dst_allpoints);
	src_points = copy(src_mesh.src_nodes[matched_inds]);
	dst_points = similar(src_points, 0);
	properties = Array{Dict{Any, Any}}(0);

	for ind in matched_inds
		push!(dst_points, convert(Point, dst_allpoints[ind][1:2]));
		push!(properties, dst_allpoints[ind][3]);
	end

	return Match(src_index, dst_index, src_points, dst_points, properties);
end
