type Match
  src_index          # source mesh index
  dst_index          # destination mesh index

  src_points::Points	    # source point coordinate in image
  dst_points::Points        # destination point coordinate in image
  correspondence_properties::Array{Dict{Any, Any}} # in the same array order as points, contains properties of correspondences

  filters::Array{Dict{Any, Any}}    # Array of filters 
  properties::Dict{Any, Any}
end

### counting
function count_correspondences(match::Match) return size(match.src_points, 1);	end
function count_filtered_correspondences(match::Match) return length(get_filtered_indices(match)); end
function count_filters(match::Match) return length(match.filters); end

function get_ratio_filtered(match::Match, min_corresps = 0) 
  if count_correspondences(match) < min_corresps return 1.0 end
return count_filtered_correspondences(match) / max(count_correspondences(match), 1); end

function get_ratio_edge_proximity(match::Match)
     if count_filtered_correspondences(match) == 0 return 0.0 end
     norms = map(norm, get_filtered_properties(match, "dv"))
return maximum(norms) / match.properties["params"]["match"]["search_r"]; end

function get_src_index(match::Match)
  return match.src_index
end

function get_dst_index(match::Match)
  return match.dst_index
end

function get_correspondences(match::Match; globalized=false)
	if globalized
	return match.src_points + fill(get_offset(match.src_index), count_correspondences(match)), match.dst_points + fill(get_offset(match.dst_index), count_correspondences(match))
	else
	return match.src_points, match.dst_points;
	end
end

function get_filtered_correspondences(match::Match; globalized=false)
	src_pts, dst_pts = get_correspondences(match; globalized = globalized);
	return src_pts[get_filtered_indices(match)], dst_pts[get_filtered_indices(match)];
end

function get_filtered_correspondence_properties(match::Match)
	return match.correspondence_properties[get_filtered_indices(match)]
end

function get_filtered_indices(match::Match)
	return setdiff(1:count_correspondences(match), union(map(getindex, match.filters, repeated("rejected"))...));
end

function get_rejected_indices(match::Match)
	return union(map(getindex, match.filters, repeated("rejected"))...)
end

### property handling
function get_properties(match::Match, property_name)
	ret = map(get_dfs, match.correspondence_properties, repeated(property_name));
	if length(ret) != 0
		ret = Array{typeof(ret[1])}(ret);
	end
end

function get_properties(match::Match, fn::Function, args...)
	return fn(match, args...)
end

function get_dfs(dict::Dict, keytofetch)
	for key in keys(dict)
		if key == keytofetch
		  return get(dict, key, nothing)
		elseif typeof(get(dict, key, nothing)) <: Dict
		  if get_dfs(get(dict, key, nothing), keytofetch) != nothing 
		  	return get_dfs(get(dict, key, nothing), keytofetch) 
		  end
		end
	end
	return nothing
end

function get_filtered_properties(match::Match, property_name)
	cp = get_filtered_correspondence_properties(match)
	return map(get_dfs, cp, repeated(property_name));
end

### reviewing
function set_reviewed!(match::Match)
	match.properties["review"]["author"] = author();
	return;
end

function is_reviewed(match::Match)
  	return match.properties["review"]["author"] == null_author()
end

function get_author(match::Match)
	author = match.properties["review"]["author"]
	return author
end

function is_flagged(match::Match)
	flagged = false
	if haskey(match.properties["review"], "flagged")
		flagged = match.properties["review"]["flagged"]
	end
	return flagged
end

function flag!(match::Match)
	match.properties["review"]["flagged"] = true;
end

"""
Flag a match based on property criteria
"""
function flag!(match::Match, property_name, compare, threshold)
	attributes = get_filtered_properties(match, property_name)
	inds_to_filter = find(i -> compare(i, threshold), attributes)
	if length(inds_to_filter) > 0
		flag!(match)
	end
end

function unflag!(match::Match)
	match.properties["review"]["flagged"] = false;
end

function get_correspondence_patches(match::Match, ind)
	src_path = get_path(match.src_index);
	dst_path = get_path(match.dst_index);
	
	#assert(src_path[end-1:end] == "h5")

	props = match.correspondence_properties[ind]

		scale = props["ranges"]["scale"];
		src_patch = h5read(src_path, "img", props["ranges"]["src_range"])
		src_pt_loc = props["ranges"]["src_pt_loc"];
		dst_patch = h5read(dst_path, "img", props["ranges"]["dst_range"])
		dst_pt_loc = props["ranges"]["dst_pt_loc"];
		if scale != 1
			src_pt_loc = ceil(Int64, scale * src_pt_loc);
			dst_pt_loc = ceil(Int64, scale * dst_pt_loc);
			src_patch = imscale(src_patch, scale)[1]
			dst_patch = imscale(dst_patch, scale)[1]
		end
		src_pt = src_pt_loc
		dst_pt = dst_pt_loc

	xc = normxcorr2(src_patch, dst_patch);
	dv = ceil(Int64, props["vects"]["dv"] * scale)

	return src_patch, src_pt, dst_patch, dst_pt, xc, dst_pt-src_pt+dv
end

### helper methods
function get_ranges(pt, src_index, src_offset, src_img_size, dst_index, dst_offset, dst_img_size, params)
	get_ranges(pt, src_index, dst_index, params["match"]["block_r"], params["match"]["search_r"], params["registry"]["global_offsets"]);
end

function get_ranges(pt, src_index, src_offset, src_img_size, dst_index, dst_offset, dst_img_size, block_r::Int64, search_r::Int64, global_offsets = true)
	# convert to local coordinates in both src / dst images, and then round up to an integer
	src_pt = ceil(Int64, pt);
	if global_offsets
	dst_pt = pt + src_offset - dst_offset
	else
	dst_pt = pt + src_offset
	end
	dst_pt = ceil(Int64, dst_pt);

	block_range = -block_r:block_r;
	search_range = -(block_r+search_r):(block_r+search_r);

	src_range_full = src_pt[1] + block_range, src_pt[2] + block_range;
	dst_range_full = dst_pt[1] + search_range, dst_pt[2] + search_range;

	
	range_in_src = intersect(src_range_full[1], 1:src_img_size[1]), intersect(src_range_full[2], 1:src_img_size[2]);
	range_in_dst = intersect(dst_range_full[1], 1:dst_img_size[1]), intersect(dst_range_full[2], 1:dst_img_size[2]);

	src_pt_locs = findfirst(range_in_src[1] .== src_pt[1]), findfirst(range_in_src[2] .== src_pt[2]);
	dst_pt_locs = findfirst(range_in_dst[1] .== dst_pt[1]), findfirst(range_in_dst[2] .== dst_pt[2]);
	dst_pt_locs_full = findfirst(dst_range_full[1] .== dst_pt[1]), findfirst(dst_range_full[2] .== dst_pt[2]);

	if src_pt_locs[1] == 0 || src_pt_locs[2] == 0 || dst_pt_locs[1] == 0 || dst_pt_locs[2] == 0 
	return nothing
	end


	return src_index, range_in_src, [src_pt_locs[1], src_pt_locs[2]], dst_index, range_in_dst, dst_range_full, [dst_pt_locs[1], dst_pt_locs[2]], [dst_pt_locs_full[1], dst_pt_locs_full[2]];
end

"""
Template match two images & record translation for source image
"""
function monoblock_match(src_index, dst_index, src_image, dst_image, params=get_params(src_index))
	if params["match"]["monoblock_match"] == false return; end
	scale = params["match"]["monoblock_scale"];
	ratio = params["match"]["monoblock_ratio"];
	src_image_scaled = imscale(src_image, scale)[1]
	dst_image_scaled = imscale(dst_image, scale)[1]

	scaled_rads = ceil(Int64, ratio * size(src_image_scaled, 1) / 2), round(Int64, ratio * size(src_image_scaled, 2) / 2)
	range_in_src = ceil(Int64, size(src_image_scaled, 1) / 2) + (-scaled_rads[1]:scaled_rads[1]), round(Int64, size(src_image_scaled, 2) / 2) + (-scaled_rads[2]:scaled_rads[2])

	range_in_dst = 1:size(dst_image_scaled, 1), 1:size(dst_image_scaled, 2);
	dst_range_full = 1:size(dst_image_scaled, 1), 1:size(dst_image_scaled, 2);
	src_pt_locs = [1, 1]; 
	dst_pt_locs = [first(range_in_src[1]), first(range_in_src[2])];
	dst_pt_locs_full = dst_pt_locs;

	ranges = src_index, range_in_src, src_pt_locs, dst_index, range_in_dst, dst_range_full, dst_pt_locs, dst_pt_locs_full

	dv = get_match(src_pt_locs, ranges, src_image_scaled, dst_image_scaled)[3]["dv"]

	#=view(src_image_scaled[range_in_src...]/255)
	view(dst_image_scaled[range_in_dst...]/255)
	println(dst_pt_locs + dv)
	img = normxcorr2(src_image_scaled[range_in_src...], dst_image_scaled[range_in_dst...]);
	view(img / maximum(img)) =#

	if params["registry"]["global_offsets"]
	update_offset(src_index, get_offset(dst_index) + dv / scale);
	else
	update_offset(src_index, (dv / scale));
	end
	print("    ")
	println("Monoblock matched... relative displacement calculated at $( dv / scale)")

end

"""
Template match two image patches to produce point pair correspondence
"""
function get_match(pt, ranges, src_image, dst_image, params)
	src_index, src_range, src_pt_loc, dst_index, dst_range, dst_range_full, dst_pt_loc, dst_pt_loc_full = ranges;

	scale = params["match"]["blockmatch_scale"]

	correspondence_properties = Dict{Any, Any}();
	correspondence_properties["ranges"] = Dict{Any, Any}();
	correspondence_properties["ranges"]["src_pt_loc"] = src_pt_loc;
	correspondence_properties["ranges"]["src_range"] = src_range;
	correspondence_properties["ranges"]["dst_pt_loc"] = dst_pt_loc;
	correspondence_properties["ranges"]["dst_range"] = dst_range;
	correspondence_properties["ranges"]["scale"] = scale;


	src_range = ceil(Int64, scale * first(src_range[1])) : ceil(Int64, scale * last(src_range[1])), ceil(Int64, scale * first(src_range[2])) : ceil(Int64, scale * last(src_range[2]))
	dst_range = ceil(Int64, scale * first(dst_range[1])) : ceil(Int64, scale * last(dst_range[1])), ceil(Int64, scale * first(dst_range[2])) : ceil(Int64, scale * last(dst_range[2]))
	dst_range_full = ceil(Int64, scale * first(dst_range_full[1])) : ceil(Int64, scale * last(dst_range_full[1])), ceil(Int64, scale * first(dst_range_full[2])) : ceil(Int64, scale * last(dst_range_full[2]))
	src_pt_loc = ceil(Int64, scale * src_pt_loc);
	dst_pt_loc = ceil(Int64, scale * dst_pt_loc);
	dst_pt_loc_full = ceil(Int64, scale * dst_pt_loc_full);

	# padding
	if dst_range != dst_range_full
		indices_within_range = findin(dst_range_full[1], dst_range[1]), findin(dst_range_full[2], dst_range[2])
		intersect_img = dst_image[dst_range...];
		avg = mean(intersect_img);
		if isnan(avg) return nothing; end
		avg = round(eltype(dst_image), avg);
		padded_img = fill(avg, length(dst_range_full[1]), length(dst_range_full[2]));
		padded_img[indices_within_range...] = intersect_img;
		xc = normxcorr2(src_image[src_range[1], src_range[2]], padded_img);
	else
	xc = normxcorr2(src_image[src_range[1], src_range[2]], dst_image[dst_range[1], dst_range[2]]);
	end

	r_max = maximum(xc)
	if isnan(r_max) return nothing end;
  	ind = findfirst(r_max .== xc)
	i_max, j_max = rem(ind, size(xc, 1)), cld(ind, size(xc, 1));
  	if i_max == 0 
    		i_max = size(xc, 1)
  	end

	di = (i_max + src_pt_loc[1] - dst_pt_loc_full[1]) / scale;
	dj = (j_max + src_pt_loc[2] - dst_pt_loc_full[2]) / scale;


	#println("src_range: $src_range, dst_range: $dst_range")
	#println("i_max, j_max: $i_max, $j_max, src_pt_loc: $src_pt_loc, dst_pt_loc: $dst_pt_loc, dst_pt_loc_full = $dst_pt_loc_full")
	#println("di: $di, dj: $dj, scale = $scale")


	correspondence_properties["patches"] = Dict{Any, Any}();
	correspondence_properties["patches"]["src_normalized_dyn_range"] = (maximum(src_image[src_range...]) - minimum(src_image[src_range...])) / typemax(eltype(src_image));
	correspondence_properties["patches"]["src_kurtosis"] = kurtosis(src_image[src_range...]);
	correspondence_properties["xcorr"] = Dict{Any, Any}();
	correspondence_properties["xcorr"]["r_max"] = r_max;
	correspondence_properties["xcorr"]["sigmas"] = Dict{Any, Any}();
	correspondence_properties["xcorr"]["sigmas"][.5] = sigma(xc, .5);
	correspondence_properties["xcorr"]["sigmas"][.6] = sigma(xc, .6);
	correspondence_properties["xcorr"]["sigmas"][.7] = sigma(xc, .7);
	correspondence_properties["xcorr"]["sigmas"][.8] = sigma(xc, .8);
	correspondence_properties["vects"] = Dict{Any, Any}();
	correspondence_properties["vects"]["dv"] = [di, dj];
	correspondence_properties["vects"]["norm"] = norm([di, dj]);


	if !params["registry"]["global_offsets"]
		return vcat(pt + get_offset(src_index) + [di, dj], correspondence_properties);
	else
		return vcat(pt + get_offset(src_index) - get_offset(dst_index) + [di, dj], correspondence_properties);
	end
end

function filter!(match::Match, property_name, compare, threshold)
	attributes = get_properties(match, property_name)
	if attributes == nothing return 0; end
	inds_to_filter = find(i -> compare(i, threshold), attributes);
	push!(match.filters, Dict{Any, Any}(
				"author" => author(),
				"type"	  => property_name,
				"threshold" => threshold,
				"rejected"  => inds_to_filter
			      ));
	#println("$(length(inds_to_filter)) / $(count_correspondences(match)) rejected.");
	return length(inds_to_filter);
end

function filter!(match::Match, filters) 
     for filter in filters
        filter!(match, filter...);
     end
end

#### HACKY
function get_residual_norms_post(match, ms)
	src_pts_after, dst_pts_after, filtered = get_globalized_correspondences_post(ms, findfirst(match_in_ms -> match_in_ms.src_index == match.src_index && match_in_ms.dst_index == match.dst_index, ms.matches));
	return(map(norm, dst_pts_after - src_pts_after))
end

function check(match::Match, function_name, compare, threshold, vars...)
     return compare(eval(function_name)(match, vars...), threshold)
end

function check!(match::Match, crits) 
     for crit in crits
       if check(match, crit...) flag!(match); return true; end
     end
     unflag!(match);
     return false
end

### ADD MANUAL FILTER
function filter_manual!(match::Match, inds_to_filter; filtertype="manual")
	push!(match.filters, Dict{Any, Any}(
				"author"	  => author(),
				"type"	  => filtertype,
				"rejected"  => inds_to_filter
			      ));
	return;
end

function clear_filters!(match::Match; filtertype=nothing)
	match.filters = match.filters[setdiff(1:length(match.filters), find(filter -> filter["type"] == filtertype, match.filters))]
	if filtertype == nothing	match.filters = Array{Dict{Any, Any}, 1}(); end

end

function undo_filter!(match::Match)
	if length(match.filters) > 0
		pop!(match.filters);
	end
end

function Match(src_mesh::Mesh, dst_mesh::Mesh, params=get_params(src_mesh); src_image=get_image(src_mesh), dst_image=get_image(dst_mesh), rotate=0)
	println("Matching $(get_index(src_mesh)) -> $(get_index(dst_mesh)):")
	if src_mesh == dst_mesh
		println("nothing at")
		println(get_index(src_mesh))
		println(get_index(dst_mesh))
		return nothing
	end

	#load at full scale for monoblock match
	SHARED_SRC_IMAGE[1:size(src_image, 1), 1:size(src_image, 2)] = src_image;
	SHARED_DST_IMAGE[1:size(dst_image, 1), 1:size(dst_image, 2)] = dst_image;
  	src_index = get_index(src_mesh);
	dst_index = get_index(dst_mesh);

	monoblock_match(src_index, dst_index, src_image, dst_image, params);

	#load at full scale for monoblock match
	scale = params["match"]["blockmatch_scale"];

	if scale != 1
	src_image = imscale(src_image, scale)[1]
	dst_image = imscale(dst_image, scale)[1]
	SHARED_SRC_IMAGE[1:size(src_image, 1), 1:size(src_image, 2)] = src_image;
	SHARED_DST_IMAGE[1:size(dst_image, 1), 1:size(dst_image, 2)] = dst_image;
	end

	ranges = pmap(get_ranges, src_mesh.src_nodes, repeated(src_index), repeated(get_offset(src_index)), repeated(get_image_sizes(src_index)), repeated(dst_index), repeated(get_offset(dst_index)), repeated(get_image_sizes(dst_index)), repeated(params["match"]["block_r"]), repeated(params["match"]["search_r"]), repeated(params["registry"]["global_offsets"]));
	ranged_inds = find(i -> i != nothing, ranges);
	#src_nodes = copy(src_mesh.src_nodes[ranged_inds]);
	ranges = ranges[ranged_inds];
	print("    ")
	println("$(length(ranged_inds)) / $(length(src_mesh.src_nodes)) nodes to check.")
	if length(ranged_inds) != 0
		ranges = Array{typeof(ranges[1]), 1}(ranges);
	end

	dst_allpoints = pmap(get_match, src_mesh.src_nodes[ranged_inds], ranges, repeated(SHARED_SRC_IMAGE), repeated(SHARED_DST_IMAGE), repeated(params));
	#dst_allpoints = map(get_match, src_nodes, ranges, repeated(SHARED_SRC_IMAGE), repeated(SHARED_DST_IMAGE), repeated(params));
	matched_inds = find(i -> i != nothing, dst_allpoints);
	src_points = copy(src_mesh.src_nodes[ranged_inds][matched_inds]);
	filters = Array{Dict{Any, Any}}(0);

	dst_points = [convert(Point, dst_allpoints[ind][1:2]) for ind in matched_inds]
	correspondence_properties = [dst_allpoints[ind][3] for ind in matched_inds]
  	properties = Dict{Any, Any}(
		"params" => params,
		"review" => Dict{Any, Any}(
				"flagged" => false,
				"flags" => Dict{Any, Any},
				"author" => null_author()
				) 
			);

	return Match(src_index, dst_index, src_points, dst_points, correspondence_properties, filters, properties);
end
