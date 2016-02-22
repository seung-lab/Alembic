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

function get_correspondences(match::Match; globalized=false)
	if globalized
	return match.src_points + fill(get_offset(match.src_index), count_correspondences(match)), match.dst_points + fill(get_offset(match.dst_index), count_correspondences(match))
	else
	return match.src_points, match.dst_points;
	end
end
#= DEPRECATED
function get_correspondence_properties(match::Match)
	return match.correspondence_properties;
end
=#

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
function get_properties(match::Match, property_name::String)
	return map(get, match.correspondence_properties, repeated(property_name), repeated(nothing));
end

### reviewing
function set_reviewed!(match::Match)
	if !haskey(match.properties, "review") match.properties["review"] = Dict{Any, Any}(); end
	match.properties["review"]["author"] = author();
	return;
end

function is_reviewed(match::Match)
	isreviewed = false
	if haskey(match.properties, "review")
		if haskey(match.properties["review"], "author")
			isreviewed = true
		end
	end
	return isreviewed
end

function get_author(match::Match)
	author = null_author()
	if haskey(match.properties, "review")
		if haskey(match.properties["review"], "author")
			author = match.properties["review"]["author"]
		end
	end
	return author
end

function is_flagged(match::Match)
	if !haskey(match.properties, "review") 
		match.properties["review"] = Dict{Any, Any}("flag" => false)
	else
		if !haskey(match.properties["review"], "flag") 
			match.properties["review"]["flag"] = false
		end
	end

	return match.properties["review"]["flag"]
end

function flag!(match::Match)
	if !haskey(match.properties, "review") match.properties["review"] = Dict{Any, Any}(); end
	match.properties["review"]["flag"] = true;
end

function unflag!(match::Match)
	if !haskey(match.properties, "review") match.properties["review"] = Dict{Any, Any}(); end
	match.properties["review"]["flag"] = false;
end

function get_correspondence_patches(match::Match, ind)
	src_path = get_path(match.src_index);
	dst_path = get_path(match.dst_index);
	
	assert(src_path[end-1:end] == "h5")

	props = match.correspondence_properties[ind]

	if haskey(props, "scale")
	scale = props["scale"];
	else
	scale = 1
	end

	# hack to support old properties
	if !haskey(props, "full")
		src_patch = h5read(src_path, "img", props["src_range"])
		src_pt = props["src_pt_loc"]
		dst_patch = h5read(dst_path, "img", props["dst_range"])
		dst_pt = props["dst_pt_loc"]
	else
		src_patch = h5read(src_path, "img", props["full"]["src_range"])
		src_pt_loc = props["full"]["src_pt_loc"];
		dst_patch = h5read(dst_path, "img", props["full"]["dst_range"])
		dst_pt_loc = props["full"]["dst_pt_loc"];
		if props["scale"] != 1
			src_pt_loc = ceil(Int64, scale * src_pt_loc);
			dst_pt_loc = ceil(Int64, scale * dst_pt_loc);
			src_patch = imscale(src_patch, scale)[1]
			dst_patch = imscale(dst_patch, scale)[1]
		end
		src_pt = src_pt_loc
		dst_pt = dst_pt_loc
	end

	xc = normxcorr2(src_patch, dst_patch);
	dv = ceil(Int64, props["dv"] * scale)

	return src_patch, src_pt, dst_patch, dst_pt, xc, dst_pt-src_pt+dv
end

### helper methods
function get_ranges(pt, src_mesh, dst_mesh, params)
	get_ranges(pt, src_mesh, dst_mesh, params["match"]["block_r"], params["match"]["search_r"], params["registry"]["global_offsets"]);
end

function get_ranges(pt, src_index, dst_index, block_r::Int64, search_r::Int64, global_offsets = true)
	# convert to local coordinates in both src / dst images, and then round up to an integer
	src_pt = ceil(Int64, pt);
	if global_offsets
	dst_pt = pt + get_offset(src_index) - get_offset(dst_index);
	else
	dst_pt = pt + get_offset(src_index)
	end
	dst_pt = ceil(Int64, dst_pt);

	block_range = -block_r:block_r;
	search_range = -(block_r+search_r):(block_r+search_r);

	src_range_full = src_pt[1] + block_range, src_pt[2] + block_range;
	dst_range_full = dst_pt[1] + search_range, dst_pt[2] + search_range;

	src_img_size = get_image_sizes(src_index);
	dst_img_size = get_image_sizes(dst_index);
	
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

# includes range for this
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

	dv = get_match(dst_pt_locs, ranges, src_image_scaled, dst_image_scaled)[3]["dv"]

	view(src_image_scaled[range_in_src...]/255)
	view(dst_image_scaled[range_in_src...]/255)

	if params["registry"]["global_offsets"]
	update_offset(src_index, get_offset(dst_index) + dv / scale);
	else
	update_offset(src_index, dv / scale);
	end
	print("    ")
	println("Monoblock matched... relative displacement calculated at $( dv / scale)")

end

#using sharedarray
function get_match(pt, ranges, src_image, dst_image, params = nothing)
	src_index, src_range, src_pt_loc, dst_index, dst_range, dst_range_full, dst_pt_loc, dst_pt_loc_full = ranges;

	if params == nothing
		scale = 1;
	else
		scale = params["match"]["blockmatch_scale"]
	end

	correspondence_properties = Dict{Any, Any}();
	correspondence_properties["full"] = Dict{Any, Any}();
	correspondence_properties["full"]["src_pt_loc"] = src_pt_loc;
	correspondence_properties["full"]["src_range"] = src_range;
	correspondence_properties["full"]["dst_pt_loc"] = dst_pt_loc;
	correspondence_properties["full"]["dst_range"] = dst_range;


	src_range = ceil(Int64, scale * first(src_range[1])) : ceil(Int64, scale * last(src_range[1])), ceil(Int64, scale * first(src_range[2])) : ceil(Int64, scale * last(src_range[2]))
	dst_range = ceil(Int64, scale * first(dst_range[1])) : ceil(Int64, scale * last(dst_range[1])), ceil(Int64, scale * first(dst_range[2])) : ceil(Int64, scale * last(dst_range[2]))
	dst_range_full = ceil(Int64, scale * first(dst_range_full[1])) : ceil(Int64, scale * last(dst_range_full[1])), ceil(Int64, scale * first(dst_range_full[2])) : ceil(Int64, scale * last(dst_range_full[2]))
	src_pt_loc = ceil(Int64, scale * src_pt_loc);
	dst_pt_loc = ceil(Int64, scale * dst_pt_loc);
	dst_pt_loc_full = ceil(Int64, scale * dst_pt_loc_full);

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


	correspondence_properties["img"] = Dict{Any, Any}();
	correspondence_properties["img"]["src_normalized_dyn_range"] = (maximum(src_image[src_range...]) - minimum(src_image[src_range...])) / typemax(eltype(src_image));
	correspondence_properties["img"]["src_kurtosis"] = kurtosis(src_image[src_range...]);
	correspondence_properties["scale"] = scale;
	correspondence_properties["dv"] = [di, dj];
	correspondence_properties["norm"] = norm([di, dj]);
	correspondence_properties["r_val"] = r_max;


	if params == nothing || !params["registry"]["global_offsets"]
		return vcat(pt + get_offset(src_index) + [di, dj], correspondence_properties);
	else
		return vcat(pt + get_offset(src_index) - get_offset(dst_index) + [di, dj], correspondence_properties);
	end
end



function filter!(match::Match, property_name, compare, threshold)
	attributes = get_properties(match, property_name)
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
#=
function eval_filter(match::Match, property_name, compare, threshold)
	attributes = get_properties(match, property_name)
	inds_to_filter = find(i -> compare(i, threshold), attributes);
	rejected_inds = get_rejected_indices(match);
	
	if length(inds_to_filter) != 0 filter_reject_match = true;
	else filter_reject_match = false; end

	if length(rejected_inds) != 0 actual_reject_match = true;
	else actual_reject_match = false; end

	false_rejections = setdiff(inds_to_filter, rejected_inds)
	false_acceptances = setdiff(rejected_inds, inds_to_filter)
	common_rejections = intersect(rejected_inds, inds_to_filter)

#=
	println("false_rejections: $false_rejections")
	println("false_acceptances: $false_acceptances")
	println("common_rejections: $common_rejections")
	=#
	return length(false_rejections), length(false_acceptances), length(common_rejections), count_correspondences(match), filter_reject_match, actual_reject_match;
end
=#

function eval_filters(match::Match, filters, conjunction=false)

	inds_to_filter = Array{Any, 1}();
	thresholds = Array{Int64, 1}();

	for filter in filters
	attributes = get_properties(match, filter[1]);
	push!(inds_to_filter, find(i -> filter[2](i, filter[3]), attributes));
	push!(thresholds, filter[4]);
	end

	rejected_inds = get_rejected_indices(match);

	if conjunction == false
		if Base.|((thresholds .< map(length, inds_to_filter))...) filter_reject_match = true
		else filter_reject_match = false; end
	else
		if Base.&((thresholds .< map(length, inds_to_filter))...) filter_reject_match = true
		else filter_reject_match = false; end
	end

	inds_to_filter = union(inds_to_filter...)

	if length(rejected_inds) > 0 actual_reject_match = true;
	else actual_reject_match = false; end

	false_rejections = setdiff(inds_to_filter, rejected_inds)
	false_acceptances = setdiff(rejected_inds, inds_to_filter)
	common_rejections = intersect(rejected_inds, inds_to_filter)

	return length(false_rejections), length(false_acceptances), length(common_rejections), count_correspondences(match), filter_reject_match, actual_reject_match;
end

function eval_filters(matches::Array{Match, 1}, filters, conjunction=false)
	evals = map(eval_filters, matches, repeated(filters), repeated(conjunction))

	total_false_rej= 0;
	total_false_acc = 0;
	total_correct = 0;
	total_corresp = 0;

	match_false_rej = 0;
	match_false_acc = 0;
	match_correct = 0;

	for i in evals
	total_false_rej = total_false_rej + i[1];
	total_false_acc = total_false_acc + i[2];
	total_correct = total_correct + i[3];
	total_corresp = total_corresp + i[4];

	if i[5] == true && i[6] == true match_correct = match_correct + 1; end
	if i[5] == true && i[6] == false match_false_rej = match_false_rej + 1; end
	if i[5] == false && i[6] == true match_false_acc = match_false_acc + 1; end

	end

	total = total_false_acc + total_correct

	println("filtering by: $(filters...)")
#=
	println("false rejections: $(total_false_rej / total)")
	println("false non-rejections: $(total_false_acc / total)")
	println("correct rejections: $(total_correct / total)")
=#
	println("Per seam:")
	println("precision: $(100 * match_correct / (match_false_rej + match_correct)) %")
	println("recall: $(100 * match_correct / (match_correct + match_false_acc)) %")
	println();

	println("Per correspondence:")
	println("precision: $(100 * total_correct / (total_false_rej + total_correct)) %")
	println("recall: $(100 * total_correct / (total_correct + total_false_acc)) %")
	return total_false_rej, total_false_acc, total_correct, total_corresp
end

function eval_filters_meshsets(mses, filters)
	evals = pmap(eval_filters, mses, repeated(filters))

	total_false_rej = 0;
	total_false_acc = 0;
	total_correct = 0;

	for i in evals
	total_false_rej = total_false_rej + i[1];
	total_false_acc = total_false_acc + i[2];
	total_correct = total_correct + i[3];
	total_corresp = total_corresp + i[4];
	end

	println("filtering by: $(filters...)")
	println("precision: $(100 * total_correct / (total_false_rej + total_correct)) %")
	println("recall: $(100 * total_correct / (total_correct + total_false_acc)) %")
	println("total correspondences: $total_corresp")
	println("total correct: $total_correct")
	println("total falsely rejected: $total_false_rej")
	println("total falsely accepted: $total_false_acc")
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
	pop!(match.filters);
end



function Match(src_mesh::Mesh, dst_mesh::Mesh, params=get_params(src_mesh); src_image=get_image(src_mesh), dst_image=get_image(dst_mesh), rotate=0)
	println("Matching $(src_mesh.index) -> $(dst_mesh.index):")
	if src_mesh == dst_mesh
		println("nothing at")
		println(src_mesh.index)
		println(dst_mesh.index)
		return nothing
	end

	#load at full scale for monoblock match
	SHARED_SRC_IMAGE[1:size(src_image, 1), 1:size(src_image, 2)] = src_image;
	SHARED_DST_IMAGE[1:size(dst_image, 1), 1:size(dst_image, 2)] = dst_image;
  	src_index = src_mesh.index;
	dst_index = dst_mesh.index;

	monoblock_match(src_index, dst_index, src_image, dst_image, params);

	#load at full scale for monoblock match
	scale = params["match"]["blockmatch_scale"];

	if scale != 1
	src_image = imscale(src_image, scale)[1]
	dst_image = imscale(dst_image, scale)[1]
	SHARED_SRC_IMAGE[1:size(src_image, 1), 1:size(src_image, 2)] = src_image;
	SHARED_DST_IMAGE[1:size(dst_image, 1), 1:size(dst_image, 2)] = dst_image;
	end

	ranges = pmap(get_ranges, src_mesh.src_nodes, repeated(src_index), repeated(dst_index), repeated(params));
	ranged_inds = find(i -> i != nothing, ranges);
	src_nodes = copy(src_mesh.src_nodes[ranged_inds]);
	ranges = copy(ranges[ranged_inds]);
	print("    ")
	println("$(length(src_nodes)) / $(length(src_mesh.src_nodes)) nodes to check.")

	dst_allpoints = pmap(get_match, src_nodes, ranges, repeated(SHARED_SRC_IMAGE), repeated(SHARED_DST_IMAGE), repeated(params));
	#dst_allpoints = map(get_match, src_nodes, ranges, repeated(SHARED_SRC_IMAGE), repeated(SHARED_DST_IMAGE), repeated(params));
	matched_inds = find(i -> i != nothing, dst_allpoints);
	src_points = copy(src_nodes[matched_inds]);
	filters = Array{Dict{Any, Any}}(0);

	dst_points = [convert(Point, dst_allpoints[ind][1:2]) for ind in matched_inds]
	correspondence_properties = [dst_allpoints[ind][3] for ind in matched_inds]
  	properties = Dict{Any, Any}(
		"review" => Dict{Any, Any}(
				"flag" => true) 
			);

	return Match(src_index, dst_index, src_points, dst_points, correspondence_properties, filters, properties);
end




function sync_images(src_image_ref, dst_image_ref)
	src_image_local = fetch(src_image_ref);
	dst_image_local = fetch(dst_image_ref);
	global SHARED_SRC_IMAGE = SharedArray(eltype(src_image_local), size(src_image_local), pids=local_procs());
	global SHARED_DST_IMAGE = SharedArray(eltype(dst_image_local), size(dst_image_local), pids=local_procs());
	SHARED_SRC_IMAGE[:, :] = src_image_local[:, :];
	SHARED_DST_IMAGE[:, :] = dst_image_local[:, :];

	for pid in local_procs()
	remotecall(pid, sync_images_subroutine, SHARED_SRC_IMAGE, SHARED_DST_IMAGE);
      end
#=	tofetch = Array{RemoteRef}(0);
	for pid in local_procs()
		push!(tofetch, remotecall(pid, sync_images_subroutine, SHARED_SRC_IMAGE, SHARED_DST_IMAGE));
	end
      	for ref in tofetch
		fetch(ref);
	end
=#
end

function sync_images_subroutine(local_src_image, local_dst_image)
	global SHARED_SRC_IMAGE = local_src_image;
	global SHARED_DST_IMAGE = local_dst_image;
end

function my_host_addr()
	return Base.Worker(myid()).bind_addr;
end

function get_local_host_addr(id)
	return remotecall_fetch(id, my_host_addr);
end

function local_procs()
	localhost = Base.Worker(myid()).bind_addr;
	remotehosts = map(get_local_host_addr, procs());
	local_procs_indices = find(p -> p == localhost, remotehosts);
	return procs()[local_procs_indices];
end

function Base.size(r::RemoteRef, args...)
      if r.where == myid()
	return size(fetch(r), args...)
		    end
	return remotecall_fetch(r.where, size, r, args...)
end

function Base.eltype(r::RemoteRef, args...)
      if r.where == myid()
	return eltype(fetch(r), args...)
		    end
	return remotecall_fetch(r.where, eltype, r, args...)
end
