type Match
  src_index          # source mesh index
  dst_index          # destination mesh index

  src_points::Points	    # source point coordinate in image
  dst_points::Points        # destination point coordinate in image
  correspondence_properties::Array{Dict{Any, Any}} # in the same array order as points, contains properties of correspondences

  filters::Array{Dict{Any, Any}}    # Array of filters 
end

function count_correspondences(match::Match) return size(match.src_points, 1);	end

function get_correspondence_patches(match::Match, ind)
	src_img = get_image(match.src_index);
	dst_img = get_image(match.dst_index);
	src_patch = src_img[match.correspondence_properties[ind]["src_range"]...]
	src_pt = match.correspondence_properties[ind]["src_pt_loc"]
	dst_patch = dst_img[match.correspondence_properties[ind]["dst_range"]...]
	dst_pt = match.correspondence_properties[ind]["dst_pt_loc"]
	xc = normxcorr2(src_patch, dst_patch);

	return src_patch, src_pt, dst_patch, dst_pt, xc
end

function get_ranges(pt, src_mesh, dst_mesh, params)
	get_ranges(pt, src_mesh, dst_mesh, params["block_size"], params["search_r"]);
end

function get_ranges(pt, src_index, dst_index, block_size::Int64, search_r::Int64)
	# convert to local coordinates in both src / dst images, and then round up to an integer
	src_pt = ceil(Int64, pt);
	dst_pt = pt + get_offset(src_index) - get_offset(dst_index);
	dst_pt = ceil(Int64, dst_pt);

	block_range = -block_size:block_size;
	search_range = -(block_size+search_r):(block_size+search_r);

	src_range = src_pt[1] + block_range, src_pt[2] + block_range;
	dst_range = dst_pt[1] + search_range, dst_pt[2] + search_range;

	src_img_size = get_image_sizes(src_index);
	dst_img_size = get_image_sizes(dst_index);
	
	range_in_src = intersect(src_range[1], 1:src_img_size[1]), intersect(src_range[2], 1:src_img_size[2]);
	range_in_dst = intersect(dst_range[1], 1:dst_img_size[1]), intersect(dst_range[2], 1:dst_img_size[2]);

	src_pt_locs = findfirst(range_in_src[1] .== src_pt[1]), findfirst(range_in_src[2] .== src_pt[2]);
	dst_pt_locs = findfirst(range_in_dst[1] .== dst_pt[1]), findfirst(range_in_dst[2] .== dst_pt[2]);

	if src_pt_locs[1] == 0 || src_pt_locs[2] == 0 || dst_pt_locs[1] == 0 || dst_pt_locs[2] == 0 
	return nothing
	end

	return src_index, range_in_src, [src_pt_locs[1][1], src_pt_locs[2][1]], dst_index, range_in_dst, [dst_pt_locs[1][1], dst_pt_locs[2][1]];
end

#using sharedarray
function get_match(pt, ranges, src_image, dst_image, params)
	src_index, src_range, src_pt_loc, dst_index, dst_range, dst_pt_loc = ranges;
	
	xc = normxcorr2(src_image[src_range[1], src_range[2]], dst_image[dst_range[1], dst_range[2]]);
	r_max = maximum(xc)
	if isnan(r_max) return nothing end;
  	ind = findfirst(r_max .== xc)
	i_max, j_max = rem(ind, size(xc, 1)), cld(ind, size(xc, 1));
  	if i_max == 0 
    		i_max = size(xc, 1)
  	end

	di = i_max + src_pt_loc[1] - dst_pt_loc[1];	
	dj = j_max + src_pt_loc[2] - dst_pt_loc[2];	

	correspondence_properties = Dict{Any, Any}();
	correspondence_properties["src_range"] = src_range;
	correspondence_properties["src_pt_loc"] = src_pt_loc;
	correspondence_properties["dst_range"] = dst_range;
	correspondence_properties["dst_pt_loc"] = dst_pt_loc;
	correspondence_properties["di"] = di;
	correspondence_properties["dj"] = dj;
	correspondence_properties["r_val"] = r_max;
	return vcat(pt + get_offset(src_index) - get_offset(dst_index) + [di; dj], correspondence_properties);
end

function get_r_val(correspondence_property)
	return correspondence_property["r_val"];
end

function get_norm(correspondence_property)
	return norm(correspondence_property["di"], correspondence_property["dj"]);
end

function filter_r_val!(match::Match, min_r)
	r_vals = map(get_r_val, match.correspondence_properties);
	filtered_inds = find(i -> i < min_r, r_vals);
	push!(match.filters, Dict{Any, Any}(
				"by"	  => ENV["USER"],
				"type"	  => "r_val, absolute",
				"threshold" => min_r,
				"timestamp" => string(now()),
				"rejected"  => filtered_inds
			      ));
	println("$(length(filtered_inds)) / $(count_correspondences(match)) rejected.");
	return;
end



function filter_norm_absolute!(match::Match, max_norm)
	norms = map(get_norm, match.correspondence_properties);
	filtered_inds = find(i -> i > max_norm, norms);
	push!(match.filters, Dict{Any, Any}(
				"by"	  => ENV["USER"],
				"type"	  => "max_norm, absolute",
				"threshold" => max_norm,
				"timestamp" => string(now()),
				"rejected"  => filtered_inds
			      ));
	println("$(length(filtered_inds)) / $(count_correspondences(match)) rejected.");
	return;
end

### ADD MANUAL FILTER
function filter_manual!(match::Match)
	filtered_inds = Points(0)
	while(true)
		#choose point and add
		ind_to_remove = 0;
		push!(filtered_inds, ind_to_remove)
		#undo by pop!
	end
	push!(match.filters, Dict{Any, Any}(
				"by"	  => ENV["USER"],
				"type"	  => "manual",
				"timestamp" => string(now()),
				"rejected"  => filtered_inds
			      ));
	return;
end

function undo_filter!(match::Match)
	pop!(match.filters);
end

function count_filtered_correspondences(match::Match)
	return length(get_filtered_indices(match));
end

function get_filtered_correspondences(match::Match)
	return match.src_points[get_filtered_indices(match)], match.dst_points[get_filtered_indices(match)];
end

function get_filtered_correspondence_properties(match::Match)
	return match.correspondence_properties[get_filtered_indices(match)]
end

function get_filtered_indices(match::Match)
	return setdiff(1:count_correspondences(match), union(map(getindex, match.filters, repeated("rejected"))...));
end

function Match(src_mesh::Mesh, dst_mesh::Mesh, src_image=nothing, dst_image=nothing, params=get_params(src_mesh))
	if src_mesh == dst_mesh
		return nothing
	end

	src_image_local = src_image;
	dst_image_local = dst_image;

	if src_image == nothing src_image_local = get_image(src_mesh); end
	if dst_image == nothing dst_image_local = get_image(dst_mesh); end

	LOCAL_SRC_IMAGE = SharedArray(eltype(src_image_local), size(src_image_local), pids=local_procs());
	LOCAL_SRC_IMAGE[:, :] = src_image_local;
	LOCAL_DST_IMAGE = SharedArray(eltype(dst_image_local), size(dst_image_local), pids=local_procs());
	LOCAL_DST_IMAGE[:, :] = dst_image_local;

  	src_index = src_mesh.index;
	dst_index = dst_mesh.index;

	@time ranges = pmap(get_ranges, src_mesh.src_nodes, repeated(src_index), repeated(dst_index), repeated(params));
	ranged_inds = find(i -> i != nothing, ranges);
	src_nodes = copy(src_mesh.src_nodes[ranged_inds]);
	ranges = copy(ranges[ranged_inds]);
	println("Matching $(src_mesh.index) -> $(dst_mesh.index): $(length(src_nodes)) / $(length(src_mesh.src_nodes)) nodes to check.")

	dst_allpoints = pmap(get_match, src_nodes, ranges, repeated(LOCAL_SRC_IMAGE), repeated(LOCAL_DST_IMAGE), repeated(params));
	matched_inds = find(i -> i != nothing, dst_allpoints);
	src_points = copy(src_nodes[matched_inds]);
	filters = Array{Dict{Any, Any}}(0);

	dst_points = [convert(Point, dst_allpoints[ind][1:2]) for ind in matched_inds]
	correspondence_properties = [dst_allpoints[ind][3] for ind in matched_inds]

	return Match(src_index, dst_index, src_points, dst_points, correspondence_properties, filters);
end

function sync_images(src_image_ref, dst_image_ref)
	src_image_local = fetch(src_image_ref);
	dst_image_local = fetch(dst_image_ref);
	global LOCAL_SRC_IMAGE = SharedArray(eltype(src_image_local), size(src_image_local), pids=local_procs());
	global LOCAL_DST_IMAGE = SharedArray(eltype(dst_image_local), size(dst_image_local), pids=local_procs());
	LOCAL_SRC_IMAGE[:, :] = src_image_local[:, :];
	LOCAL_DST_IMAGE[:, :] = dst_image_local[:, :];

	for pid in local_procs()
	remotecall(pid, sync_images_subroutine, LOCAL_SRC_IMAGE, LOCAL_DST_IMAGE);
      end
#=	tofetch = Array{RemoteRef}(0);
	for pid in local_procs()
		push!(tofetch, remotecall(pid, sync_images_subroutine, LOCAL_SRC_IMAGE, LOCAL_DST_IMAGE));
	end
      	for ref in tofetch
		fetch(ref);
	end
=#
end

function sync_images_subroutine(local_src_image, local_dst_image)
	global LOCAL_SRC_IMAGE = local_src_image;
	global LOCAL_DST_IMAGE = local_dst_image;
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
