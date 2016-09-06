global SRC_PATCH_FULL = Array{Float64, 2}(10,10);
global SRC_PATCH = Array{Float64, 2}(10,10);
global SRC_PATCH_G = Array{Float64, 2}(10,10);
global DST_PATCH_FULL = Array{Float64, 2}(10,10);
global DST_PATCH = Array{Float64, 2}(10,10);
global DST_PATCH_G = Array{Float64, 2}(10,10);

type Match
  src_index          # source mesh index
  dst_index          # destination mesh index

  src_points::Points	    # source point coordinate in image
  dst_points::Points        # destination point coordinate in image
  correspondence_properties::Array{Dict{Any, Any}} # in the same array order as points, contains properties of correspondences

  filters::Array{Dict{Any, Any}}    # Array of filters 
  properties::Dict{Any, Any}
end

### index related
function get_src_index(match::Match)
  return match.src_index
end
function get_dst_index(match::Match)
  return match.dst_index
end
function get_index(match::Match)
	return get_src_and_dst_indices(match::Match)
end
function get_src_and_dst_indices(match::Match)
  return match.src_index, match.dst_index
end

### counting
function count_correspondences(match::Match) return size(match.src_points, 1);	end
function count_filtered_correspondences(match::Match) return length(get_filtered_indices(match)); end
function count_rejected_correspondences(match::Match) return length(get_rejected_indices(match)); end

function count_filters(match::Match) return length(match.filters); end
function count_filtered_properties(match::Match, property_name, compare, threshold)
 return sum(compare(get_filtered_properties(match::Match, property_name), threshold))
end

### correspondences related
function get_filtered_indices(match::Match)
	return setdiff(1:count_correspondences(match), union(map(getindex, match.filters, repeated("rejected"))...));
end
function get_rejected_indices(match::Match)
	return union(map(getindex, match.filters, repeated("rejected"))...)
end
function get_filtered_correspondence_properties(match::Match)
	return match.correspondence_properties[get_filtered_indices(match)]
end

function get_correspondences(match::Match; globalized::Bool=false, global_offsets::Bool=true, filtered::Bool=false, use_post::Bool=false, src_mesh = nothing, dst_mesh = nothing)
  	if filtered
	  src_pts = match.src_points[get_filtered_indices(match)]; dst_pts = match.dst_points[get_filtered_indices(match)];
	else
  	  src_pts = copy(match.src_points); dst_pts = copy(match.dst_points);
	end
	if use_post
	  src_pts = deform(src_pts, src_mesh);
	  dst_pts = deform(dst_pts, dst_mesh);
	end
	if globalized
	  globalize!(src_pts, get_offset(match.src_index)); 
	  global_offsets ? globalize!(dst_pts, get_offset(match.dst_index)) : nothing
	end
	filtered ? ret = (src_pts, dst_pts) : ret = (src_pts, dst_pts, get_filtered_indices(match));
	return ret;
end

### property handling
function get_dfs(dict::Dict, keytofetch)
	for key in keys(dict)
		if key == keytofetch
		  return Base.get(dict, key, nothing)
		elseif typeof(Base.get(dict, key, nothing)) <: Dict
		  if get_dfs(Base.get(dict, key, nothing), keytofetch) != nothing 
		  	return get_dfs(Base.get(dict, key, nothing), keytofetch) 
		  end
		end
	end
	return nothing
end

function get_properties(match::Match, property_name)
	ret = map(get_dfs, match.correspondence_properties, repeated(property_name));
	if length(ret) != 0
		ret = Array{typeof(ret[1])}(ret);
	end
end

function get_properties(match::Match, fn::Function, args...)
	return fn(match, args...)
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

function is_flagged(match::Match, filter_name = nothing)
  	if filter_name == nothing return match.properties["review"]["flagged"];
	else
	return haskey(match.properties["review"]["flags"], filter_name); end
end

function flag!(match::Match, crit = nothing)
	match.properties["review"]["flagged"] = true;
	if crit != nothing
	match.properties["review"]["flags"][crit[1]] = crit[2:end];
      end
end

"""
Flag a match based on property criteria
"""
function flag!(match::Match, property_name, compare, threshold)
	attributes = get_filtered_properties(match, property_name)
	inds_to_filter = find(i -> compare(i, threshold), attributes)
	if length(inds_to_filter) > 0
		flag!(match, (property_name, compare, threshold))
	end
end

function unflag!(match::Match, filter_name = nothing)
  	if filter_name == nothing
	match.properties["review"]["flagged"] = false;
	for key in keys(match.properties["review"]["flags"])
	  pop!(match.properties["review"]["flags"], key);
	end
      else
	  if haskey(match.properties["review"]["flags"], filter_name)
	  pop!(match.properties["review"]["flags"], filter_name);
	end
	  if length(match.properties["review"]["flags"]) == 0
	  match.properties["review"]["flagged"] = false;
	end
      end
end

function get_correspondence_patches(match::Match, ind)
	src_path = match.src_index;
	dst_path = match.dst_index;
	#assert(src_path[end-1:end] == "h5")

	props = match.correspondence_properties[ind]

		scale = props["ranges"]["scale"];
		highpass_sigma = match.properties["params"]["match"]["highpass_sigma"]
		#src_patch = h5read(src_path, "img", props["ranges"]["src_range"])
		src_pt_loc = props["ranges"]["src_pt_loc"];
		#dst_patch = h5read(dst_path, "img", props["ranges"]["dst_range"])
		dst_pt_loc = props["ranges"]["dst_pt_loc"];

		src_pt = src_pt_loc
		dst_pt = dst_pt_loc

	prepare_patches(src_path, dst_path, props["ranges"]["src_range"], props["ranges"]["dst_range"], props["ranges"]["dst_range_full"], scale, highpass_sigma; from_disk = true)

	function rescale(img)
	  img_new = copy(img)
	  min = minimum(img)
	  for i in 1:length(img_new)
	    img_new[i] = img_new[i] - min
	  end
	  max = maximum(img_new)
	  for i in 1:length(img_new)
	    img_new[i] = img_new[i] / max
	  end
	  return round(UInt8, 255 * img_new)
	end

	src_patch = rescale(SRC_PATCH);
	dst_patch = rescale(DST_PATCH);

	#return SRC_PATCH, DST_PATCH
	xc = normxcorr2_preallocated(SRC_PATCH, DST_PATCH);
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
	  rel_offset = src_offset - dst_offset;
	else
	  rel_offset = src_offset
	end
	dst_pt = src_pt + rel_offset
	dst_pt = ceil(Int64, dst_pt);

	block_range = -block_r:block_r;
	search_range = -(block_r+search_r):(block_r+search_r);

	src_range_full = src_pt[1] + block_range, src_pt[2] + block_range;
	dst_range_full = dst_pt[1] + search_range, dst_pt[2] + search_range;

	range_in_src = intersect(src_range_full[1], 1:src_img_size[1]), intersect(src_range_full[2], 1:src_img_size[2]);
	if length(range_in_src[1]) % 2 != 1
	  range_in_src = range_in_src[1][1:end-1], range_in_src[2];
	end
	if length(range_in_src[2]) % 2 != 1
	  range_in_src = range_in_src[1], range_in_src[2][1:end-1];
	end
	range_in_dst = intersect(dst_range_full[1], 1:dst_img_size[1]), intersect(dst_range_full[2], 1:dst_img_size[2]);

	src_pt_locs = findfirst(range_in_src[1] .== src_pt[1]), findfirst(range_in_src[2] .== src_pt[2]);
	dst_pt_locs = findfirst(range_in_dst[1] .== dst_pt[1]), findfirst(range_in_dst[2] .== dst_pt[2]);
	dst_pt_locs_full = findfirst(dst_range_full[1] .== dst_pt[1]), findfirst(dst_range_full[2] .== dst_pt[2]);

	if src_pt_locs[1] == 0 || src_pt_locs[2] == 0 || dst_pt_locs[1] == 0 || dst_pt_locs[2] == 0 
	return nothing
	end


	return src_index, range_in_src, [src_pt_locs[1], src_pt_locs[2]], dst_index, range_in_dst, dst_range_full, [dst_pt_locs[1], dst_pt_locs[2]], [dst_pt_locs_full[1], dst_pt_locs_full[2]], rel_offset;
end

"""
Template match two images & record translation for source image - already scaled
"""
function prematch(src_index, dst_index, src_image, dst_image, params=get_params(src_index))
  println("prematch:")
	if params["match"]["prematch"] == false return; end
	scale = params["match"]["prematch_scale"];
	ratio = params["match"]["prematch_template_ratio"];
	angles = params["match"]["prematch_angles"];
	highpass_sigma = params["match"]["highpass_sigma"];
	if angles > 1
	angles_to_try = linspace(0, 360, angles + 1)[1:end-1];
	else angles_to_try = [0]	
      	end

	#=
	cur_max_r = 0.0
	cur_rot = 0.0
	cur_offset = 0.0
	=#

	scaled_rads = floor(Int64, ratio * size(src_image, 1) / 2), floor(Int64, ratio * size(src_image, 2) / 2)

	range_in_src = ceil(Int64, size(src_image, 1) / 2) + (-scaled_rads[1]:scaled_rads[1]), ceil(Int64, size(src_image, 2) / 2) + (-scaled_rads[2]:scaled_rads[2])
	dst_range_full = ceil(Int64, size(src_image, 1) / 2) + (-3 * scaled_rads[1]:3 * scaled_rads[1]), ceil(Int64, size(src_image, 2) / 2) + (-3 * scaled_rads[2]: 3 * scaled_rads[2])
	range_in_dst = intersect(dst_range_full[1], 1:size(dst_image, 1)), intersect(dst_range_full[2], 1:size(dst_image, 2));


	# median of the patch from above
	#src_pt_locs = [round(Int64, length(range_in_src[1]) / 2), round(Int64, length(range_in_src[2]) / 2)]
	#src_pt_global = range_in_src[1][src_pt_locs[1]], range_in_src[2][src_pt_locs[2]]
	#src_pt_locs = (1,1)
	#src_pt_global = range_in_src[1][src_pt_locs[1]], range_in_src[2][src_pt_locs[2]]
	# the location of the median in global coordinates
	#src_pt_global = range_in_src[1][src_pt_locs[1]], range_in_src[2][src_pt_locs[2]]

	#src_image_patch = src_image[range_in_src...]


	function try_angle(angle, src_image, dst_image, range_in_src, range_in_dst, dst_range_full)
		print("try angle: $angle deg...")
		src_image_patch = src_image[range_in_src...]
		src_image_rotated = imrotate(src_image_patch, angle; parallel = false)
		range_in_src_rotated = 1:size(src_image_rotated, 1), 1:size(src_image_rotated, 2)

		tform = make_rotation_matrix(angle, size(src_image))
		tform_patch = make_rotation_matrix(angle, size(src_image_patch))

		# set the top left of the original patch to be the location of the point; the global location is just the first point in the src
		src_pt_locs = [1,1]
		src_pt_g = first(range_in_src[1]), first(range_in_src[2])

		# location of the median from above, in patch coordinate
		src_pt_locs_rotated = floor(Int64, ([src_pt_locs..., 1]' * tform_patch)[1:2])
		dst_pt = floor(Int64, ([src_pt_g..., 1]' * tform)[1:2])

		# dst_pt_locs is unused in this, so may be set to zeros; dst_pt_locs_full may actually lie outside, so we use direct calculation and not findfirst
		dst_pt_locs = [0,0];
		dst_pt_locs_full = [dst_pt[1] - dst_range_full[1][1] + 1, dst_pt[2] - dst_range_full[2][1] + 1];

		  #println("src_pt_locs: $src_pt_locs_rotated, dst_pt: $dst_pt, dst_pt_locs: $dst_pt_locs, dst_pt_locs_full: $dst_pt_locs_full")
		rel_offset = [0,0]

		ranges = src_index, range_in_src_rotated, src_pt_locs_rotated, dst_index, range_in_dst, dst_range_full, dst_pt_locs, dst_pt_locs_full, rel_offset
		#if angle == 90
		#	if angle == 0	ImageView.view(dst_image / 255) end
		#if angle == 0 || angle == 180 || angle == 300
		# ImageView.view(src_image_rotated[range_in_src...] / 255)
		#end=#

		match = get_match(src_pt_locs, ranges, src_image_rotated, dst_image, scale, highpass_sigma; full = true)

		if match == nothing return 0, 0, [0, 0] end

		r = match[3]["xcorr"]["r_max"]
		dv = match[3]["vects"]["dv"]
		offset = round(Int64, dv) #/ scale

		println("trying $angle degrees... r: $r, dv: $offset")
		#    println("trying $angle degrees... r: $r, dv: $dv")
		#	if params["registry"]["global_offsets"]
		#	offset = get_offset(dst_index) + round(Int64, dv) #/ scale

		return r, angle, offset
	end
	prematches = pmap(try_angle, angles_to_try, repeated(src_image), repeated(dst_image), repeated(range_in_src), repeated(range_in_dst), repeated(dst_range_full))
	prematches = prematches[prematches .!= nothing]
	cur_max_r = 0.0
	cur_rot = 0.0
	cur_offset = 0.0
	for prematch in prematches
		if prematch[1] > cur_max_r
		  	cur_max_r = prematch[1]
			cur_rot = prematch[2]
			cur_offset = prematch[3]
		end
	end
	update_registry(src_index; rotation = cur_rot, offset = cur_offset);
	print("    ")
	println("Prematch complete... rotation: $cur_rot, offset: $cur_offset")
end

function zeropad_to_meanpad!(img)
  running_sum = 0;
  count = 0;
  zero_entries = Array{Int64, 1}();
 	@simd for i in 1:length(img)
	  @inbounds px = img[i];
	  if px != 0
	    @fastmath running_sum += px;
	    @fastmath count += 1;
	  end
	  if px == 0
	    push!(zero_entries, i)
	  end
	end
	if length(zero_entries) == 0 return end
 	@fastmath nonzero_mean = running_sum / count
	if eltype(img) != Float64
 		@fastmath nonzero_mean = round(eltype(img), nonzero_mean)
	end
 	for i in zero_entries
	  @inbounds px = img[i];
	  if px == 0
	    @fastmath @inbounds img[i] = nonzero_mean;
#=	    if i%2 != 0
	    @fastmath @inbounds img[i] = nonzero_mean;
	  else
	    @fastmath @inbounds img[i] = nonzero_mean + eps::Float64;
	    end =#
	  end
	end
end

# if from_disk src_image / dst_image are indices
function prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, scale, highpass_sigma; from_disk = false)

  	# the following two if statements should only ever be called together
	if size(SRC_PATCH_FULL) != map(length,src_range)
		global SRC_PATCH_FULL = Array{Float64, 2}(map(length, src_range)...);
	end
	if size(DST_PATCH_FULL) != map(length,dst_range_full)
		global DST_PATCH_FULL = Array{Float64, 2}(map(length, dst_range_full)...);
	      else
		DST_PATCH_FULL[:] = 0;
	end

	indices_within_range = findin(dst_range_full[1], dst_range[1]), findin(dst_range_full[2], dst_range[2])
	if !from_disk
	@inbounds SRC_PATCH_FULL[:] = src_image[src_range...]
	@inbounds DST_PATCH_FULL[indices_within_range...] = dst_image[dst_range...]
      else
	if get_rotation(src_image) != 0 @inbounds SRC_PATCH_FULL[:] = imrotate(h5read(get_path(src_image), "img"), get_rotation(src_image); parallel = true)[src_range...];
	else
	@inbounds SRC_PATCH_FULL[:] = h5read(get_path(src_image), "img", src_range);
        end
	@inbounds DST_PATCH_FULL[indices_within_range...] = h5read(get_path(dst_image), "img", dst_range)
      end

#	if isnan(mean(DST_PATCH_FULL)) || isnan(mean(SRC_PATCH_FULL)) return nothing; end

#=	SRC_PATCH_FULL[SRC_PATCH_FULL .== 0] = mean(SRC_PATCH_FULL[SRC_PATCH_FULL .!= 0])
	DST_PATCH_FULL[DST_PATCH_FULL .== 0] = mean(DST_PATCH_FULL[DST_PATCH_FULL .!= 0])=#

	function imscale_src_patch!(img, scale_factor)
  	if scale == 1.0
  		if size(SRC_PATCH) != size(img)
  			bb = ImageRegistration.BoundingBox{Int64}(0,0, size(img, 1), size(img, 2))
   			global SRC_PATCH = zeros(Float64, bb.h, bb.w)
   			global SRC_PATCH_G = zeros(Float64, bb.h, bb.w)
   		end
		@inbounds SRC_PATCH[:] = img
	else
  		tform = make_scale_matrix(scale_factor)
  		bb = ImageRegistration.BoundingBox{Float64}(0,0, size(img, 1), size(img, 2))
  		wbb = tform_bb(bb, tform)
  		tbb = snap_bb(wbb)
  		if size(SRC_PATCH) != (tbb.h, tbb.w)
   			global SRC_PATCH = zeros(Float64, tbb.h, tbb.w)
   			global SRC_PATCH_G = zeros(Float64, tbb.h, tbb.w)
   		end
  		ImageRegistration.imwarp!(SRC_PATCH, img, tform);
	end
	zeropad_to_meanpad!(SRC_PATCH);
  	@inbounds SRC_PATCH_G[:] = SRC_PATCH
	end

	function imscale_dst_patch!(img, scale_factor)
		if scale == 1.0
			if size(DST_PATCH) != size(img)
				bb = ImageRegistration.BoundingBox{Int64}(0,0, size(img, 1), size(img, 2))
				global DST_PATCH = zeros(Float64, bb.h, bb.w)
				global DST_PATCH_G = zeros(Float64, bb.h, bb.w)
			end
			@inbounds DST_PATCH[:] = img
		else
			tform = [scale_factor 0 0; 0 scale_factor 0; 0 0 1];
			bb = ImageRegistration.BoundingBox{Float64}(0,0, size(img, 1), size(img, 2))
			wbb = tform_bb(bb, tform)
			tbb = snap_bb(wbb)
			if size(DST_PATCH) != (tbb.h, tbb.w)
				global DST_PATCH = zeros(Float64, tbb.h, tbb.w)
				global DST_PATCH_G = zeros(Float64, tbb.h, tbb.w)
			end
			ImageRegistration.imwarp!(DST_PATCH, img, tform);
		end
		zeropad_to_meanpad!(DST_PATCH);
		@inbounds DST_PATCH_G[:] = DST_PATCH
	end

	imscale_src_patch!(SRC_PATCH_FULL, scale);	
	imscale_dst_patch!(DST_PATCH_FULL, scale);	
	
	if highpass_sigma != 0
		highpass_sigma = highpass_sigma * scale
		@fastmath Images.imfilter_gaussian_no_nans!(SRC_PATCH_G, [highpass_sigma, highpass_sigma])
		elwise_sub!(SRC_PATCH, SRC_PATCH_G);
		@fastmath Images.imfilter_gaussian_no_nans!(DST_PATCH_G, [highpass_sigma, highpass_sigma])
		elwise_sub!(DST_PATCH, DST_PATCH_G);
      end

      return SRC_PATCH, DST_PATCH
end

function get_match(pt, ranges, src_image, dst_image, params)
	return get_match(pt, ranges, src_image, dst_image, params["match"]["blockmatch_scale"], params["match"]["highpass_sigma"]);
end
"""
Template match two image patches to produce point pair correspondence
"""
function get_match(pt, ranges, src_image, dst_image, scale = 1.0, highpass_sigma = 0; full = false)
	src_index, src_range, src_pt_loc, dst_index, dst_range, dst_range_full, dst_pt_loc, dst_pt_loc_full, rel_offset = ranges;
#=	if sum(src_image[src_range[1], first(src_range[2])]) == 0 && sum(src_image[src_range[1], last(src_range[2])]) == 0 &&
			sum(src_image[first(src_range[1]), src_range[2]]) == 0 && sum(src_image[last(src_range[1]), src_range[2]]) == 0 return nothing end
	if sum(dst_image[dst_range[1], first(dst_range[2])]) == 0 && sum(dst_image[dst_range[1], last(dst_range[2])]) == 0 &&
			sum(dst_image[first(dst_range[1]), dst_range[2]]) == 0 && sum(dst_image[last(dst_range[1]), dst_range[2]]) == 0 return nothing end=#

	#see if any of the edges in the template are fully padded
#=	if sum(src_image[src_range[1], first(src_range[2])]) == 0 return nothing end
	if sum(src_image[src_range[1], last(src_range[2])]) == 0 return nothing end
	if sum(src_image[first(src_range[1]), src_range[2]]) == 0 return nothing end
	if sum(src_image[last(src_range[1]), src_range[2]]) == 0 return nothing end=#
	if sum(src_image[src_range[1], round(Int64,median(src_range[2]))]) == 0 return nothing end
	if sum(src_image[round(Int64,median(src_range[1])), src_range[2]]) == 0 return nothing end
	dst_quart_range_i = linspace(dst_range[1][1], dst_range[1][end], 5)
	dst_quart_range_j = linspace(dst_range[2][1], dst_range[2][end], 5)
	cent_sum = 0;
	cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]):round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[2])]);
	cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]):round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[4])]);
	cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]), round(Int64, dst_quart_range_j[2]):round(Int64, dst_quart_range_j[4])]);
	cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[2]):round(Int64, dst_quart_range_j[4])]);
	if cent_sum == 0 return nothing end
#	if sum(dst_image[dst_range[1], round(Int64,linspace(dst_range[2][1], dst_range[2][end], 5)[2])]) == 0 return nothing end
#	if sum(dst_image[dst_range[1], round(Int64,linspace(dst_range[2][1], dst_range[2][end], 5)[4])]) == 0 return nothing end
#	if sum(dst_image[round(Int64,linspace(dst_range[1][1], dst_range[1][end], 5)[2]), dst_range[2]]) == 0 return nothing end
#	if sum(dst_image[round(Int64,linspace(dst_range[1][1], dst_range[1][end], 5)[4]), dst_range[2]]) == 0 return nothing end

	correspondence_properties = Dict{Any, Any}();
	correspondence_properties["ranges"] = Dict{Any, Any}();
	correspondence_properties["ranges"]["src_pt_loc"] = src_pt_loc;
	correspondence_properties["ranges"]["src_range"] = src_range;
	correspondence_properties["ranges"]["dst_pt_loc"] = dst_pt_loc;
	correspondence_properties["ranges"]["dst_range"] = dst_range;
	correspondence_properties["ranges"]["dst_range_full"] = dst_range_full;
	correspondence_properties["ranges"]["scale"] = scale;

#=
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
		xc = normxcorr2_preallocated(src_image[src_range[1], src_range[2]], padded_img);
	else
	xc = normxcorr2_preallocated(src_image[src_range[1], src_range[2]], dst_image[dst_range[1], dst_range[2]]);
	end
	=#
 if (prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, scale, highpass_sigma) == nothing) return nothing end;
#	prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, scale, highpass_sigma)
	xc = normxcorr2_preallocated(SRC_PATCH, DST_PATCH; shape = full ? "full" : "valid");
#=
	if dst_range != dst_range_full
		indices_within_range = findin(dst_range_full[1], dst_range[1]), findin(dst_range_full[2], dst_range[2])
		intersect_img = dst_image[dst_range...];
		avg = mean(intersect_img);
		if isnan(avg) return nothing; end
		avg = round(eltype(dst_image), avg);
		padded_img = fill(avg, length(dst_range_full[1]), length(dst_range_full[2]));
		padded_img[indices_within_range...] = intersect_img;
		if scale == 1.0 
		xc = normxcorr2_preallocated(src_image[src_range[1], src_range[2]], padded_img);
	      else
		xc = normxcorr2_preallocated(imscale_register_a!(src_image[src_range[1], src_range[2]], scale)[1], imscale_register_b!(padded_img, scale)[1])
	      end
	else
		if scale == 1.0 
		xc = normxcorr2_preallocated(src_image[src_range[1], src_range[2]], dst_image[dst_range[1], dst_range[2]]);
	      else
		xc = normxcorr2_preallocated(imscale_register_a!(src_image[src_range[1], src_range[2]], scale)[1], imscale_register_b!(dst_image[dst_range[1], dst_range[2]], scale)[1])
	      end
	end=#

	r_max = maximum(xc)
	if isnan(r_max) return nothing end;
#	if r_max > 1.0 println("rounding error") end
  	ind = findfirst(r_max .== xc)
	i_max, j_max = ind2sub(xc, ind)
  	if i_max == 0 
    		i_max = size(xc, 1)
  	end

	#=if i_max != 1 && j_max != 1 && i_max != size(xc, 1) && j_max != size(xc, 2)
		xc_w = sum(xc[i_max-1:i_max+1, j_max-1:j_max+1])
		xc_i = sum(xc[i_max+1, j_max-1:j_max+1]) - sum(xc[i_max-1, j_max-1:j_max+1])
		xc_j = sum(xc[i_max-1:i_max+1, j_max+1]) - sum(xc[i_max-1:i_max+1, j_max-1])
		i_max += xc_i / xc_w 
		j_max += xc_j / xc_w 
		if isnan(i_max) || isnan(j_max) 
		  println(xc_i)
		  println(xc_j)
		  println(xc_w)
		  return nothing
		end
	end=#
	if 2 < i_max < size(xc, 1) - 1 && 2 < j_max < size(xc, 2) - 1
		xc_w = sum(xc[i_max-2:i_max+2, j_max-2:j_max+2])
		xc_i = sum(xc[i_max+2, j_max-2:j_max+2]) * 2 + sum(xc[i_max+1, j_max-2:j_max+2]) - sum(xc[i_max-1, j_max-2:j_max+2]) - sum(xc[i_max-2, j_max-2:j_max+2]) * 2
		xc_j = sum(xc[i_max-1:i_max+1, j_max+1]) - sum(xc[i_max-1:i_max+1, j_max-1])
		i_max += xc_i / xc_w 
		j_max += xc_j / xc_w 
		if isnan(i_max) || isnan(j_max) 
		  println(xc_i)
		  println(xc_j)
		  println(xc_w)
		  return nothing
		end
	end

	if scale != 1.0
	#length of the dst_range_full is always odd, so only need to care about the oddity / evenness of the source to decide the size of the xc in full
	if full
	  xc_i_len = length(dst_range_full[1]) + length(src_range[1]) - 1
	  xc_j_len = length(dst_range_full[2]) + length(src_range[2]) - 1
	else
	if isodd(length(src_range[1])) 
	  xc_i_len = length(dst_range_full[1]) - length(src_range[1]) + 1
          else
	  xc_i_len = length(dst_range_full[1]) - length(src_range[1]) 
	end
	if isodd(length(src_range[2])) 
	  xc_j_len = length(dst_range_full[2]) - length(src_range[2]) + 1
          else
	  xc_j_len = length(dst_range_full[2]) - length(src_range[2]) 
	end
      end
	xc_i_locs = linspace(1, xc_i_len, size(xc, 1))
	xc_j_locs = linspace(1, xc_j_len, size(xc, 2))
#=	if myid() == 2
	println(xc_i_len)
	println(size(xc, 1))
	println(step(xc_i_locs))
	println(i_max)
      end=#

	i_max = (i_max - 1) * step(xc_i_locs) + 1
	j_max = (j_max - 1) * step(xc_j_locs) + 1
      end

      if full
	di = Float64(i_max + src_pt_loc[1] - dst_pt_loc_full[1] - length(src_range[1]))
	dj = Float64(j_max + src_pt_loc[2] - dst_pt_loc_full[2] - length(src_range[2]))
      else
	di = Float64(i_max - 1 + src_pt_loc[1] - dst_pt_loc_full[1])
	dj = Float64(j_max - 1 + src_pt_loc[2] - dst_pt_loc_full[2])
      end

	#= # version that treats each value as a bin and sends it to the centre
	xc_i_locs = linspace(1, xc_i_len, size(xc, 1) + 1)
	xc_j_locs = linspace(1, xc_j_len, size(xc, 2) + 1)

	i_max = (i_max - 0.5) * step(xc_i_locs) + 1
	j_max = (j_max - 0.5) * step(xc_j_locs) + 1
      end
	di = Float64(i_max - (dst_pt_loc_full[1] - src_pt_loc[1] + 1))
	dj = Float64(j_max - (dst_pt_loc_full[2] - src_pt_loc[2] + 1))
	=#



	#println("src_range: $src_range, dst_range: $dst_range")
	#println("i_max, j_max: $i_max, $j_max, src_pt_loc: $src_pt_loc, dst_pt_loc: $dst_pt_loc, dst_pt_loc_full = $dst_pt_loc_full")
	#println("di: $di, dj: $dj, scale = $scale")


	correspondence_properties["patches"] = Dict{Any, Any}();
	correspondence_properties["patches"]["src_normalized_dyn_range"] = (maximum(src_image[src_range...]) - minimum(src_image[src_range...])) / typemax(eltype(src_image));
	correspondence_properties["patches"]["src_kurtosis"] = kurtosis(src_image[src_range...]);
	correspondence_properties["xcorr"] = Dict{Any, Any}();
	correspondence_properties["xcorr"]["r_max"] = r_max;
	correspondence_properties["xcorr"]["sigmas"] = Dict{Float64, Any}();
	for beta in [0.5, 0.75, 0.95]
	correspondence_properties["xcorr"]["sigmas"][beta] = sigma(xc, beta) / scale;
        end
	correspondence_properties["vects"] = Dict{Any, Any}();
	correspondence_properties["vects"]["dv"] = [di, dj];
	correspondence_properties["vects"]["norm"] = norm([di, dj]);
	correspondence_properties["posts"] = Dict{Any, Any}();

	return vcat(pt + rel_offset + [di, dj], correspondence_properties);
end

function filter!(match::Match, function_name, compare, threshold, vars...)
	# attributes = get_properties(match, property_name)
	attributes = eval(function_name)(match, vars...)
	if attributes == nothing return 0; end
	inds_to_filter = find(i -> compare(i, threshold), attributes);
	type_name = function_name
	if length(vars) > 0
		type_name = vars[1]
	end	
	push!(match.filters, Dict{Any, Any}(
				"author" => author(),
				"type"	  => type_name,
				"threshold" => threshold,
				"rejected"  => inds_to_filter,
				"function" => function_name
			      ));
	#println("$(length(inds_to_filter)) / $(count_correspondences(match)) rejected.");
	return length(inds_to_filter);
end

function filter!(match::Match, filter)
	return filter!(match, filter...)
end

function get_residual_norms_post(match, ms)
	src_pts_after, dst_pts_after, filtered = get_globalized_correspondences_post(ms, findfirst(match_in_ms -> match_in_ms.src_index == match.src_index && match_in_ms.dst_index == match.dst_index, ms.matches));
	return(map(norm, dst_pts_after - src_pts_after))
end

function check(match::Match, function_name, compare, threshold, vars...)
     return compare(eval(function_name)(match, vars...), threshold)
end

function check!(match::Match, crits) 
	unflag!(match);
	for crit in crits
		if check(match, crit...) 
			flag!(match, crit);
		end
	end
	return is_flagged(match);
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

function Match(src_mesh::Mesh, dst_mesh::Mesh, params=get_params(src_mesh); rotate=0)
	println("Matching $(get_index(src_mesh)) -> $(get_index(dst_mesh)):")
	# if src_mesh == dst_mesh
	# 	println("nothing at")
	# 	println(get_index(src_mesh))
	# 	println(get_index(dst_mesh))
	# 	return nothing
	# end

  	src_index = get_index(src_mesh); dst_index = get_index(dst_mesh);
	src_img = get_image(src_index);  dst_img = get_image(dst_index);

	if params["match"]["prematch"]
	prematch(src_index, dst_index, src_img, dst_img, params);
        end

	src_offset = get_offset(src_index); dst_offset = get_offset(dst_index);
	src_size = get_image_size(src_index); dst_size = get_image_size(dst_index);

	if get_rotation(src_index) != 0
	  src_img = imrotate(src_img, get_rotation(src_index); parallel = true);
	  src_size = get_image_size(src_index; rotated = true);
	  remesh!(src_mesh);
	end

	if params["registry"]["global_offsets"] && get_rotation(dst_index) != 0
	  println("rotation in the dst image detected with global offsets - aborting...")
	  return nothing
	end

	print("computing ranges:")
	@time ranges = map(get_ranges, src_mesh.src_nodes, repeated(src_index), repeated(src_offset), repeated(src_size), repeated(dst_index), repeated(dst_offset), repeated(dst_size), repeated(params["match"]["block_r"]), repeated(params["match"]["search_r"]), repeated(params["registry"]["global_offsets"]));
	ranged_inds = find(i -> i != nothing, ranges);
	ranges = ranges[ranged_inds];
	print("    ")

	println("$(length(ranged_inds)) / $(length(src_mesh.src_nodes)) nodes to check.")

	if length(ranged_inds) != 0
		ranges = Array{typeof(ranges[1]), 1}(ranges);
	end


	print("computing matches:")
	print("    ")
        @time dst_allpoints = pmap(get_match, src_mesh.src_nodes[ranged_inds], ranges, repeated(src_img), repeated(dst_img), repeated(params["match"]["blockmatch_scale"]), repeated(params["match"]["highpass_sigma"])) 
        #@time dst_allpoints = map(get_match, src_mesh.src_nodes[ranged_inds], ranges, repeated(src_img), repeated(dst_img), repeated(params["match"]["blockmatch_scale"]), repeated(params["match"]["highpass_sigma"])) 

	matched_inds = find(i -> i != nothing, dst_allpoints);

	#=for ap in dst_allpoints
	  if typeof(ap) == RemoteException
	    println(ap)
	  end
	end=#

	src_points = copy(src_mesh.src_nodes[ranged_inds][matched_inds]);
	dst_points = [convert(Point, dst_allpoints[ind][1:2]) for ind in matched_inds]
	correspondence_properties = [dst_allpoints[ind][3] for ind in matched_inds]

	filters = Array{Dict{Any, Any}}(0);
  	properties = Dict{Any, Any}(
		"params" => params,
		"review" => Dict{Any, Any}(
				"flagged" => false,
				"flags" => Dict{Any, Any}(),
				"author" => null_author()
				) 
			);

@everywhere gc();

	return Match(src_index, dst_index, src_points, dst_points, correspondence_properties, filters, properties);
end