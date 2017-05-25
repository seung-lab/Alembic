### ratios
function get_ratio_filtered(match::Match, min_corresps = 0) 
  if count_correspondences(match) < min_corresps return 1.0 end # ignore low match cases
return count_filtered_correspondences(match) / max(count_correspondences(match), 1); end

function get_ratio_rejected(match::Match, min_corresps = 0) 
  if count_correspondences(match) < min_corresps return 0.0 end # ignore low match cases
return count_rejected_correspondences(match) / max(count_correspondences(match), 1); end

function get_ratio_edge_proximity(match::Match)
     if count_filtered_correspondences(match) == 0 return 0.0 end
     norms = map(norm, get_filtered_properties(match, "dv"))
return maximum(norms) / match.properties["params"]["match"]["search_r"]; end

function count_outlier_norms(match::Match, sigma=3)
	return sum(get_norms_std_sigmas(match) .> sigma)
end

function get_median_dv(match::Match)
	if count_filtered_correspondences(match) == 0 
		return 0.0
	end
	dvs = get_filtered_properties(match, "dv")
	x, y = [dv[1] for dv in dvs], [dv[2] for dv in dvs]
	return [median(x), median(y)]
end

function get_maximum_centered_norm(match::Match)
	if count_filtered_correspondences(match) == 0 
		return 0.0
	end
	dvs = get_filtered_properties(match, "dv")
	x, y = [dv[1] for dv in dvs], [dv[2] for dv in dvs]
	med = [median(x), median(y)]
	norms = map(norm, [dv - med for dv in dvs])
	return maximum(norms)
end

function get_centered_norms(match::Match)
	if count_correspondences(match) == 0 
		return nothing
	end
	dvs = get_properties(match, "dv")
	x, y = [dv[1] for dv in dvs], [dv[2] for dv in dvs]
	med = [median(x), median(y)]
	norms = map(norm, [dv - med for dv in dvs])
	return norms
end

function get_properties_ratios(match::Match, args...)
  	props_den = get_properties(match, args[end])
  	for arg in args[(end-1):-1:1]
	props_num = get_properties(match, arg)
	props_den = props_num ./ props_den
      	end
	return props_den
end

function get_norm_std(match::Match)
	if count_filtered_correspondences(match) == 0 
		return 0.0
	end
	norms = convert(Array{Float64}, map(norm, get_filtered_properties(match, "dv")))
	return std(convert(Array{Float64}, norms))
end

function get_norms_std_sigmas(match::Match)
	if count_filtered_correspondences(match) == 0 
		return 0.0
	end
	norms = convert(Array{Float64}, map(norm, get_properties(match, "dv")))
	filtered_norms = convert(Array{Float64}, map(norm, get_filtered_properties(match, "dv")))
	mu = mean(filtered_norms)
	stdev = std(filtered_norms)
	return (norms - mu) / stdev
end

function get_normalized_norm_from_consensus(m::Match, r)

#std1s = similar(m.src_points, Float64);
#std2s = similar(m.src_points, Float64);

dv1s = zeros(Float64, length(m.src_points));
dv2s = zeros(Float64, length(m.src_points));


normalized_norms = similar(m.src_points, Float64);
r_sq = r * r;

for i in 1:count_correspondences(m)
	accepted = 0
	@inbounds p1 = m.src_points[i][1];
	@inbounds p2 = m.src_points[i][2];

	for k in 1:count_correspondences(m)
		@fastmath @inbounds d1 = m.src_points[k][1] - p1
		@fastmath @inbounds d2 = m.src_points[k][2] - p2
		@fastmath @inbounds d_norm_sq = d1^2 + d2^2
			if d_norm_sq == 0 
				@inbounds dv1s[k] = 0;
				@inbounds dv2s[k] = 0;
			continue end;
		@fastmath if d_norm_sq < r_sq
			if d_norm_sq == 0 continue end;
			@fastmath @inbounds dv1s[k] = m.dst_points[k][1] - m.src_points[k][1];
			@fastmath @inbounds dv2s[k] = m.dst_points[k][2] - m.src_points[k][2];
			accepted += 1;
		else
			@inbounds dv1s[k] = 0;
			@inbounds dv2s[k] = 0;
		end
	      
	end

	if accepted == 0 normalized_norms[i] = Inf; continue; end

	@fastmath m1 = sum(dv1s)/accepted;
	@fastmath m2 = sum(dv2s)/accepted;

	for k in 1:count_correspondences(m)
	  	if dv1s[k] == 0 continue; end
		@fastmath @inbounds dv1s[k] = (dv1s[k] - m1)^2;
	end
	for k in 1:count_correspondences(m)
	  	if dv2s[k] == 0 continue; end
		@fastmath @inbounds dv2s[k] = (dv2s[k] - m2)^2;
	end

	sig1 = sqrt(sum(dv1s) / accepted)
	sig2 = sqrt(sum(dv2s) / accepted)

	@fastmath @inbounds std1 = (m.dst_points[i][1] - p1 - m1) / sig1;
	@fastmath @inbounds std2 = (m.dst_points[i][2] - p2 - m2) / sig2;
	@fastmath @inbounds normalized_norms[i] = sqrt(std1^2 + std2^2);
#=
	@fastmath @inbounds diff_sq1 = (m.dst_points[i][1] - p1 - m1) ^ 2;
	@fastmath @inbounds diff_sq2 = (m.dst_points[i][2] - p2 - m2) ^ 2;
	@fastmath @inbounds norms[i] = sqrt(diff_sq1 + diff_sq2)
	=#
end
	return normalized_norms;
end

function get_normalized_norm_from_filtered_consensus(m::Match, r)

#std1s = similar(m.src_points, Float64);
#std2s = similar(m.src_points, Float64);

dv1s = zeros(Float64, length(m.src_points));
dv2s = zeros(Float64, length(m.src_points));

normalized_norms = similar(m.src_points, Float64);
r_sq = r * r;

filtered_inds = Array{Int64, 1}(get_filtered_indices(m));
filtered_len = length(filtered_inds)

src_points = Points(length(filtered_inds))
dst_points = Points(length(filtered_inds))
src_points[:] = m.src_points[filtered_inds]
dst_points[:] = m.dst_points[filtered_inds]

for i in 1:count_correspondences(m)
	accepted = 0
	@inbounds p1 = m.src_points[i][1];
	@inbounds p2 = m.src_points[i][2];

	for k in 1:filtered_len
		@fastmath @inbounds d1 = src_points[k][1] - p1
		@fastmath @inbounds d2 = src_points[k][2] - p2
		@fastmath @inbounds d_norm_sq = d1^2 + d2^2
		@fastmath if d_norm_sq < r_sq
			if d_norm_sq == 0 
				@inbounds dv1s[k] = 0;
				@inbounds dv2s[k] = 0;
			continue end;
			@fastmath @inbounds dv1s[k] = dst_points[k][1] - src_points[k][1];
			@fastmath @inbounds dv2s[k] = dst_points[k][2] - src_points[k][2];
			accepted += 1;
		else
			@inbounds dv1s[k] = 0;
			@inbounds dv2s[k] = 0;
		end
	      
	end

	if accepted == 0 normalized_norms[i] = Inf; continue; end

	@fastmath m1 = sum(dv1s)/accepted;
	@fastmath m2 = sum(dv2s)/accepted;

	for k in 1:count_correspondences(m)
	  	if dv1s[k] == 0 continue; end
		@fastmath @inbounds dv1s[k] = (dv1s[k] - m1)^2;
	end
	for k in 1:count_correspondences(m)
	  	if dv2s[k] == 0 continue; end
		@fastmath @inbounds dv2s[k] = (dv2s[k] - m2)^2;
	end

	sig1 = sqrt(sum(dv1s) / accepted)
	sig2 = sqrt(sum(dv2s) / accepted)

	@fastmath @inbounds std1 = (m.dst_points[i][1] - p1 - m1) / sig1;
	@fastmath @inbounds std2 = (m.dst_points[i][2] - p2 - m2) / sig2;
	@fastmath @inbounds normalized_norms[i] = sqrt(std1^2 + std2^2);
end
	return normalized_norms;
end

function get_filtered_norm_diff_from_filtered_consensus(ms, r, threshold)
	normdiffs = (map(get_filtered_norm_diff_from_filtered_consensus, ms.matches, repeated(r))...)
	for i in 1:count_matches(ms)
	  if maximum(normdiffs[i]) > threshold
	    println(i)
	  end
	end
	return vcat(normdiffs...)
end

function get_filtered_norm_diff_from_filtered_consensus(m::Match, r)

#std1s = similar(m.src_points, Float64);
#std2s = similar(m.src_points, Float64);


norms_diff = similar(m.src_points, Float64);
r_sq = r * r;

filtered_inds = Array{Int64, 1}(get_filtered_indices(m));
filtered_len = length(filtered_inds)

src_points = Points(length(filtered_inds))
dst_points = Points(length(filtered_inds))
src_points[:] = m.src_points[filtered_inds]
dst_points[:] = m.dst_points[filtered_inds]
dv1s = zeros(Float64, filtered_len);
dv2s = zeros(Float64, filtered_len);

for i in 1:count_correspondences(m)
	accepted = 0
	@inbounds p1 = m.src_points[i][1];
	@inbounds p2 = m.src_points[i][2];

	for k in 1:filtered_len
		@fastmath @inbounds d1 = src_points[k][1] - p1
		@fastmath @inbounds d2 = src_points[k][2] - p2
		@fastmath @inbounds d_norm_sq = d1^2 + d2^2
		@fastmath if d_norm_sq < r_sq
			if d_norm_sq == 0 
				@inbounds dv1s[k] = 0;
				@inbounds dv2s[k] = 0;
			continue end;
			@fastmath @inbounds dv1s[k] = dst_points[k][1] - src_points[k][1];
			@fastmath @inbounds dv2s[k] = dst_points[k][2] - src_points[k][2];
			accepted += 1;
			#=if i == 2
			  println("dv1: $(dv1s[k]), dv2: $(dv2s[k])")
			end =#
		else
			@inbounds dv1s[k] = 0;
			@inbounds dv2s[k] = 0;
		end
	      
	end

	if accepted == 0 norms_diff[i] = Inf; continue; end


	@fastmath m1 = sum(dv1s)/accepted;
	@fastmath m2 = sum(dv2s)/accepted;


	@fastmath @inbounds d1 = (m.dst_points[i][1] - p1 - m1);
	@fastmath @inbounds d2 = (m.dst_points[i][2] - p2 - m2);
	@fastmath @inbounds norms_diff[i] = sqrt(d1^2 + d2^2);
end
	return norms_diff[filtered_inds];
end

function get_filtered_normalized_norm_from_filtered_consensus(m::Match, r)

#std1s = similar(m.src_points, Float64);
#std2s = similar(m.src_points, Float64);


normalized_norms = similar(m.src_points, Float64);
r_sq = r * r;

filtered_inds = Array{Int64, 1}(get_filtered_indices(m));
filtered_len = length(filtered_inds)

src_points = Points(length(filtered_inds))
dst_points = Points(length(filtered_inds))
src_points[:] = m.src_points[filtered_inds]
dst_points[:] = m.dst_points[filtered_inds]
dv1s = zeros(Float64, filtered_len);
dv2s = zeros(Float64, filtered_len);

for i in 1:count_correspondences(m)
	accepted = 0
	@inbounds p1 = m.src_points[i][1];
	@inbounds p2 = m.src_points[i][2];

	for k in 1:filtered_len
		@fastmath @inbounds d1 = src_points[k][1] - p1
		@fastmath @inbounds d2 = src_points[k][2] - p2
		@fastmath @inbounds d_norm_sq = d1^2 + d2^2
		@fastmath if d_norm_sq < r_sq
			if d_norm_sq == 0 
				@inbounds dv1s[k] = 0;
				@inbounds dv2s[k] = 0;
			continue end;
			@fastmath @inbounds dv1s[k] = dst_points[k][1] - src_points[k][1];
			@fastmath @inbounds dv2s[k] = dst_points[k][2] - src_points[k][2];
			accepted += 1;
		else
			@inbounds dv1s[k] = 0;
			@inbounds dv2s[k] = 0;
		end
	      
	end

	#if accepted == 0 normalized_norms[i] = Inf; continue; end
	if accepted == 0 normalized_norms[i] = Inf; continue; end

	@fastmath m1 = sum(dv1s)/accepted;
	@fastmath m2 = sum(dv2s)/accepted;

	for k in 1:filtered_len
	  	if dv1s[k] == 0 continue; end
		@fastmath @inbounds dv1s[k] = (dv1s[k] - m1)^2;
	end
	for k in 1:filtered_len
	  	if dv2s[k] == 0 continue; end
		@fastmath @inbounds dv2s[k] = (dv2s[k] - m2)^2;
	end

	sig1 = sqrt(sum(dv1s) / accepted)
	sig2 = sqrt(sum(dv2s) / accepted)

	@fastmath @inbounds std1 = (m.dst_points[i][1] - p1 - m1) / sig1;
	@fastmath @inbounds std2 = (m.dst_points[i][2] - p2 - m2) / sig2;
	@fastmath @inbounds normalized_norms[i] = sqrt(std1^2 + std2^2);
end
	return normalized_norms[filtered_inds];
end
