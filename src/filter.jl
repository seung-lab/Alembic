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

