function migrate!(meshset)
  if !haskey(meshset.properties["params"]["solve"], "max_iters") println("MIGRATION: ADDED MAX_ITERS IN PARAMS"); meshset.properties["params"]["solve"]["max_iters"] = 500; end
  return meshset
end

function clean_correspondence_sigmas!(match::Match)
	for i in 1:count_correspondences(match)
		for (k, v) in match.correspondence_properties[i]
			if contains(k, "sigma")
				if isnan(v)
					match.correspondence_properties[i][k] = Inf
				end
			end
		end
	end
	return match
end
