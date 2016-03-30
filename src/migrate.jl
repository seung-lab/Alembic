function migrate!(meshset)
  match_ind = 0;
  mesh_ind = 1;
  for (ind, match) in enumerate(meshset.matches)
    if count_filtered_correspondences(match) != 0 match_ind = ind end
  end
  # migrate matches with sigmas at .5,.6,.7,.8
  if !haskey(meshset.matches[match_ind].correspondence_properties[1], "sigma_5") && !haskey(meshset.matches[match_ind].correspondence_properties[1], "xcorr")
  	println("MIGRATION: 2016-03-20 Match: computing sigmas for correspondences"); 
	migrate_correspondence_sigmas!(meshset);
  end

  # MAX_ITERS for solving
  if !haskey(meshset.properties["params"]["solve"], "max_iters") 
  	println("MIGRATION: 2016-03-21 Params: added MAX_ITERS in params"); 
        meshset.properties["params"]["solve"]["max_iters"] = 500;
  end

  if !haskey(meshset.matches[match_ind].correspondence_properties[1], "xcorr")
  	println("MIGRATION: 2016-03-22 Match: move sigmas / r_max to xcorr dict"); 
	migrate_to_xcorr_dict!(meshset);
  end

  if !haskey(meshset.matches[match_ind].correspondence_properties[1], "ranges")
  	println("MIGRATION: 2016-03-22 Match: move full / scale to ranges dict"); 
	migrate_to_ranges_dict!(meshset);
  end

  if !haskey(meshset.matches[match_ind].correspondence_properties[1], "vects")
  	println("MIGRATION: 2016-03-22 Match: move dv / norm to vects"); 
	migrate_to_vects_dict!(meshset);
  end

  if !haskey(meshset.matches[match_ind].correspondence_properties[1], "patches")
  	println("MIGRATION: 2016-03-22 Match: move img to patches"); 
	migrate_to_patches_dict!(meshset);
  end

  if !haskey(meshset.matches[match_ind].properties, "review") || !haskey(meshset.matches[match_ind].properties["review"], "author") || !haskey(meshset.matches[match_ind].properties["review"], "flags")
  	println("MIGRATION: 2016-03-22 Match: making review / flags dict"); 
	migrate_to_review_dict!(meshset);
  end


  if !haskey(meshset.properties, "meta") || !haskey(meshset.properties["meta"], "parent")
  	println("MIGRATION: 2016-03-22 MeshSet: adding meta dict to properties"); 
	meshset.properties["meta"] = Dict{Any, Any}(
					"parent" => nothing,
					"split_index" => 0
					)
  end

  if !haskey(meshset.meshes[mesh_ind].properties, "params") || !haskey(meshset.matches[match_ind].properties, "params")
  	println("MIGRATION: 2016-03-30 MeshSet: adding params to matches / meshes"); 
	for mesh in meshset.meshes
	  mesh.properties["params"] = meshset.properties["params"]
	end
	for match in meshset.matches
	  match.properties["params"] = meshset.properties["params"]
	end
  end

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
end

#MIGRATION ONLY
function migrate_correspondence_sigmas!(match::Match)
	for (i, props) in enumerate(match.correspondence_properties)
		xc = get_correspondence_patches(match, i)[5]
		props["sigma_5"] = sigma(xc, 0.5)
		props["sigma_6"] = sigma(xc, 0.6)
		props["sigma_7"] = sigma(xc, 0.7)
		props["sigma_8"] = sigma(xc, 0.8)
	end
	return match
end

#MIGRATION ONLY
function migrate_correspondence_sigmas!(ms::MeshSet)
	ms.matches = pmap(migrate_correspondence_sigmas!, ms.matches)
end

function migrate_to_xcorr_dict!(match::Match)
	for props in match.correspondence_properties
	  	props["xcorr"] = Dict{Any, Any}();
	  	props["xcorr"]["r_max"] = props["r_val"]
	  	props["xcorr"]["sigmas"] = Dict{Any, Any}();
	  	props["xcorr"]["sigmas"][0.5] = props["sigma_5"];
	  	props["xcorr"]["sigmas"][0.6] = props["sigma_6"];
	  	props["xcorr"]["sigmas"][0.7] = props["sigma_7"];
	  	props["xcorr"]["sigmas"][0.8] = props["sigma_8"];

		delete!(props, "r_val");
		delete!(props, "sigma_5");
		delete!(props, "sigma_6");
		delete!(props, "sigma_7");
		delete!(props, "sigma_8");
	end
	return match
end

function migrate_to_xcorr_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_xcorr_dict!, ms.matches)
end

function migrate_to_ranges_dict!(match::Match)
	for props in match.correspondence_properties
	  	props["ranges"] = props["full"];
		delete!(props, "full");
		if haskey(props, "scale")
	  	props["ranges"]["scale"] = props["scale"];
		delete!(props, "scale")
	        else
	  	props["ranges"]["scale"] = 1.0;
	end
      end
	return match
end

function migrate_to_ranges_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_ranges_dict!, ms.matches)
end

function migrate_to_vects_dict!(match::Match)
	for props in match.correspondence_properties
	  	props["vects"] = Dict{Any, Any}();
	  	props["vects"]["dv"] = props["dv"]
	  	props["vects"]["norm"] = props["norm"]
		delete!(props, "dv")
		delete!(props, "norm")
	end
	return match
end

function migrate_to_vects_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_vects_dict!, ms.matches)
end

function migrate_to_patches_dict!(match::Match)
	for props in match.correspondence_properties
	  	props["patches"] = props["img"]
		delete!(props, "img")
	end
	return match
end

function migrate_to_patches_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_patches_dict!, ms.matches)
end

function migrate_to_review_dict!(match::Match)
  	match.properties["review"] = Dict{Any, Any}()
	match.properties["review"]["flagged"] = false;
	match.properties["review"]["flags"] = Dict{Any, Any}();
	match.properties["review"]["author"] = null_author();
	return match
end

function migrate_to_review_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_review_dict!, ms.matches)
end

