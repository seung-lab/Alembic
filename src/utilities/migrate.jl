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
  if !haskey(meshset.properties[:params][:solve], :max_iters) 
  	println("MIGRATION: 2016-03-21 Params: added MAX_ITERS in params"); 
        meshset.properties[:params][:solve][:max_iters] = 500;
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

  if !haskey(meshset.matches[match_ind].properties, :review) || !haskey(meshset.matches[match_ind].properties[:review], "author") || !haskey(meshset.matches[match_ind].properties[:review], "flags")
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

  if !haskey(meshset.meshes[mesh_ind].properties, :params) || !haskey(meshset.matches[match_ind].properties, :params)
  	println("MIGRATION: 2016-03-30 MeshSet: adding params to matches / meshes"); 
	for mesh in meshset.meshes
	  mesh.properties[:params] = meshset.properties[:params]
	end
	for match in meshset.matches
	  match.properties[:params] = meshset.properties[:params]
	end
  end

  if length(meshset.properties[:params][:review]) == 0
  	println("MIGRATION: 2016-03-30 MeshSet: adding review criteria"); 
	meshset.properties[:params][:review][:filtered_ratio] = (:get_ratio_filtered, <, 0.2, 20);
	meshset.properties[:params][:review][:ratio_edge_proximity] = (:get_ratio_edge_proximity, >, 0.95)
  end

  if typeof(meshset.matches[match_ind].properties[:review][:flags]) == DataType;
  	println("MIGRATION: 2016-04-11 MeshSet: fixing flags from ::Dict{Any, Any} to Dict{Any, Any}()"); 
	for match in meshset.matches
	  match.properties[:review][:flags] = Dict{Any, Any}();
	end
  end

  if !haskey(meshset.matches[match_ind].correspondence_properties, "posts")
  	println("MIGRATION: 2016-05-11 MeshSet: adding posts to correspondence_properties / moving dv_post and norm_post"); 
	for match in meshset.matches
	  for prop in match.correspondence_properties
	    prop["posts"] = Dict{Any, Any}();
  		if haskey(prop["vects"], "dv_post")
		  prop["posts"]["dv_post"] = prop["vects"]["dv_post"];
		  prop["posts"]["norm_post"] = prop["vects"]["norm_post"];
		  delete!(prop["vects"], "dv_post")
		  delete!(prop["vects"], "norm_post")
		end
	  end
	end
  end

  if !haskey(meshset.properties[:params][:match], :reflexive)
  	println("MIGRATION: 2016-05-19 MeshSet, Mesh, Match: adding reflexive as a field under params"); 
	for mesh in meshset.meshes
	  mesh.properties[:params][:match][:reflexive] = true;
	end
	for match in meshset.matches
	  match.properties[:params][:match][:reflexive] = true;
	end
	  meshset.properties[:params][:match][:reflexive] = true;
  end

  return meshset
end

function migrate_registries!(tile_height=3840)
    premontaged_reg_fn = get_registry_path(premontaged(0,0))
    montaged_reg_fn = get_registry_path(montaged(0,0))
    prealigned_reg_fn = get_registry_path(prealigned(0,0))
    aligned_reg_fn = get_registry_path(aligned(0,0))

    if isfile(premontaged_reg_fn)
      	premontaged_reg = readdlm(premontaged_reg_fn);

    	if size(premontaged_reg, 2) == 4
    	  println("MIGRATION: 2016-04-08: adding sizes to piriform premontages");
    	  premontaged_reg = hcat(premontaged_reg, fill(tile_height, size(premontaged_reg, 1), 2));
    	  writedlm(get_registry_path(premontaged(0,0)), premontaged_reg)
    	end

    	global REGISTRY_PREMONTAGED = parse_registry(get_registry_path(premontaged(0,0)));

    	if size(REGISTRY_PREMONTAGED, 2) == 6
    	  println("MIGRATION: 2016-04-09: adding needs_render to premontaged registry");
    	  premontaged_reg = hcat(premontaged_reg, fill(false, size(premontaged_reg, 1)));
    	  writedlm(get_registry_path(premontaged(0,0)), premontaged_reg)
    	end

    	global REGISTRY_PREMONTAGED = parse_registry(get_registry_path(premontaged(0,0)));

        if size(REGISTRY_PREMONTAGED, 2) == 7
          println("MIGRATION: 2016-08-04: changing needs_render to is_rendered / added rotation to premontaged registry");
          premontaged_reg = hcat(premontaged_reg[:, 1], fill(0, size(premontaged_reg, 1)), premontaged_reg[:, 2:5], !(Array{Bool, 1}(premontaged_reg[:, 6])))
          writedlm(get_registry_path(premontaged(0,0)), premontaged_reg)
        end
        global REGISTRY_PREMONTAGED = parse_registry(get_registry_path(premontaged(0,0)));

    end
    if isfile(montaged_reg_fn)
        montaged_reg = readdlm(montaged_reg_fn);

    	if size(REGISTRY_MONTAGED, 2) == 6
    	  println("MIGRATION: 2016-04-09: adding needs_render to montaged registry");
    	  montaged_reg = hcat(montaged_reg, fill(false, size(montaged_reg, 1)));
    	  writedlm(get_registry_path(montaged(0,0)), montaged_reg)
    	end

    	global REGISTRY_MONTAGED = parse_registry(get_registry_path(montaged(0,0)));

        if size(REGISTRY_MONTAGED, 2) == 7
          println("MIGRATION: 2016-08-04: changing needs_render to is_rendered / added rotation to montaged registry");
          montaged_reg = hcat(montaged_reg[:, 1], fill(0, size(montaged_reg, 1)), montaged_reg[:, 2:5], !(Array{Bool, 1}(montaged_reg[:, 6])))
          writedlm(get_registry_path(montaged(0,0)), montaged_reg)
        end
        global REGISTRY_MONTAGED = parse_registry(get_registry_path(montaged(0,0)));

    end
    if isfile(prealigned_reg_fn)
        prealigned_reg = readdlm(prealigned_reg_fn);

    	if size(REGISTRY_PREALIGNED, 2) == 6
    	  println("MIGRATION: 2016-04-09: adding needs_render to prealigned registry");
    	  prealigned_reg = hcat(prealigned_reg, fill(false, size(prealigned_reg, 1)));
    	  writedlm(get_registry_path(prealigned(0,0)), prealigned_reg)
    	end

    	global REGISTRY_PREALIGNED = parse_registry(get_registry_path(prealigned(0,0)));


    	if size(REGISTRY_PREALIGNED, 2) == 7
    	  println("MIGRATION: 2016-08-04: changing needs_render to is_rendered / added rotation to prealigned registry");
    	  prealigned_reg = hcat(prealigned_reg[:, 1], fill(0, size(prealigned_reg, 1)), prealigned_reg[:, 2:5], !(Array{Bool, 1}(prealigned_reg[:, 6])))
    	  writedlm(get_registry_path(prealigned(0,0)), prealigned_reg)
    	end
    	global REGISTRY_PREALIGNED = parse_registry(get_registry_path(prealigned(0,0)));

    end
    if isfile(aligned_reg_fn)
        aligned_reg = readdlm(aligned_reg_fn);
        
        if size(REGISTRY_ALIGNED, 2) == 6
          println("MIGRATION: 2016-04-09: adding needs_render to aligned registry");
          aligned_reg = hcat(aligned_reg, fill(false, size(aligned_reg, 1)));
          writedlm(get_registry_path(aligned(0,0)), aligned_reg)
        end

        global REGISTRY_ALIGNED = parse_registry(get_registry_path(aligned(0,0)));

    	if size(REGISTRY_ALIGNED, 2) == 7
    	  println("MIGRATION: 2016-08-04: changing needs_render to is_rendered / added rotation to aligned registry");
    	  aligned_reg = hcat(aligned_reg[:, 1], fill(0, size(aligned_reg, 1)), aligned_reg[:, 2:5], !(Array{Bool, 1}(aligned_reg[:, 6])))
    	  writedlm(get_registry_path(aligned(0,0)), aligned_reg)
    	end
    	global REGISTRY_ALIGNED = parse_registry(get_registry_path(aligned(0,0)));

    end
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
	  	if haskey(props, "full")  	
		  	props["ranges"] = props["full"];
			delete!(props, "full");
			if haskey(props, "scale")
			  	props["ranges"]["scale"] = props["scale"];
				delete!(props, "scale")
	        else
			  	props["ranges"]["scale"] = 1.0;
			end
		else
			props["ranges"] = Dict()
			props["ranges"]["scale"] = 1.0
			props["ranges"]["src_range"] = props["src_range"]
			props["ranges"]["src_pt_loc"] = props["src_pt_loc"]
			props["ranges"]["dst_range"] = props["dst_range"]
			props["ranges"]["dst_pt_loc"] = props["dst_pt_loc"]
			delete!(props, "src_range")
			delete!(props, "src_pt_loc")
			delete!(props, "dst_range")
			delete!(props, "dst_pt_loc")
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
	  	props["vects"][:vects_dv] = props["dv"]
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
		if haskey(props, "img")
		  	props["patches"] = props["img"]
			delete!(props, "img")
		else
			props["patches"] = Dict()
			props["patches"][:patches_src_normalized_dyn_range] = props[:patches_src_normalized_dyn_range]
			props["patches"][:patches_src_kurtosis] = props[:patches_src_kurtosis]
			delete!(props, :patches_src_normalized_dyn_range)
			delete!(props, :patches_src_kurtosis)
		end
	end
	return match
end

function migrate_to_patches_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_patches_dict!, ms.matches)
end

function migrate_to_review_dict!(match::Match)
  	match.properties[:review] = Dict{Any, Any}()
	match.properties[:review][:flagged] = false;
	match.properties[:review][:flags] = Dict{Any, Any}();
	match.properties[:review][:author] = null_author();
	return match
end

function migrate_to_review_dict!(ms::MeshSet)
	ms.matches = map(migrate_to_review_dict!, ms.matches)
end

function check_match_flags(ms::MeshSet)
	for match in ms.matches
		if typeof(match.properties[:review][:flags]) == DataType;
			println("MIGRATION: 2016-04-11 MeshSet: fixing flags from ::Dict{Any, Any} to Dict{Any, Any}()"); 
			match.properties[:review][:flags] = Dict{Any, Any}()
		end
	end
	return ms
end

function check_match_flags!(firstindex::Index, lastindex::Index)
	parent_name = get_name(firstindex, lastindex)
	for i in 1:count_children(parent_name)
		ms = load_split(parent_name, i)
		ms = check_match_flags(ms)
		save(ms)
	end
end
