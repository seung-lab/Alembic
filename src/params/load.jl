global PARAMS = Dict()
global BUCKET = ""
global DATASET = ""

function load_params(fn)
	params_dict = JSON.parsefile(fn)
	d = params_dict["dirs"]
	p = params_dict["params"]
	global BUCKET = d["bucket"]
	global DATASET = d["dataset"]
	root = joinpath(BUCKET, DATASET)
	global PARAMS = Dict(
		 :dirs => Dict(
		 	:bucket => root,
		 	:src_image => joinpath(root, d["src_image"]),
		 	:dst_image => joinpath(root, d["dst_image"]),
		 	# :match_image => joinpath(root, d["src_image"], d["match_image"]),
		 	# :mask => joinpath(root, d["src_image"], d["mask"]),
		 	# :mesh => joinpath(root, d["src_image"], d["mesh"]),
		 	# :match => joinpath(root, d["src_image"], d["match"])
		 	:match_image => joinpath(root, d["match_image"]),
		 	:mask => joinpath(root, d["mask"]),
		 	:mesh => joinpath(root, d["mesh"]),
		 	:match => joinpath(root, d["match"]),
		 	:cache => d["cache"]
		  ),
	     :mesh => Dict(
			:mesh_length => p["mesh_length"]), 
	     :match => Dict(
	     	:z_start => p["z_start"],
	     	:z_stop => p["z_stop"],
			:prematch => false,
			:mip => p["mip"],
			:block_r => p["block_r"], 
			:search_r => p["search_r"],
			:bandpass_sigmas => p["bandpass_sigmas"],
			:depth => 1,
			:reflexive => false),
	     :solve => Dict(
			:mesh_spring_coeff => p["mesh_spring_coeff"],
			:match_spring_coeff => p["match_spring_coeff"],
			:ftol_cg => p["ftol_cg"],
			:max_iters => p["max_iters"],
	     	:use_conjugate_gradient => p["use_cg"],
	     	:eta_gd => 0,
	     	:ftol_gd => 0,
	     	:eta_newton => 0,
	     	:ftol_newton => 0,
	     	:method => p["method"]),
	     :filter => Dict(
	     		:sigma_filter_high => (1,:get_correspondence_properties, >, 6.5, Symbol("xcorr_sigma_0.95")),
	     		:sigma_filter_mid => (2,:get_correspondence_properties, >, 20, Symbol("xcorr_sigma_0.75")),
	     		:sigma_filter_low => (3,:get_correspondence_properties, >, 40, Symbol("xcorr_sigma_0.5")),
	     		:dyn_range_filter => (4,:get_correspondence_properties, <, 0.50, :patches_src_normalized_dyn_range),
	     		:r_filter => (5,:get_correspondence_properties, <, 0.0275, :xcorr_r_max),
	     		:kurtosis_filter => (6,:get_correspondence_properties, >, 25, :patches_src_kurtosis),
	     		:kurtosis_filter_dst => (7,:get_correspondence_properties, >, 50, "dst_kurtosis"),
	     		:kurtosis_filter_edge => (8,:get_correspondence_properties, <, -1.60, :patches_src_kurtosis),
			      ),
	     :render => Dict(
			      ),
	     :review => Dict(
	     		:too_few_corresps => (:count_correspondences, <, 100),
				:rejected_ratio => (:get_ratio_rejected, >, 0.35),
				:ratio_edge_proximity => (:get_ratio_edge_proximity, >, 0.80)
		)
		)
end