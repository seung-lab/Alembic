global PARAMS = Dict()

function load_params(fn::AbstractString)
	params = JSON.parsefile(fn)
	return load_params(params)
end

function load_params(params::Dict)
	dirs = params["dirs"]
	mesh = params["mesh"]
	match = params["match"]
	filter = params["filter"]
	solve = params["solve"]
	render = params["render"]
	root = joinpath(dirs["bucket"], dirs["dataset"])
	global PARAMS = Dict(
		 :name => params["name"],
		 :task_method => params["task_method"],
		 :dirs => Dict(
		 	:dataset => dirs["dataset"],
		 	:bucket => dirs["bucket"],
		 	:src_image => dirs["src_image"],
		 	:dst_image => dirs["dst_image"],
		 	:match_image => dirs["match_image"],
		 	:mask => dirs["mask"] != "" ? dirs["mask"] : "",
		 	:mesh => dirs["mesh"],
		 	:match => dirs["match"],
		 	:meshset => dirs["meshset"],
		 	:roi => dirs["roi"],
		 	:src_patch => dirs["src_patch"],
		 	:dst_patch => dirs["dst_patch"],
		 	:src_match => dirs["src_match"],
		 	:dst_match => dirs["dst_match"],
		 	:xc => dirs["xc"]
			),
	     :mesh => Dict(
	     	:z_start => mesh["z_start"],
	     	:z_stop => mesh["z_stop"],
			:mesh_length => mesh["mesh_length"]
			), 
	     :match => Dict(
			:prematch => false,
			:mip => match["mip"],
			:ignore_value => match["ignore_value"],
			:block_r => match["block_r"], 
			:search_r => match["search_r"],
			:bandpass_sigmas => match["bandpass_sigmas"],
			:depth => match["depth"],
			:symmetric => match["symmetric"]
			),
	     :solve => Dict(
			:mesh_spring_coeff => solve["mesh_spring_coeff"],
			:match_spring_coeff => solve["match_spring_coeff"],
			:ftol_cg => solve["ftol_cg"],
			:max_iters => solve["max_iters"],
	     	:use_cg => solve["use_cg"],
	     	:eta_gd => solve["eta_gd"],
	     	:ftol_gd => solve["ftol_gd"],
	     	:eta_newton => solve["eta_newton"],
	     	:ftol_newton => solve["ftol_newton"],
	     	:method => solve["method"]
			),
	     :filter => Dict(
     		:sigma_filter_high => (1,:get_correspondence_properties, Symbol(>), filter["sigma_filter_high"], Symbol("xcorr_sigma_0.95")),
     		:sigma_filter_mid => (2,:get_correspondence_properties, Symbol(>), filter["sigma_filter_mid"], Symbol("xcorr_sigma_0.75")),
     		:sigma_filter_low => (3,:get_correspondence_properties, Symbol(>), filter["sigma_filter_low"], Symbol("xcorr_sigma_0.5")),
     		:dyn_range_filter => (4,:get_correspondence_properties, Symbol(<), filter["dyn_range_filter"], :patches_src_normalized_dyn_range),
     		:r_filter => (5,:get_correspondence_properties, Symbol(<), filter["r_filter"], :xcorr_r_max),
     		:kurtosis_filter => (6,:get_correspondence_properties, Symbol(>), filter["kurtosis_filter"], :patches_src_kurtosis),
     		:kurtosis_filter_dst => (7,:get_correspondence_properties, Symbol(>), filter["kurtosis_filter_dst"], :patches_dst_kurtosis),
     		:kurtosis_filter_edge => (8,:get_correspondence_properties, Symbol(<), filter["kurtosis_filter_edge"], :patches_src_kurtosis),
     		:r_delta_low => (9,:get_correspondence_properties, Symbol(<), filter["r_delta_low"], Symbol("xcorr_delta_5"))
			),
	     :render => Dict(
     		:mip => render["mip"]
		    ),
	     :review => Dict(
     		:too_few_corresps => (:count_correspondences, Symbol(<), 100),
			:rejected_ratio => (:get_ratio_rejected, Symbol(>), 0.35),
			:ratio_edge_proximity => (:get_ratio_edge_proximity, Symbol(>), 0.80)
			)
		)
end
