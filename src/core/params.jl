global PARAMS = Dict()

function load_params(fn::AbstractString)
	params = JSON.parsefile(fn)
	return load_params(params)
end

function load_params(params::Dict)
	data = params["data"]
	mesh = params["mesh"]
	match = params["match"]
	filter = params["filter"]
	solve = params["solve"]
	global PARAMS = Dict(
		 :name => params["name"],
		 :task_method => params["task_method"],
		 :src_image => Dict(
		 	:path => data["src_image"]["path"],
		 	:mip => data["src_image"]["mip"]
		 	),
		 :dst_image => Dict(
		 	:path => data["dst_image"]["path"],
		 	:mip => data["dst_image"]["mip"],
		 	:interpolation => data["dst_image"]["interpolation"]
		 	),
		 :match_image => Dict(
		 	:path => data["match_image"]["path"],
		 	:mip => data["match_image"]["mip"],
		 	:value => data["match_image"]["value"]
		 	),
		 :defect_mask => Dict(
		 	:path => data["defect_mask"]["path"],
		 	:mip => data["defect_mask"]["mip"],
		 	:value => data["defect_mask"]["value"],
		 	:apply => data["defect_mask"]["apply"]
		 	),
		 :defect_split => Dict(
		 	:path => data["defect_split"]["path"],
		 	:mip => data["defect_split"]["mip"],
		 	:value => data["defect_split"]["value"],
		 	:apply => data["defect_split"]["apply"]
		 	),
		 :roi_mask => Dict(
		 	:path => data["roi_mask"]["path"],
		 	:mip => data["roi_mask"]["mip"],
		 	:value => data["roi_mask"]["value"],
		 	:apply => data["roi_mask"]["apply"]
		 	),
		 :src_match => Dict(
		 	:path => data["src_match"]["path"],
		 	:mip => data["src_match"]["mip"]
		 	),
		 :dst_match => Dict(
		 	:path => data["dst_match"]["path"],
		 	:mip => data["dst_match"]["mip"]
		 	),
		 :src_patch => Dict(
		 	:path => data["src_patch"]["path"],
		 	:mip => data["src_patch"]["mip"]
		 	),
		 :dst_patch => Dict(
		 	:path => data["dst_patch"]["path"],
		 	:mip => data["dst_patch"]["mip"]
		 	),
		 :xc => Dict(
		 	:path => data["xc"]["path"],
		 	:mip => data["xc"]["mip"]
		 	),
		 :mesh => Dict(
		 	:path => mesh["path"],
	     	:z_start => mesh["z_start"],
	     	:z_stop => mesh["z_stop"],
			:mesh_length => mesh["mesh_length"]
		 	),
		 :match => Dict(
		 	:path => match["path"],
			:prematch => false,
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
	     :review => Dict(
     		:too_few_corresps => (:count_correspondences, Symbol(<), 100),
			:rejected_ratio => (:get_ratio_rejected, Symbol(>), 0.35),
			:ratio_edge_proximity => (:get_ratio_edge_proximity, Symbol(>), 0.80)
			)
		)
end
