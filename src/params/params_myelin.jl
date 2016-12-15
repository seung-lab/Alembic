MESH_LENGTH_MONTAGE = 200
GLOBAL_OFFSETS_MONTAGE = true
BLOCKMATCH_SCALE_MONTAGE = 1.0
BLOCK_R_MONTAGE = 140
SEARCH_R_MONTAGE = 400
PREMATCH_MONTAGE = false
MESH_SPRING_COEFF_MONTAGE = 1.0
MATCH_SPRING_COEFF_MONTAGE = 100.0 
FTOL_CG_MONTAGE = 1e-8
MAX_ITERS_MONTAGE = 2000
USE_CONJUGATE_GRADIENT_MONTAGE = true
ETA_GD_MONTAGE = 0.6
FTOL_GD_MONTAGE = 1e-8
ETA_NEWTON_MONTAGE = 0.6
FTOL_NEWTON_MONTAGE = 1e-16


MESH_LENGTH_PREALIGNMENT = 4000
GLOBAL_OFFSETS_PREALIGNMENT = false
BLOCKMATCH_SCALE_PREALIGNMENT = 0.125
BLOCK_R_PREALIGNMENT = 600
SEARCH_R_PREALIGNMENT = 3000
PREMATCH_PREALIGNMENT = false
PREMATCH_TEMPLATE_RATIO_PREALIGNMENT = 0.15
PREMATCH_SCALE_PREALIGNMENT = 0.10
PREMATCH_ANGLES_PREALIGNMENT = 0;

MESH_LENGTH_ALIGNMENT = 200
GLOBAL_OFFSETS_ALIGNMENT = true
BLOCKMATCH_SCALE_ALIGNMENT = 0.20
BLOCK_R_ALIGNMENT = 800
#SEARCH_R_ALIGNMENT = 250
#SEARCH_R_ALIGNMENT = 350
SEARCH_R_ALIGNMENT = 1200
#SEARCH_R_ALIGNMENT = 550
PREMATCH_ALIGNMENT = false
MESH_SPRING_COEFF_ALIGNMENT = 1.0
MATCH_SPRING_COEFF_ALIGNMENT = 20.0
FTOL_CG_ALIGNMENT = 1e-6
MAX_ITERS_ALIGNMENT = 2000
USE_CONJUGATE_GRADIENT_ALIGNMENT = true
ETA_GD_ALIGNMENT = 0.01
FTOL_GD_ALIGNMENT = 3e-3
ETA_NEWTON_ALIGNMENT = 0.5
FTOL_NEWTON_ALIGNMENT = 1e-8

global PARAMS_MONTAGE = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_MONTAGE), 
			     "match" => Dict(
					"prematch" => PREMATCH_MONTAGE,
					"blockmatch_scale" => BLOCKMATCH_SCALE_MONTAGE,
					"block_r" => BLOCK_R_MONTAGE, 
					"search_r" => SEARCH_R_MONTAGE,
					"bandpass_sigmas" => (0,20),
					"depth" => 1,
					"reflexive" => false),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_MONTAGE,
					"match_spring_coeff" => MATCH_SPRING_COEFF_MONTAGE,
					"ftol_cg" => FTOL_CG_MONTAGE,
					"max_iters" => MAX_ITERS_MONTAGE,
			     	"use_conjugate_gradient" => USE_CONJUGATE_GRADIENT_MONTAGE,
			     	"eta_gd" => ETA_GD_MONTAGE,
			     	"ftol_gd" => FTOL_GD_MONTAGE,
			     	"eta_newton" => ETA_NEWTON_MONTAGE,
			     	"ftol_newton" => FTOL_NEWTON_MONTAGE),
			     "filter" => Dict(
			     		"sigma_filter_high" => (1,:get_properties, >, 2, 0.95),
			     		"sigma_filter_mid" => (2,:get_properties, >, 120, 0.75),
			     		"sigma_filter_low" => (2,:get_properties, >, 350, 0.50),
			     		"r_filter_min" => (2,:get_properties, <, 0.03, "r_max"),
			     		"r_filter_max" => (2,:get_properties, >, 1, "r_max"),
					"centered_norm_filter" => (2,:get_centered_norms, >, 75)
						# "norm_filter" => (:get_norms_std_sigmas, >, 5)
			     		# "norm_filter" => (:get_norms_std_sigmas, >, 2.5)
					),
			     "render" => Dict(
			     		"crop" => [0, 0],
			     		"thumbnail_scale" => 0.02
					      ),
			     "review" => Dict(
						# "too_few_corresps" => (:count_correspondences, <, 10),
						"rejected_ratio" => (:get_ratio_rejected, >, 0.66, 20),
						"ratio_edge_proximity" => (:get_ratio_edge_proximity, >, 0.95),
						#"norm_outliers" => (:count_outlier_norms, >, 0, 3), # too useless because they're so close to each other to begin with
						"centered_norm" => (:get_maximum_centered_norm, >, 80)
					      ),
			     "registry" => Dict(
					"global_offsets" => GLOBAL_OFFSETS_MONTAGE
					)
			     )


global PARAMS_MONTAGE_FALLBACK = PARAMS_MONTAGE;

#PARAMS_MONTAGE_FALLBACK["match"]["search_r"] = 128;

global PARAMS_PREALIGNMENT = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_PREALIGNMENT), 
			     "match" => Dict(
					"prematch" => PREMATCH_PREALIGNMENT,
					"prematch_template_ratio" => PREMATCH_TEMPLATE_RATIO_PREALIGNMENT,
					"prematch_scale" => PREMATCH_SCALE_PREALIGNMENT, 
					"prematch_angles" => PREMATCH_ANGLES_PREALIGNMENT, 
					"blockmatch_scale" => BLOCKMATCH_SCALE_PREALIGNMENT,
					"block_r" => BLOCK_R_PREALIGNMENT, 
					"search_r" => SEARCH_R_PREALIGNMENT,
					"bandpass_sigmas" => (0, 0),
					"depth" => 1,
					"reflexive" => false),
			     "solve" => Dict(
					"method" => "regularized",
					"lambda" => 0.9),
					# "mesh_spring_coeff" => MESH_SPRING_COEFF_PREALIGNMENT,
					# "match_spring_coeff" => MATCH_SPRING_COEFF_PREALIGNMENT,
					# "ftol_cg" => FTOL_CG_PREALIGNMENT,
					# "max_iters" => MAX_ITERS_PREALIGNMENT,
					# "use_conjugate_gradient" => USE_CONJUGATE_GRADIENT_PREALIGNMENT,
					# "eta_gd" => ETA_GD_PREALIGNMENT,
					# "ftol_gd" => FTOL_GD_PREALIGNMENT,
					# "eta_newton" => ETA_NEWTON_PREALIGNMENTE,
					# "ftol_newton" => FTOL_NEWTON_PREALIGNMENT)
			     "filter" => Dict(
			     		"sigma_filter_high" => (2,:get_properties, >, 5, 0.95),
			     		"sigma_filter_mid" => (2,:get_properties, >, 40, 0.75),
			     		"sigma_filter_low" => (2,:get_properties, >, 150, 0.50)
			     		# "norm_filter" => (:get_norms_std_sigmas, >, 5)
					      ),
			     "render" => Dict(
			     		"thumbnail_scale" => 0.05
					      ),
			     "review" => Dict(
				     	# "r_below" => (:count_filtered_properties, >, 0, "r_max", <, 0.2),
			     		"too_few_corresps" => (:count_correspondences, <, 3),
						"rejected_ratio" => (:get_ratio_rejected, >, 0.25, 0),
						"ratio_edge_proximity" => (:get_ratio_edge_proximity, >, 0.95),
						# "norm_outliers" => (:count_outlier_norms, >, 0, 3),
						"centered_norm" => (:get_maximum_centered_norm, >, 1000)
					      ),
			     "registry" => Dict(
					"global_offsets" => GLOBAL_OFFSETS_PREALIGNMENT
					)
			     )
global PARAMS_ALIGNMENT = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_ALIGNMENT), 
			     "match" => Dict(
					"prematch" => PREMATCH_ALIGNMENT,
					"blockmatch_scale" => BLOCKMATCH_SCALE_ALIGNMENT,
					"block_r" => BLOCK_R_ALIGNMENT, 
					"search_r" => SEARCH_R_ALIGNMENT,
					"bandpass_sigmas" => (2.5, 15),
#					"highpass_sigma" => 20,
#					"lowpass_sigma" => 2.5,
					"depth" => 2,
					"reflexive" => false),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_ALIGNMENT,
					"match_spring_coeff" => MATCH_SPRING_COEFF_ALIGNMENT,
					"ftol_cg" => FTOL_CG_ALIGNMENT,
					"max_iters" => MAX_ITERS_ALIGNMENT,
			     	"use_conjugate_gradient" => USE_CONJUGATE_GRADIENT_ALIGNMENT,
			     	"eta_gd" => ETA_GD_ALIGNMENT,
			     	"ftol_gd" => FTOL_GD_ALIGNMENT,
			     	"eta_newton" => ETA_NEWTON_ALIGNMENT,
			     	"ftol_newton" => FTOL_NEWTON_ALIGNMENT),
#=			     "filter" => Dict(
			     		"sigma_filter_high" => (:get_properties, >, 6, 0.95),
			     		"sigma_filter_mid" => (:get_properties, >, 60, 0.75),
			     		"sigma_filter_low" => (:get_properties, >, 360, 0.50),
			     		"dyn_range_filter" => (:get_properties, <, 0.75, "src_normalized_dyn_range"),
			     		"r_filter" => (:get_properties, <, 0.0275, "r_max"),
			     		# "norm_filter" => (:get_norms_std_sigmas, >, 5),
			     		"kurtosis_filter" => (:get_properties, >, 25, "src_kurtosis"),
			     		"kurtosis_filter_edge" => (:get_properties, <, -1.60, "src_kurtosis"),
					"centered_norm_filter" => (:get_centered_norms, >, 600)
					      ), =#
			     "filter" => Dict(
			     		"sigma_filter_high" => (1,:get_properties, >, 5, 0.95),
			     		"sigma_filter_mid" => (2,:get_properties, >, 75, 0.75),
			     		"sigma_filter_low" => (3,:get_properties, >, 150, 0.50),
			     		"dyn_range_filter" => (4,:get_properties, <, 0.45, "src_normalized_dyn_range"),
			     		"r_filter" => (5,:get_properties, <, 0.0275, "r_max"),
			     		# "norm_filter" => (:get_norms_std_sigmas, >, 5),
			     		"kurtosis_filter" => (6,:get_properties, >, 25, "src_kurtosis"),
			     		"kurtosis_filter_edge" => (7,:get_properties, <, -1.60, "src_kurtosis"),
					#"centered_norm_filter" => (:get_centered_norms, >, 200)
					"centered_norm_filter" => (8,:get_centered_norms, >, 500)
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
			     		"too_few_corresps" => (:count_correspondences, <, 100),
						"rejected_ratio" => (:get_ratio_rejected, >, 0.10),
						"ratio_edge_proximity" => (:get_ratio_edge_proximity, >, 0.95),
						# "norm_outliers" => (:count_outlier_norms, >, 0, 4),
						"centered_norm" => (:get_maximum_centered_norm, >, 250)
					      ),
			     "registry" => Dict(
					"global_offsets" => GLOBAL_OFFSETS_ALIGNMENT
					)
)

global PARAMS_ALIGNMENT_SKIPPED = deepcopy(PARAMS_ALIGNMENT);
PARAMS_ALIGNMENT_SKIPPED["match"]["search_r"] = 1400;
PARAMS_ALIGNMENT_SKIPPED["filter"] = Dict(
			     		"sigma_filter_high" => (:get_properties, >, 7, 0.95),
			     		"sigma_filter_mid" => (:get_properties, >, 100, 0.75),
			     		"sigma_filter_low" => (:get_properties, >, 250, 0.50),
			     		"r_filter" => (:get_properties, <, 0.02, "r_max"),
			     		# "norm_filter" => (:get_norms_std_sigmas, >, 5),
			     		"kurtosis_filter" => (:get_properties, >, 25, "src_kurtosis"),
			     		"kurtosis_filter_edge" => (:get_properties, <, -1.60, "src_kurtosis"),
					"centered_norm_filter" => (:get_centered_norms, >, 700)
					      );

PARAMS_ALIGNMENT_SKIPPED["review"] = Dict(
			     		"too_few_corresps" => (:count_correspondences, <, 100),
						"rejected_ratio" => (:get_ratio_rejected, >, 0.45),
						"ratio_edge_proximity" => (:get_ratio_edge_proximity, >, 0.95),
						"centered_norm" => (:get_maximum_centered_norm, >, BLOCK_R_ALIGNMENT * 1.25)

)

