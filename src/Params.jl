MESH_LENGTH_MONTAGE = 175
GLOBAL_OFFSETS_MONTAGE = true
BLOCKMATCH_SCALE_MONTAGE = 1.0
BLOCK_R_MONTAGE = 60
SEARCH_R_MONTAGE = 88
MONOBLOCK_SCALE_MONTAGE = 1
MONOBLOCK_MATCH_MONTAGE = false
MONOBLOCK_RATIO_MONTAGE = 0.4
MONOBLOCK_PADDING_MONTAGE = 0.0
MESH_SPRING_COEFF_MONTAGE = 1.0
MATCH_SPRING_COEFF_MONTAGE = 3.0 
FTOL_CG_MONTAGE = 1/1000

MESH_LENGTH_PREALIGNMENT = 3000
GLOBAL_OFFSETS_PREALIGNMENT = false
BLOCKMATCH_SCALE_PREALIGNMENT = 0.25
BLOCK_R_PREALIGNMENT = 2000
SEARCH_R_PREALIGNMENT = 800
MONOBLOCK_SCALE_PREALIGNMENT = 0.2
MONOBLOCK_MATCH_PREALIGNMENT = false
MONOBLOCK_RATIO_PREALIGNMENT = 0.4
MONOBLOCK_PADDING_PREALIGNMENT = 0.2
MESH_SPRING_COEFF_PREALIGNMENT = 1.0
MATCH_SPRING_COEFF_PREALIGNMENT = 3.0 
FTOL_CG_PREALIGNMENT = 1/1000

MESH_LENGTH_ALIGNMENT = 500
GLOBAL_OFFSETS_ALIGNMENT = true
BLOCKMATCH_SCALE_ALIGNMENT = 1.0
BLOCK_R_ALIGNMENT = 400
SEARCH_R_ALIGNMENT = 207
MONOBLOCK_SCALE_ALIGNMENT = 1
MONOBLOCK_MATCH_ALIGNMENT = false
MONOBLOCK_RATIO_ALIGNMENT = 0.4
MONOBLOCK_PADDING_ALIGNMENT = 0.0
MESH_SPRING_COEFF_ALIGNMENT = 1.0
MATCH_SPRING_COEFF_ALIGNMENT = 3.0 
FTOL_CG_ALIGNMENT = 1/1000

global PARAMS_MONTAGE = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_MONTAGE), 
			     "match" => Dict(
					"blockmatch_scale" => BLOCKMATCH_SCALE_MONTAGE,
					"block_r" => BLOCK_R_MONTAGE, 
					"search_r" => SEARCH_R_MONTAGE,
					"monoblock_scale" => MONOBLOCK_SCALE_MONTAGE, 
			#		"monoblock_padding" => MONOBLOCK_PADDING_MONTAGE, 
					"monoblock_ratio" => MONOBLOCK_RATIO_MONTAGE, 
					"monoblock_match" => MONOBLOCK_MATCH_MONTAGE),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_MONTAGE,
					"match_spring_coeff" => MATCH_SPRING_COEFF_MONTAGE,
					"ftol_cg" => FTOL_CG_MONTAGE),
			     "filter" => Dict(
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
					      ),
			     "registry" => Dict(
					"global_offsets" => GLOBAL_OFFSETS_MONTAGE
					)
			     )
global PARAMS_PREALIGNMENT = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_PREALIGNMENT), 
			     "match" => Dict(
					"blockmatch_scale" => BLOCKMATCH_SCALE_PREALIGNMENT,
					"block_r" => BLOCK_R_PREALIGNMENT, 
					"search_r" => SEARCH_R_PREALIGNMENT,
					"monoblock_scale" => MONOBLOCK_SCALE_PREALIGNMENT, 
					"monoblock_ratio" => MONOBLOCK_RATIO_PREALIGNMENT, 
					"monoblock_match" => MONOBLOCK_MATCH_PREALIGNMENT),
			     "solve" => Dict(
					"method" => "regularized",
					"lambda" => 0.9,
					"mesh_spring_coeff" => MESH_SPRING_COEFF_PREALIGNMENT,
					"match_spring_coeff" => MATCH_SPRING_COEFF_PREALIGNMENT,
					"ftol_cg" => FTOL_CG_PREALIGNMENT),
			     "filter" => Dict(
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
					      ),
			     "registry" => Dict(
					"global_offsets" => GLOBAL_OFFSETS_PREALIGNMENT
					)
			     )
global PARAMS_ALIGNMENT = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_ALIGNMENT), 
			     "match" => Dict(
					"blockmatch_scale" => BLOCKMATCH_SCALE_ALIGNMENT,
					"block_r" => BLOCK_R_ALIGNMENT, 
					"search_r" => SEARCH_R_ALIGNMENT,
					"monoblock_scale" => MONOBLOCK_SCALE_ALIGNMENT, 
			#		"monoblock_padding" => MONOBLOCK_PADDING_ALIGNMENT, 
					"monoblock_ratio" => MONOBLOCK_RATIO_ALIGNMENT, 
					"monoblock_match" => MONOBLOCK_MATCH_ALIGNMENT),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_ALIGNMENT,
					"match_spring_coeff" => MATCH_SPRING_COEFF_ALIGNMENT,
					"ftol_cg" => FTOL_CG_ALIGNMENT),
			     "filter" => Dict(
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
					      ),
			     "registry" => Dict(
					"global_offsets" => GLOBAL_OFFSETS_ALIGNMENT
					)
			     )
