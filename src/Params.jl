SCALE_MONTAGE = 1.0
MESH_LENGTH_MONTAGE = 175
BLOCK_R_MONTAGE = 60
SEARCH_R_MONTAGE = 88
MESH_SPRING_COEFF_MONTAGE = 1.0
MATCH_SPRING_COEFF_MONTAGE = 3.0 
FTOL_CG_MONTAGE = 1/1000

SCALE_PREALIGNMENT = 0.5
MESH_LENGTH_PREALIGNMENT = 2500
BLOCK_R_PREALIGNMENT = 322
SEARCH_R_PREALIGNMENT = 1500
MESH_SPRING_COEFF_PREALIGNMENT = 1.0
MATCH_SPRING_COEFF_PREALIGNMENT = 3.0 
FTOL_CG_PREALIGNMENT = 1/1000
GLOBAL_OFFSETS_PREALIGNMENT = false

SCALE_ALIGNMENT = 1.0
MESH_LENGTH_ALIGNMENT = 500
BLOCK_R_ALIGNMENT = 400
SEARCH_R_ALIGNMENT = 207
MESH_SPRING_COEFF_ALIGNMENT = 1.0
MATCH_SPRING_COEFF_ALIGNMENT = 3.0 
FTOL_CG_ALIGNMENT = 1/1000
GLOBAL_OFFSETS_ALIGNMENT = true

global PARAMS_MONTAGE = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_MONTAGE), 
			     "match" => Dict(
					"scale" => SCALE_MONTAGE,
					"block_r" => BLOCK_R_MONTAGE, 
					"search_r" => SEARCH_R_MONTAGE),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_MONTAGE,
					"match_spring_coeff" => MATCH_SPRING_COEFF_MONTAGE,
					"ftol_cg" => FTOL_CG_MONTAGE),
			     "filter" => Dict(
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
					      )
			     )
global PARAMS_PREALIGNMENT = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_PREALIGNMENT), 
			     "match" => Dict(
					"scale" => SCALE_PREALIGNMENT,
					"block_r" => BLOCK_R_PREALIGNMENT, 
					"search_r" => SEARCH_R_PREALIGNMENT),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_PREALIGNMENT,
					"match_spring_coeff" => MATCH_SPRING_COEFF_PREALIGNMENT,
					"ftol_cg" => FTOL_CG_PREALIGNMENT),
			     "filter" => Dict(
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
					      )
			     )
global PARAMS_ALIGNMENT = Dict(
			     "mesh" => Dict(
					"mesh_length" => MESH_LENGTH_ALIGNMENT), 
			     "match" => Dict(
					"scale" => SCALE_ALIGNMENT,
					"block_r" => BLOCK_R_ALIGNMENT, 
					"search_r" => SEARCH_R_ALIGNMENT),
			     "solve" => Dict(
					"mesh_spring_coeff" => MESH_SPRING_COEFF_ALIGNMENT,
					"match_spring_coeff" => MATCH_SPRING_COEFF_ALIGNMENT,
					"ftol_cg" => FTOL_CG_ALIGNMENT),
			     "filter" => Dict(
					      ),
			     "render" => Dict(
					      ),
			     "review" => Dict(
					      )
			     )
