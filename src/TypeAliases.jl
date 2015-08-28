module TypeAliases

export Index
export Triangle, Triangles
export Weight, Weights
export Pair, Pairs
export Point, Points
export Edges
export BinaryProperty, FloatProperty


typealias Index Tuple{Int64, Int64, Int64};			# (wafer, section, tile)

typealias Triangle Tuple{Int64, Int64, Int64};			# index of three points of the triangle for some point
typealias Triangles Array{Triangle, 1};				# index of three points of the triangle for some point

typealias Weight Tuple{Float64, Float64, Float64};		# weights for respective triangle
typealias Weights Array{Weight, 1};				# weights for respective triangle

typealias Pair Tuple{Int64, Int64};				# useful for abstraction
typealias Pairs Array{Pair, 1};				# useful for abstraction

typealias Point Array{Float64, 1};				# [i; j]
typealias Points Array{Point, 1};				# array of points
typealias BinaryProperty Array{Bool, 1};			# array of bools

typealias Edges SparseMatrixCSC{Float64, Int64}			# sparse array for edges - columns represent edges and the rows represent the nodes
typealias FloatProperty Array{Float64, 1}

end
