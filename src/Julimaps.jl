#module Julimaps
#push!(LOAD_PATH, "./mesh")

### TypeAliases ###############################

export Index
export Triangle, Triangles
export Weight, Weights
export Pairing, Pairings
export Point, Points
export Edges
export BinaryProperty, FloatProperty

typealias Index Tuple{Int64, Int64, Int64, Int64};    # (wafer, section, row, column)

typealias Triangle Tuple{Int64, Int64, Int64};      # index of three points of the triangle for some point
typealias Triangles Array{Triangle, 1};       # index of three points of the triangle for some point

typealias Weight Tuple{Float64, Float64, Float64};    # weights for respective triangle
typealias Weights Array{Weight, 1};       # weights for respective triangle

typealias Pairing Tuple{Int64, Int64};        # useful for abstraction
typealias Pairings Array{Pairing, 1};       # useful for abstraction

typealias Point Array{Float64, 1};        # [i; j]
typealias Points Array{Point, 1};       # array of points
typealias BinaryProperty Array{Bool, 1};      # array of bools

typealias Edges SparseMatrixCSC{Float64, Int64}     # sparse array for edges - columns represent edges and the rows represent the nodes
typealias FloatProperty Array{Float64, 1}


### GLOBAL VARIABLES ###########################
 
export noMatch, noTriangle;

global noMatch = [0; 0; -1];
global noTriangle = (0, 0, 0);

global OVERVIEW_INDEX = -1;
global MONTAGED_INDEX = -2;
global PRE_ALIGNED_INDEX = -3;
global ALIGNED_INDEX = -4;


using Images
using HDF5
using JLD
include("Params.jl")
include("IO.jl")
include("convolve.jl")
include("Mesh.jl")
include("Matches.jl")
include("MeshSet.jl")
include("MeshSolve.jl")
include("MeshSession.jl")

#end
