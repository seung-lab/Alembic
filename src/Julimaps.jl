

# TypeAliases


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
typealias BinaryProperty Array{Bool, 1};    	  # array of bools

typealias Edges SparseMatrixCSC{Float64, Int64}     # sparse array for edges - columns represent edges and the rows represent the nodes
typealias FloatProperty Array{Float64, 1}   	# array of floats

# global constants, independent of deployment

global const NO_MATCH = [0; 0; -1];
global const NO_TRIANGLE = (0, 0, 0);
global const NO_RANGE = (0:0, 0:0);
global const NO_INDEX = (0, 0, 0, 0);

global const OVERVIEW_INDEX = -1;
global const PREMONTAGED_INDEX = 1;
global const MONTAGED_INDEX = -2;
global const PREALIGNED_INDEX = -3;
global const ALIGNED_INDEX = -4;

global const IMG_ELTYPE = UInt8
global const SUP_SIZE = (40000, 40000)

global const SHARED_SRC_IMAGE = SharedArray(IMG_ELTYPE, SUP_SIZE)
global const SHARED_DST_IMAGE = SharedArray(IMG_ELTYPE, SUP_SIZE)


# dependencies
using HDF5
using JLD
using Images
using ImageView
using Colors
using FixedPointNumbers
using Base.Test
using Cairo
using IterativeSolvers
using ImageRegistration
using Optim
using Distributions
using PyPlot

include("Index.jl")
include("registry.jl")
include("filesystem.jl")
include("IO.jl")
include("Params.jl")
include("convolve.jl")
include("Mesh.jl")
include("Match.jl")
include("MeshSet.jl")
include("MeshSolve.jl")
include("MeshConjGrad.jl")
include("MeshSession.jl")
include("TileToOverview.jl")
include("prealign.jl")
include("imageprocessing.jl")
include("render.jl")
include("review.jl")
include("solve.jl")
include("visualize.jl")
include("utilities.jl")
include("transforms.jl")
include("draw.jl")
include("visualize.jl")
include("player.jl")
include("inspect.jl")

#end
