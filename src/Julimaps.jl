

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
global const NO_WEIGHTS = (0, 0, 0);
global const NO_POINT = [typemin(Int64), typemin(Int64)];
global const NO_RANGE = (0:0, 0:0);
global const NO_INDEX = (0, 0, 0, 0);

global const OVERVIEW_INDEX = -1;
global const PREMONTAGED_INDEX = 1;
global const MONTAGED_INDEX = -2;
global const PREALIGNED_INDEX = -3;
global const ALIGNED_INDEX = -4;
global const FINISHED_INDEX = -5;


if !haskey(ENV, "USER")
  ENV["USER"] = "ubuntu"
end


if ENV["USER"] != "ubuntu"
global const ON_AWS = false;
else
global const ON_AWS = true;
end

global const IMG_ELTYPE = UInt8
global const SUP_SIZE = (50000, 50000)

global const SHARED_SRC_IMAGE = SharedArray(IMG_ELTYPE, SUP_SIZE)
global const SHARED_DST_IMAGE = SharedArray(IMG_ELTYPE, SUP_SIZE)

PKGS_USED = ["HDF5", "JLD", "Images", "ImageView", "Colors", "FixedPointNumbers", "Cairo", "IterativeSolvers", "Optim", "Distributions", "RegERMs", "PyPlot"]

PKGS_USED_CLONABLE = ["https://github.com/JuliaSparse/MKLSparse.jl.git", ""]


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
using RegERMs
#if ENV["USER"] != "dih" && !ON_AWS && !contains(gethostname(), "seunglab")
using PyPlot
using MKLSparse

include("author.jl")
include("Index.jl")
include("registry.jl")
if ON_AWS
  #include("filesystem_formyelin.jl")
    include("dataset_retina.jl")
      include("params_retina.jl")
      #  using AWS
      #  using AWS.S3
      #include("filesystem_aws.jl")
      #include("aws_credentials.jl")
    else
      include("dataset_default.jl")
      include("params_default.jl")
    end
    include("IO.jl")
    include("convolve.jl")
    include("imagecovariance.jl")
    include("Mesh.jl")
    include("Match.jl")
    include("MeshSet.jl")
    include("migrate.jl")
    include("evaluate.jl")
    include("parallelism.jl")
    include("meshconjgrad.jl")
    include("meshsession.jl")
    include("tiletooverview.jl")
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
    include("check.jl")
    include("brushtool.jl")

