#module Alembic

if !haskey(ENV, "USER")
  ENV["USER"] = "ubuntu"
end

if ENV["USER"] != "ubuntu"
  global const ON_CLOUD = false;
else
  global const ON_CLOUD = true;
end
#=
if contains(gethostname(), "seunglab") || contains(gethostname(), "seungom") || ENV["USER"] == "dih"
  global const USE_PYPLOT = false;
else
  global const USE_PYPLOT = true;
end
=#

#global const USE_PYPLOT = true;
global const USE_PYPLOT = false;


PKGS_USED = ["HDF5", "JLD", "Images", "ImageView", "Colors", "FixedPointNumbers", "Cairo", "IterativeSolvers", "Optim", "Distributions", "RegERMs", "PyPlot"]

PKGS_USED_CLONABLE = ["https://github.com/JuliaSparse/MKLSparse.jl.git", 
                      "https://github.com/seung-lab/ImageRegistration.git", 
		      "https://github.com/madeleineudell/ParallelSparseMatMul.jl.git",
                      "https://github.com/macrintr/ImageView.jl.git"]

using HDF5
using JLD
using Colors
using FixedPointNumbers
using Base.Test
using Cairo
using IterativeSolvers
using ImageRegistration
using Optim
using Distributions
using Compat
using Images
using Graphics
using StatsBase
using JSON
using SimpleTasks
if VERSION != v"0.4.6"
using Primes
end
if USE_PYPLOT
  using PyPlot
end
if !ON_CLOUD #!(contains(gethostname(), "seunglab") || contains(gethostname(), "seungom"))
  using Tk
  using ImageView
end
if !contains(gethostname(), "seung")
  using PyCall
  using MKLSparse
end

import Base.filter!

# TypeAliases
export Index
export Triangle, Triangles
export Weight, Weights
export Pairing, Pairings
export Point, Points
export Edges
export BinaryProperty, FloatProperty
export Match, Mesh

typealias Index Tuple{Int64, Int64, Int64, Int64};    # (wafer, section, row, column)
typealias Indices Tuple{Index, Index};    # (wafer, section, row, column)

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
global const NO_WEIGHTS = (0.0, 0.0, 0.0);
global const NO_POINT = [typemin(Int64), typemin(Int64)];
global const NO_RANGE = (0:0, 0:0);
global const NO_INDEX = (0, 0, 0, 0);

global const OVERVIEW_INDEX = 0;
global const PREMONTAGED_INDEX = 1;
global const MONTAGED_INDEX = -2;
global const PREALIGNED_INDEX = -3;
global const ALIGNED_INDEX = -4;
global const FINISHED_INDEX = -5;

global const eps = 1e-12;
global const eps_large = 1e-4;
global const eps_rec = 1 / eps;

blas_set_num_threads(4);



include("math/meshconjgrad.jl")
include("math/meshgradnewton.jl")
#include("math/convolve.jl") # legacy convolution code
include("math/convolve_inplace.jl")
include("math/imagecovariance.jl")

include("utilities/author.jl")
# include("utilities/parallelism.jl")
include("utilities/utilities.jl")

include("core/registry.jl")
include("core/IO.jl")
include("core/Mesh.jl")
include("core/Match.jl")
include("core/filter.jl")
include("core/MeshSet.jl")
include("core/solve.jl")

#include("datasets/dataset_myelin.jl")
#include("params/params_myelin.jl")
#include("datasets/dataset_pinky.jl")
#include("params/params_pinky.jl")
include("datasets/dataset_zebrafish.jl")
#  include("datasets/dataset_default.jl")
  include("params/params_default.jl")
#=
if ON_CLOUD
  include("datasets/dataset_zebrafish.jl")
  include("datasets/dataset_pinky.jl")
  #include("params_default.jl")
  include("params/params_pinky.jl")
else
  include("datasets/dataset_default.jl")
  # include("dataset_zebrafish.jl")
  include("params/params_pinky.jl")
end
=#
include("datasets/dataset_common.jl")

include("render/render.jl")
include("render/imageprocessing.jl")

include("archive/evaluate.jl")
include("archive/check.jl")
include("archive/meshsession.jl")
include("archive/tiletooverview.jl")
include("archive/data_export.jl")

include("import/premontage.jl")
include("import/import.jl")
#include("import/old_import.jl")

include("review/review.jl")
include("review/visualize.jl")
include("review/draw.jl")
if !ON_CLOUD #!(contains(gethostname(), "seunglab") || contains(gethostname(), "seungom"))
  include("review/inspect.jl")
  include("review/player.jl")
  include("review/brushtool.jl")
  include("review/cpselect.jl")
end
  include("tasks/tasks_env.jl")
  include("tasks/ImportTask.jl")
  include("tasks/BlockMatchTask.jl")
  include("tasks/RenderTask.jl")
  include("tasks/SolveTask.jl")
  include("tasks/awsscheduler.jl")

include("utilities/migrate.jl")

#end
