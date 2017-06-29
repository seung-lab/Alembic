module Alembic

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
using DataFrames
using Colors
using FixedPointNumbers
using Base.Test
using Cairo
using IterativeSolvers
using ImageRegistration
importall ImageRegistration
using Optim
using Distributions
using Compat
using Images
using Graphics
using StatsBase
using JSON
#using SimpleTasks
using Primes
if USE_PYPLOT
  using PyPlot
end
  using Tk
  using ImageView
if !contains(gethostname(), "seung") && !contains(gethostname(), "MacBook")
  using PyCall
  using MKLSparse
end

import Base.Iterators.repeated

abstract type AbstractMesh 		end
abstract type AbstractMatch 		end
#abstract type AbstractProperties 	end

# Aliases for Points - a Point is a 1D array (with two elements), and a Points is a 2xN array where each column is a Point. 
#const Point{T <: Number} = Array{T,1};
#const Points{T <: Number} = Array{T,2};
Point{T <: Number} = Array{T,1};
Points{T <: Number} = Array{T,2};

function columnviews{T}(pts::Points{T})
	@inbounds return map(view, repeated(pts), repeated(:), 1:size(pts, 2))
end

# sparse array for edges - columns represent edges and the rows represent the nodes
#const Edges{T <: Number} = SparseMatrixCSC{T, Int64}
Edges{T <: Number} = SparseMatrixCSC{T, Int64}

#=

# TypeAliases
export Index
export Triangle, Triangles
export Weight, Weights
export Pairing, Pairings
export Point, Points
export Edges
export BinaryProperty, FloatProperty
export Match, Mesh
=#

FourTupleIndex = NTuple{4, Int64};    # (wafer, section, row, column)
FourTupleIndices = Tuple{FourTupleIndex, FourTupleIndex};    # (wafer, section, row, column)

Triangle = NTuple{3, Int64};      # index of three points of the triangle for some point
Triangles = Array{Triangle, 1};       # index of three points of the triangle for some point

Weight = Tuple{Float64, Float64, Float64};    # weights for respective triangle
Weights = Array{Weight, 1};       # weights for respective triangle

Pairing = Tuple{Int64, Int64};        # useful for abstraction
Pairings = Array{Pairing, 1};       # useful for abstraction

#=

const FourTupleIndex = NTuple{4, Int64};    # (wafer, section, row, column)
const FourTupleIndices = Tuple{FourTupleIndex, FourTupleIndex};    # (wafer, section, row, column)

const Triangle = NTuple{3, Int64};      # index of three points of the triangle for some point
const Triangles = Array{Triangle, 1};       # index of three points of the triangle for some point

const Weight = Tuple{Float64, Float64, Float64};    # weights for respective triangle
const Weights = Array{Weight, 1};       # weights for respective triangle

const Pairing = Tuple{Int64, Int64};        # useful for abstraction
const Pairings = Array{Pairing, 1};       # useful for abstraction

  =#


#=
const Point = Array{Float64, 1};        # [i; j]
const Points = Array{Point, 1};       # array of points
=#

global const BinaryProperty = Array{Bool, 1};    	  # array of bools
global const FloatProperty = Array{Float64, 1}   	# array of floats

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

BLAS.set_num_threads(4);


include("core/include.jl")
include("math/include.jl")


include("utilities/author.jl")
# include("utilities/parallelism.jl")
include("utilities/utilities.jl")


#include("datasets/dataset_cremi.jl")
#include("params/params_default.jl")

#include("datasets/dataset_s1.jl")
include("datasets/dataset_myelin.jl")
include("params/params_myelin.jl")
# include("datasets/dataset_pinky.jl")
# include("params/params_pinky.jl")
# include("datasets/dataset_wholefish.jl")
#include("datasets/dataset_zebrafish.jl")
#include("datasets/dataset_default.jl")
#include("params/params_zebrafish.jl")
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

#include("archive/evaluate.jl")
#include("archive/check.jl")
#include("archive/meshsession.jl")
#include("archive/tiletooverview.jl")
#include("archive/data_export.jl")

include("import/premontage.jl")
include("import/import_AIBS_TEM.jl")
include("import/import_AIBS_SEM.jl")
#include("import/old_import.jl")

#include("review/review.jl")
#include("review/visualize.jl")
#include("review/draw.jl")
#include("review/inspect.jl")
#include("review/player.jl")
#include("review/brushtool.jl")
#include("review/cpselect.jl")
#=
include("tasks/tasks_env.jl")
include("tasks/ImportTask.jl")
include("tasks/BlockMatchTask.jl")
include("tasks/RenderTask.jl")
include("tasks/SolveTask.jl")
include("tasks/MaskTask.jl")
include("tasks/ThumbnailTask.jl")
include("tasks/RenderReviewTask.jl")
include("tasks/SaveStackTask.jl")
include("tasks/CubeStackTask.jl")
include("tasks/awsscheduler.jl")

include("utilities/migrate.jl")
=#

for s in filter(x->string(x)[1]!='#' && x!=:eval, names(Alembic,true))

      println("Exporting $s")

          eval(Expr(:export, s))

	end

end
