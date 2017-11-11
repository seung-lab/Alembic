module Alembic

    function setblas()
	BLAS.set_num_threads(round(Int64, Sys.CPU_CORES / nprocs()))
    end
    setblas()
    # using OffsetArrays
    using DataFrames
    using FixedPointNumbers
    using IterativeSolvers
    using ImageRegistration
    importall ImageRegistration
    using Optim
    using Distributions
    using Compat
    using StatsBase
    using JSON
    using Images
    using CloudVolume
    using Primes
    if haskey(ENV, "MKLROOT")
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

    FourTupleIndex = NTuple{4, Int64};    # (wafer, section, row, column)
    FourTupleIndices = Tuple{FourTupleIndex, FourTupleIndex};    # (wafer, section, row, column)

    Triangle = NTuple{3, Int64};      # index of three points of the triangle for some point
    Triangles = Array{Triangle, 1};       # index of three points of the triangle for some point

    Weight = Tuple{Float64, Float64, Float64};    # weights for respective triangle
    Weights = Array{Weight, 1};       # weights for respective triangle

    Pairing = Tuple{Int64, Int64};        # useful for abstraction
    Pairings = Array{Pairing, 1};       # useful for abstraction

    global const BinaryProperty = Array{Bool, 1};    	  # array of bools
    global const FloatProperty = Array{Float64, 1}   	# array of floats

    # global constants, independent of deployment
    global const NO_MATCH = [0; 0; -1];
    global const NO_TRIANGLE = (0, 0, 0);
    global const NO_WEIGHTS = (0.0, 0.0, 0.0);
    global const NO_POINT = [typemin(Int64), typemin(Int64)];
    global const NO_RANGE = (0:0, 0:0);

    global const eps = 1e-12;
    global const eps_large = 1e-4;
    global const eps_rec = 1 / eps;

    global const SPLIT_MESH_BASIS = 1000

    # blas_set_num_threads(4); #

    include("core/include.jl")
    include("math/include.jl")
    include("tasks/include.jl")
    include("utilities/author.jl")
    include("utilities/utilities.jl")

    for s in filter(x->string(x)[1]!='#' && x!=:eval, names(Alembic,true))
      # println("Exporting $s")
      eval(Expr(:export, s))
    end

end
