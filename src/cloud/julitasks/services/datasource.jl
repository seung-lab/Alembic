module Datasource

using ...Julitasks.Types

export get, put

"""
    get(datasource::DatasourceService, key::AbstractString; force::Bool=false)

Get in input from datasource.
"""
function get(datasource::DatasourceService, key::AbstractString;
        force::Bool=false)
    error("get is not implemented for $datasource")
end

"""
    get{String <: AbstractString}(datasource::DatasourceService,
        key::Array{AbstractString, 1};

Get in multiple inputs from datasource.
"""
# Using parametrics because as of 0.4.6 can not promote Array{ASCIIString, 1}
# to Array{AbstractString, 1}
function get{String <: AbstractString}(datasource::DatasourceService,
        key::Array{String, 1}; force::Bool=false)
    error("get is not implemented for $datasource")
end

"""
    put!(datasource::DatasourceService, key::AbstractString)

Put the output we have in our datasource out
"""
function put!(datasource::DatasourceService, key::AbstractString)
    error("put! is not implemented for $datasource")
end

"""
    put!{String <: AbstractString}(datasource::DatasourceService,
        key::AbstractString)

Put multiple outputs we have in our datasource out
"""
# Using parametrics because as of 0.4.6 can not promote Array{ASCIIString, 1}
# to Array{AbstractString, 1}
function put!{String <: AbstractString}(datasource::DatasourceService,
        key::Array{String, 1})
    error("put! is not implemented for $datasource")
end

end # module Datasource
