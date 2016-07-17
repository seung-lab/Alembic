module Datasource

using ...Julitasks.Types

export pull!, push!

"""
    pull!(datasource::DatasourceService, key::AbstractString; force::Bool=false)

Pull in input from datasource.
"""
function pull!(datasource::DatasourceService, key::AbstractString; force::Bool=false)
    error("pull! is not implemented for $datasource")
end

"""
    pull!(datasource::DatasourceService, key::Array{AbstractString, 1};

Pull in multiple inputs from datasource.
"""
function pull!(datasource::DatasourceService, key::Array{AbstractString, 1};
        force::Bool=false)
    error("pull! is not implemented for $datasource")
end

"""
    push!(datasource::DatasourceService, key::AbstractString)

Push the output we have in our datasource out
"""
function push!(datasource::DatasourceService, key::AbstractString)
    error("push! is not implemented for $datasource")
end

"""
    push!(datasource::DatasourceService, key::AbstractString)

Push multiple outputs we have in our datasource out
"""
function push!(datasource::DatasourceService, key::Array{AbstractString, 1})
    error("push! is not implemented for $datasource")
end

end # module Datasource
