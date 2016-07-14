module Datasource

export Service, pull!, push!

"""
    Datasource.Service

Service to ensure that inputs are loaded in properly. Allows pushing
of output data also
"""
abstract Service

"""
    pull!(datasource::Service, key::AbstractString; force::Bool=false)

Pull in input from datasource.
"""
function pull!(datasource::Service, key::AbstractString; force::Bool=false)
    error("pull! is not implemented for $datasource")
end

"""
    pull!(datasource::Service, key::Array{AbstractString, 1};

Pull in multiple inputs from datasource.
"""
function pull!(datasource::Service, key::Array{AbstractString, 1};
        force::Bool=false)
    error("pull! is not implemented for $datasource")
end

"""
    push!(datasource::Service, key::AbstractString)

Push the output we have in our datasource out
"""
function push!(datasource::Service, key::AbstractString)
    error("push! is not implemented for $datasource")
end

"""
    push!(datasource::Service, key::AbstractString)

Push multiple outputs we have in our datasource out
"""
function push!(datasource::Service, key::Array{AbstractString, 1})
    error("push! is not implemented for $datasource")
end

end # module Datasource
