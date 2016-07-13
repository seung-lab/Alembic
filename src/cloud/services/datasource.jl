module Datasource

export Service, pull!, push!

abstract Service

function pull!(datasource::Service, key::Union{AbstractString,
        Array{AbstractString, 1}}; force::Bool=false)
    error("pull! is not implemented for $datasource")
end

function push!(datasource::Service, key::Union{AbstractString,
        Array{AbstractString, 1}})
    error("push! is not implemented for $datasource")
end

end # module Datasource
