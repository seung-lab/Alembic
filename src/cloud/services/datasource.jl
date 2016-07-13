module Datasource

export Service, pull!, push!

abstract Service

function pull!(datasource::Service, key::AbstractString; force::Bool=false)
    error("pull! is not implemented for $datasource")
end

function push!(datasource::Service, key::AbstractString)
    error("push! is not implemented for $datasource")
end

end # module Datasource
