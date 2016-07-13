module DataSource

import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache

export Service, warm!, get!, put!
type Service
    remote::Bucket.Service
    cache::Cache.Service
end

function pull!(datasource::Service, key::AbstractString, force::Bool=false)
    if force || !Cache.exists(datasource.cache, key)
        buffer = PipeBuffer()
        Bucket.download(datasource.remote, key, buffer)
        Cache.put!(datasource.cache, key, buffer)
        close(buffer)
    end
end

function push!(datasource::Service, key::AbstractString)
    if !Cache.exists(datasource.cache, key)
        return false
    end

    Bucket.upload(datasource.remote,
        Cache.get(datasource.cache, key), key)
    return true
end

end # module DataSource
