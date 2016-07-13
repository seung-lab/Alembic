module BucketCacheDatasource

import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache
import Julimaps.Cloud.Services.Datasource


export BucketCacheDatasourceService

type BucketCacheDatasourceService <: Datasource.Service
    remote::Bucket.Service
    cache::Cache.Service
end

function Datasource.pull!(datasource::BucketCacheDatasourceService,
        keys::Union{AbstractString, Array{AbstractString, 1}}; force::Bool=false)
    if !(typeof(keys) <: Array)
        keys = [keys]
    end
    map((key) -> pull!(datasource, key; force=force), keys)
end

function pull!(datasource::BucketCacheDatasourceService,
    key::AbstractString; force::Bool=false)
    if force || !Cache.exists(datasource.cache, key)
        buffer = PipeBuffer()
        Bucket.download(datasource.remote, key, buffer)
        Cache.put!(datasource.cache, key, buffer)
        close(buffer)
    end
end

function Datasource.push!(datasource::BucketCacheDatasourceService,
        keys::Union{AbstractString, Array{AbstractString, 1}})
    if !(typeof(keys) <: Array)
        keys = [keys]
    end

    return map((key) -> push!(datasource, key), keys)
end

function push!(datasource::BucketCacheDatasourceService,
        key::AbstractString)
    if !Cache.exists(datasource.cache, key)
        return false
    end

    Bucket.upload(datasource.remote,
        Cache.get(datasource.cache, key), key)
    return success
end

end # module BucketCacheDatasource
