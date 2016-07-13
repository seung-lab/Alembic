module BucketCacheDatasource

import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache
import Julimaps.Cloud.Services.Datasource


export BucketCacheDatasourceService

type BucketCacheDatasourceService <: Datasource.Service
    remote::Bucket.Service
    cache::Cache.Service
end

function Datasource.pull!(datasource::BucketCacheDatasourceService, key::AbstractString; force::Bool=false)
    if force || !Cache.exists(datasource.cache, key)
        buffer = PipeBuffer()
        Bucket.download(datasource.remote, key, buffer)
        Cache.put!(datasource.cache, key, buffer)
        close(buffer)
    end
end

function Datasource.push!(datasource::BucketCacheDatasourceService, key::AbstractString)
    if !Cache.exists(datasource.cache, key)
        return false
    end

    Bucket.upload(datasource.remote,
        Cache.get(datasource.cache, key), key)
    return true
end

end # module BucketCacheDatasource
