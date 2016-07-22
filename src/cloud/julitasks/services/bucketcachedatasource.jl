module BucketCacheDatasource

using ...Julitasks.Types

import Julimaps.Cloud.Julitasks.Services.Bucket
import Julimaps.Cloud.Julitasks.Services.Cache
import Julimaps.Cloud.Julitasks.Services.Datasource

export BucketCacheDatasourceService

type BucketCacheDatasourceService <: DatasourceService
    remote::BucketService
    cache::CacheService
end

function Datasource.pull!(datasource::BucketCacheDatasourceService,
        keys::Array{AbstractString, 1}; force::Bool=false)
    return map((key) -> Datasource.pull!(datasource, key; force=force), keys)
end

function Datasource.pull!(datasource::BucketCacheDatasourceService,
        key::AbstractString; force::Bool=false)
    if force || !Cache.exists(datasource.cache, key)
        stream = Bucket.download(datasource.remote, key)
        Cache.put!(datasource.cache, key, stream)
        close(stream)
    end
    return Cache.get(datasource.cache, key)
end

function Datasource.push!(datasource::BucketCacheDatasourceService,
        keys::Array{AbstractString, 1})
    return map((key) -> Datasource.push!(datasource, key), keys)
end

function Datasource.push!(datasource::BucketCacheDatasourceService,
        key::AbstractString)
    if !Cache.exists(datasource.cache, key)
        return false
    end

    Bucket.upload(datasource.remote,
        Cache.get(datasource.cache, key), key)
    return true
end

end # module BucketCacheDatasource
