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

# Using parametrics because as of 0.4.6 can not promote Array{ASCIIString, 1}
# to Array{AbstractString, 1}
function Datasource.get{String <: AbstractString}(
        datasource::BucketCacheDatasourceService,
        keys::Array{String, 1}; force::Bool=false)
    return map((key) -> Datasource.get(datasource, key; force=force), keys)
end

function Datasource.get(datasource::BucketCacheDatasourceService,
        key::AbstractString; force::Bool=false)
    if force || !Cache.exists(datasource.cache, key)
        stream = Bucket.download(datasource.remote, key)
        Cache.put!(datasource.cache, key, stream)
        close(stream)
    end
    return Cache.get(datasource.cache, key)
end

# Using parametrics because as of 0.4.6 can not promote Array{ASCIIString, 1}
# to Array{AbstractString, 1}
function Datasource.put!{String <: AbstractString}(
        datasource::BucketCacheDatasourceService,
        keys::Array{String, 1})
    return map((key) -> Datasource.put!(datasource, key), keys)
end

function Datasource.put!(datasource::BucketCacheDatasourceService,
        key::AbstractString)
    if !Cache.exists(datasource.cache, key)
        return false
    end

    Bucket.upload(datasource.remote,
        Cache.get(datasource.cache, key), key)
    return true
end

end # module BucketCacheDatasource
