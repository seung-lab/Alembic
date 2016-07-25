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
        datasource::BucketCacheDatasourceService, keys::Array{String, 1};
        override_cache::Bool=false)
    return map((key) -> Datasource.get(datasource, key;
    override_cache=override_cache), keys)
end

function Datasource.get(datasource::BucketCacheDatasourceService,
        key::AbstractString; override_cache::Bool=false)
    if override_cache || !Cache.exists(datasource.cache, key)
        stream = Bucket.download(datasource.remote, key)
        Cache.put!(datasource.cache, key, stream)
        close(stream)
    end
    return Cache.get(datasource.cache, key)
end

# Using parametrics because as of 0.4.6 can not promote Array{ASCIIString, 1}
# to Array{AbstractString, 1}
function Datasource.put!{String <: AbstractString, I <: IO}(
        datasource::BucketCacheDatasourceService, keys::Array{String, 1},
        new_values::Array{I, 1}; only_cache::Bool=false)
    return map((index) -> Datasource.put!(datasource, keys[index],
        new_values[index]; only_cache=only_cache), 1:length(keys))
end

function Datasource.put!{String <: AbstractString}(
        datasource::BucketCacheDatasourceService, keys::Array{String, 1};
        only_cache::Bool=false)
    return map((index) -> Datasource.put!(datasource, keys[index],
        nothing; only_cache=only_cache), 1:length(keys))
end

function Datasource.put!(datasource::BucketCacheDatasourceService,
        key::AbstractString, new_value::Union{IO, Void}=nothing;
        only_cache::Bool=false)
    if new_value != nothing
        seekstart(new_value)
        Cache.put!(datasource.cache, key, new_value)
    end

    if !Cache.exists(datasource.cache, key)
        return false
    end

    if !only_cache
        Bucket.upload(datasource.remote,
            Cache.get(datasource.cache, key), key)
    end
    return true
end

end # module BucketCacheDatasource
