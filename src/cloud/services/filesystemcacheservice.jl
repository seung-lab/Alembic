module FileSystemCache

import Julimaps.Cloud.Services.Bucket
import Cache

type FileSystemCacheService <: Cache.Service
    baseDirectory::AbstractString
end

function exists(cache::FileSystemCacheService, key::AbstractString)
    return isreadable(baseDirectory * "/" * key)
end

function fetch(cache::FileSystemCacheService, key::AbstractString)
end

function clear(cache::FileSystemCacheService, key::AbstractString)
end

end # module FileSystemCache
