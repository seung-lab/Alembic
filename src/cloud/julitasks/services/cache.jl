module Cache

using ...Julitasks.Types

export Service, exists, put!, get, delete!

"""
    exists(cache::CacheService, key::AbstractString)

check to make sure cache contains the given key.
"""
function exists(cache::CacheService, key::AbstractString)
    error("exists is unimplemented for $cache")
end

"""
    put!(cache::CacheService, key::AbstractString)

Get an IO stream to put an object into
"""
function put!(cache::CacheService, key::AbstractString, value_io::IO)
    error("put is unimplemented for $cache")
end

"""
    get(cache::CacheService, key::AbstractString)

Get an IO stream to the object that is cached with this key.
Returns nothing if key is not found
"""
function get(cache::CacheService, key::AbstractString)
    error("get is unimplemented for $cache")
end

function delete!(cache::CacheService, key::AbstractString)
    error("delete! is unimplemented for $cache")
end

end #module Cache
