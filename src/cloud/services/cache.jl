module Cache

export Service, exists, put!, get, delete!

"""
    CacheService

Fetches cached information to local
"""

abstract Service

"""
    exists(cache::Service, key::AbstractString)

check to make sure cache contains the given key.
"""
function exists(cache::Service, key::AbstractString)
    error("exists is unimplemented for $cache")
end

"""
    put!(cache::Service, key::AbstractString)

Get an IO stream to put an object into
"""
function put!(cache::Service, key::AbstractString, value_buffer::IO)
    error("put is unimplemented for $cache")
end

"""
    get(cache::Service, key::AbstractString)

Get an IO stream to the object that is cached with this key.
Returns nothing if key is not found
"""
function get(cache::Service, key::AbstractString)
    error("clear is unimplemented for $cache")
end

function delete!(cache::Service, key::AbstractString)
    error("clear is unimplemented for $cache")
end

end #module Cache
