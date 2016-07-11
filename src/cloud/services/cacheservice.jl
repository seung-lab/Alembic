module Cache

"""
    CacheService

Fetches cached information to local
"""

abstract Service

function exists(cache::Service, key::AbstractString)
    error("exists is unimplemented for $cache")
end

function fetch(cache::Service, key::AbstractString)
    error("clear is unimplemented for $cache")
end

function clear(cache::Service, key::AbstractString)
    error("clear is unimplemented for $cache")
end

end #module Cache
