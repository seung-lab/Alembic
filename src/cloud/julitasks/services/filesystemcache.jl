module FileSystemCache

import Julimaps.Cloud.Julitasks.Services.Bucket
import Julimaps.Cloud.Julitasks.Services.Cache

const FOLDER_SEPARATOR = "/"

type FileSystemCacheService <: Cache.Service
    baseDirectory::AbstractString

    FileSystemCacheService(baseDirectory::AbstractString) =
        isempty(strip(baseDirectory)) ?  throw(ArgumentError("Base directory can
            not be empty")) : new(baseDirectory)
end

function sanitize(key::AbstractString)
    return replace(key, r"\.\.", "")
end

function to_filename(cache::FileSystemCacheService, key::AbstractString)
    return "$(cache.baseDirectory)$FOLDER_SEPARATOR$(sanitize(key))"
end

function create_path(full_path_file_name::AbstractString)
    path_end = (length(full_path_file_name) + 1 -
        search(reverse(full_path_file_name), "/")).stop
    if path_end > 0
        mkpath(full_path_file_name[1:path_end])
    end
end

function Cache.exists(cache::FileSystemCacheService, key::AbstractString)
    return isfile(to_filename(cache, key))
end

function Cache.put!(cache::FileSystemCacheService, key::AbstractString,
        value_buffer::IO)
    filename = to_filename(cache, key)
    create_path(filename)
    filestream = open(filename, "w")
    if !iswritable(filestream)
        error("Unable to write to $filename")
    end
    write(filestream, takebuf_array(value_buffer))
    close(filestream)
end

function Cache.get(cache::FileSystemCacheService, key::AbstractString)
    if Cache.exists(cache, key)
        return open(to_filename(cache, key), "r")
    else
        return nothing
    end
end

function Cache.delete!(cache::FileSystemCacheService, key::AbstractString)
    rm(to_filename(cache, key))
end

end # module FileSystemCache
