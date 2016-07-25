module MockServices

using Julimaps.Cloud.Julitasks.Types
using CloudTest.JulitasksTests.Utils.TestTasks

import Julimaps.Cloud.Julitasks.Services.Queue
import Julimaps.Cloud.Julitasks.Services.Bucket
import Julimaps.Cloud.Julitasks.Services.Cache
import Julimaps.Cloud.Julitasks.Services.Datasource

export MockBucketService
type MockBucketService <: BucketService
    mockFiles::Dict{AbstractString, Any}
end
MockBucketService() = MockBucketService(Dict())

# local_file can also be a file location but we're not testing that for now
function Bucket.download(bucket::MockBucketService, key::ASCIIString,
        local_file::Union{ASCIIString, IO, Void})
    data = bucket.mockFiles[key]
    if local_file != nothing
        write(local_file, readall(data))
        return local_file
    end
    seekstart(data)

    return IOBuffer(data)
end
# local_file can also be a file location but we're not testing that for now
function Bucket.upload(bucket::MockBucketService,
    local_file::Union{ASCIIString, IO}, key::ASCIIString)
    data = IOBuffer()
    write(data, readbytes(local_file))
    seekstart(data)
    bucket.mockFiles[key] = data
end

export MockQueueService
type MockQueueService <: QueueService
end
function Queue.pop_message(mock_queue::MockQueueService)
end

export MockCacheService
type MockCacheService <: CacheService
    mockValues::Dict{AbstractString, IO}
end
MockCacheService() = MockCacheService(Dict())
function Cache.exists(cache::MockCacheService, key::AbstractString)
    return haskey(cache.mockValues, key)
end
function Cache.put!(cache::MockCacheService, key::AbstractString,
        value_buffer::IO)
    data = IOBuffer()
    write(data, readbytes(value_buffer))
    seekstart(data)
    cache.mockValues[key] = data
end
function Cache.get(cache::MockCacheService, key::AbstractString)
    if haskey(cache.mockValues, key)
        data = cache.mockValues[key]
        seekstart(data)
        return data
    else
        return nothing
    end
end
function Cache.delete!(cache::MockCacheService, key::AbstractString)
    delete!(cache.mockValues, key)
end

export MockDatasourceService
type MockDatasourceService <: DatasourceService
end
function Datasource.get{String <: AbstractString}(
        datasource::MockDatasourceService, keys::Array{String, 1})
    println("getting $keys")
end
function Datasource.get(datasource::MockDatasourceService,
        key::AbstractString)
    println("getting $key")
end
function Datasource.put!{String <: AbstractString}(
        datasource::MockDatasourceService, keys::Array{String, 1})
    println("putting $keys")
end
function Datasource.put!(datasource::MockDatasourceService,
        key::AbstractString)
    println("putting $key")
end

end # module MockServices
