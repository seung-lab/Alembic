module MockServices

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache
import Julimaps.Cloud.Services.Datasource
import Julimaps.Cloud.Tasks.DaemonTask

export MockTaskNoExecute, MockTaskWithExecute
type MockTaskNoExecute <: DaemonTask.Details
end
type MockTaskWithExecute <: DaemonTask.Details
end
function DaemonTask.execute(task::MockTaskWithExecute)
    println("executing")
end

export MockBucketService
type MockBucketService <: Bucket.Service
    mockFiles::Dict{AbstractString, Any}
end
MockBucketService() = MockBucketService(Dict())

# local_file can also be a file location but we're not testing that for now
function Bucket.download(bucket::MockBucketService, key::ASCIIString,
        local_file::Union{ASCIIString, IO})
    data = bucket.mockFiles[key]
    write(local_file, readall(data))
    seekstart(data)
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
type MockQueueService <: Queue.Service
end
function Queue.pop_message(mock_queue::MockQueueService)
end

export MockCacheService
type MockCacheService <: Cache.Service
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
type MockDatasourceService <: Datasource.Service
end

end # module MockServices
