module MockServices

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache
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

function Bucket.download(bucket::MockBucketService, key::AbstractString,
        local_file::IO)
    data = bucket.mockFiles[key]
    write(local_file, data)
end
function Bucket.upload(bucket::MockBucketService, local_file::IO,
        key::AbstractString)
    bucket.mockFiles[key] = readbytes(local_file; all=true)
end

export MockQueueService
type MockQueueService <: Queue.Service
end
function Queue.pop_message(mock_queue::MockQueueService)
end

export MockCacheService
type MockCacheService <: Cache.Service
end

end # module MockServices
