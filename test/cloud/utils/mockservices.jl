module MockServices

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Daemon

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
end

export MockQueueService
type MockQueueService <: Queue.Service
end
function Queue.pop_message(mock_queue::MockQueueService)
end

end # module MockServices
