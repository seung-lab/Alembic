module TestDaemonService

using Base.Test

import Julimaps.Cloud.Queue
import Julimaps.Cloud.Bucket
import Julimaps.Cloud.DaemonTask
import Julimaps.Cloud.Daemon

import AWS
export pop_message

type MockTaskType <: DaemonTask.DaemonTaskDetails
end

function DaemonTask.execute(::MockTaskType)
    println("Running Mock Task")
end

type MockQueueService <: Queue.QueueService
    _module::Module
    MockQueueService() = new(TestDaemonService)
end
type MockBucketService <: Bucket.BucketService
    _module::Module
    MockBucketService() = new(TestDaemonService)
end

function Queue.pop_message(mock_queue::MockQueueService)
    return "hi"
end

function test_register()
    mock_queue = MockQueueService()
    mock_bucket = MockBucketService()
    daemon = Daemon.DaemonService(mock_queue, mock_bucket, 10)
    Daemon.register(daemon, MockTaskType)
    println("test_register")
    Daemon.run(daemon)
end
function __init__()
    @test test_register()
end

end #module TestDaemon
