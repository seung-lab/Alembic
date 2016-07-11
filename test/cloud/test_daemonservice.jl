module TestDaemonService

using Base.Test
using CloudTest.MockTasks
using CloudTest.MockServices

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Daemon

function test_register_no_execute_method()
    mock_queue = MockQueueService()
    mock_bucket = MockBucketService()
    daemon = Daemon.DaemonService(mock_queue, mock_bucket, 10)
    @test_throws Exception Daemon.register!(daemon,
        TEST_TASK_NAME, MockTaskNoExecute)
end

function test_register_with_execute_method()
    mock_queue = MockQueueService()
    mock_bucket = MockBucketService()
    daemon = Daemon.DaemonService(mock_queue, mock_bucket, 10)
    @test Daemon.register!(daemon, TEST_TASK_NAME, MockTaskWithExecute) !=
        nothing
end

function __init__()
    test_register_no_execute_method()
    test_register_with_execute_method()
end

end #module TestDaemon
