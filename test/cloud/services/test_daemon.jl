module TestDaemonService

using Base.Test
using CloudTest.TestTasks
using CloudTest.MockServices

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.DataSource
import Julimaps.Cloud.Services.Daemon
import Julimaps.Cloud.Tasks.DaemonTask

function test_register_no_execute_method()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    @test_throws Exception Daemon.register!(daemon,
        TEST_TASK_NAME, MockTaskNoExecute)
end

function test_register_with_execute_method()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    @test Daemon.register!(daemon, TEST_TASK_NAME, MockTaskWithExecute) !=
        nothing
end

function __init__()
    test_register_no_execute_method()
    test_register_with_execute_method()
end

end #module TestDaemon
