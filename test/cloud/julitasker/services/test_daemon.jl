module TestDaemonService

using Base.Test
using CloudTest.JulitasksTest.Utils.TestTasks
using CloudTest.JulitasksTest.Utils.MockServices

import Julimaps.Cloud.Julitasks.Services.Queue
import Julimaps.Cloud.Julitasks.Services.Bucket
import Julimaps.Cloud.Julitasks.Services.DataSource
import Julimaps.Cloud.Julitasks.Services.Daemon
import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import JSON

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
    @test Daemon.register!(daemon, TEST_TASK_NAME, MockTaskExecute) !=
        nothing
end

function test_parse_empty()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    @test_throws ArgumentError Daemon.parse(daemon, " ")
end

function test_parse_no_basic_info()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    task = make_valid_task_execute()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "basicInfo")
    @test_throws KeyError Daemon.parse(daemon, JSON.json(dict))
end

function test_parse_no_payload_info()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    task = make_valid_task_execute()
    Daemon.register!(daemon, TEST_TASK_NAME, typeof(task))
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "payloadInfo")
    @test_throws KeyError Daemon.parse(daemon, JSON.json(dict))
end

function test_parse_task_not_registered()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    task = make_valid_task_execute()
    dict = JSON.parse(JSON.json(task))

    @test_throws ErrorException Daemon.parse(daemon, JSON.json(dict))
end

function test_parse_good()
    daemon = Daemon.Service(MockQueueService(),
        MockBucketService(),
        MockDatasourceService(), 10)
    task = make_valid_task_execute()
    Daemon.register!(daemon, TEST_TASK_NAME, typeof(task))
    dict = JSON.parse(JSON.json(task))
    parsed_task = Daemon.parse(daemon, JSON.json(dict))
    @test parsed_task != nothing
    @test parsed_task.basicInfo.name == task.basicInfo.name
    @test parsed_task.basicInfo.id == task.basicInfo.id
    @test parsed_task.payloadInfo == task.payloadInfo

end

function __init__()
    test_register_no_execute_method()
    test_register_with_execute_method()

    test_parse_empty()
    test_parse_no_basic_info()
    test_parse_task_not_registered()
    test_parse_good()
end

end #module TestDaemon
