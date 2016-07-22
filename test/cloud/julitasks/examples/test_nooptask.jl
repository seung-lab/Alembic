module TestNoOpTask

using Base.Test
using CloudTest.JulitasksTests.Utils.TestTasks
using CloudTest.JulitasksTests.Utils.MockServices
using Julimaps.Cloud.Julitasks.Examples.NoOpTask

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask

function test_create()
    task = NoOpTaskDetails(make_valid_basic_info(), "TEST")
    task.basicInfo.name = "NOOP_TASK"
    @test task != nothing
end

function test_run()
    task = NoOpTaskDetails(make_valid_basic_info(), "TEST")
    task.basicInfo.name = "NOOP_TASK"

    datasource = MockDatasourceService()
    DaemonTask.run(task, datasource)
end

function __init__()

    test_create()
    test_run()
end

end # module TestNoOpTask
