module TestDaemonTask

using Base.Test
using Julimaps.Cloud.Julitasks.Types
using CloudTest.JulitasksTests.Utils.MockServices

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.DaemonTask

type NewTask <: DaemonTaskDetails end
function test_execute_undefined_task()
    task = NewTask()
    datasource = MockDatasourceService()
    @test_throws ErrorException DaemonTask.execute(task, datasource)
end

function __init__()
    test_execute_undefined_task()
end

end # module TestDaemonTask
