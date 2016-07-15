module TestBlockMatchTask

using Base.Test
using CloudTest.JulitasksTests.Utils.TestTasks
using CloudTest.JulitasksTests.Utils.MockServices
using CloudTest.RunTests.Utils.TestTasks

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Run.Tasks.BlockMatchTask

function test_execute_undefined_task()
    basic_info = make_valid_basic_info()
    alignment_task_info = make_valid_alignment_task_info()
    block_match_task = BlockMatchTask.BlockMatchTaskDetails(
        basic_info, alignment_task_info)

    datasource = MockDatasourceService()

    found_exception = nothing
    try
        DaemonTask.execute(block_match_task, datasource)
    catch exception
        found_exception = exception
    end
    @test found_exception == nothing
end

function __init__()
    test_execute_undefined_task()
end

end # module TestBlockMatchTask
