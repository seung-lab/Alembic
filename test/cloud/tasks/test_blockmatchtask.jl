module TestBlockMatchTask

using Base.Test
using CloudTest.Utils.TestTasks

import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BlockMatchTask

function test_execute_undefined_task()
    basic_info = make_valid_basic_info()
    alignment_task_info = make_valid_alignment_task_info()
    block_match_task = BlockMatchTask.BlockMatchTaskDetails(
        basic_info, alignment_task_info)
    found_exception = nothing
    try
        DaemonTask.execute(block_match_task)
    catch exception
        found_exception = exception
    end
    @test found_exception == nothing
end

function __init__()
    test_execute_undefined_task()
end

end # module TestBlockMatchTask
