module TestDaemonTask

using Base.Test

import Julimaps.Cloud.JulitasksTests.Tasks.DaemonTask

type NewTask <: DaemonTask.Details end
function test_execute_undefined_task()
    task = NewTask()
    @test_throws ErrorException DaemonTask.execute(task)
end

function __init__()
    test_execute_undefined_task()
end

end # module TestDaemonTask
