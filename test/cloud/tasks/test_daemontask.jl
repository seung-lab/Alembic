module TestDaemonTask

using Base.Test
using MockTasks

import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BlockMatchTask
import JSON

function test_daemon_task_details_no_id()
    task = make_valid_daemon_task()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "id")
    @test_throws KeyError DaemonTask.Details(dict)
end

function test_daemon_task_details_empty_id()
    task = make_valid_daemon_task()
    dict = JSON.parse(JSON.json(task))
    dict["id"] = ""
    # fail ArgumentError: premature end of integer
    @test_throws ArgumentError DaemonTask.Details(dict)
end

function test_daemon_task_details_bad_id()
    task = make_valid_daemon_task()
    dict = JSON.parse(JSON.json(task))
    dict["id"] = "a"
    #@test_throws ArgumentError: trying to convert "a"
    @test_throws ArgumentError DaemonTask.Details(dict)
end

function test_daemon_task_details_no_name()
    task = make_valid_daemon_task()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "name")
    @test_throws KeyError DaemonTask.Details(dict)
end

function test_daemon_task_details_empty_name()
    task = make_valid_daemon_task()
    dict = JSON.parse(JSON.json(task))
    dict["name"] = ""
    @test_throws ArgumentError DaemonTask.Details(dict)
end

function test_daemon_task_details_good()
    task = make_valid_daemon_task()
    dict = JSON.parse(JSON.json(task))
    details = DaemonTask.Details(dict)
    @test details != nothing
    @test details.id == TEST_ID
    @test details.name == TEST_TASK_NAME
end

type NewTask <: DaemonTask.DaemonTaskDetails end

function test_execute_undefined_task()
    task = NewTask()
    @test_throws ErrorException DaemonTask.execute(task)
end

function __init__()
    test_daemon_task_details_no_id()
    test_daemon_task_details_empty_id()
    test_daemon_task_details_bad_id()
    test_daemon_task_details_no_name()
    test_daemon_task_details_empty_name()
    test_daemon_task_details_good()

    test_execute_undefined_task()
end

end # module TestDaemonTask
