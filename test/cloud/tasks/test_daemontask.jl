module TestDaemonTask

using Base.Test

import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BlockMatchTask
import JSON

function test_daemon_task_details_no_id()
    message = "{\"name\": \"FAKE_TASK\"}"
    dict = JSON.parse(message)
    @test_throws KeyError DaemonTask.Details(dict)
end

function test_daemon_task_details_empty_id()
    message = "{\"id\": \"\",\"name\": \"FAKE_TASK\"}"
    dict = JSON.parse(message)
    # fail ArgumentError: premature end of integer
    @test_throws ArgumentError DaemonTask.Details(dict)
end

function test_daemon_task_details_bad_id()
    message = "{\"id\": \"a\",\"name\": \"FAKE_TASK\"}"
    dict = JSON.parse(message)
    #@test_throws ArgumentError: trying to convert a
    @test_throws ArgumentError DaemonTask.Details(dict)
end

function test_daemon_task_details_no_name()
    message = "{\"id\": \"2\"}"
    dict = JSON.parse(message)
    @test_throws KeyError DaemonTask.Details(dict)
end

function test_daemon_task_details_empty_name()
    message = "{\"id\": \"2\",\"name\": \"\"}"
    dict = JSON.parse(message)
    @test_throws ArgumentError DaemonTask.Details(dict)
end

function test_daemon_task_details_good()
    message = "{\"id\": \"2\",\"name\": \"FAKE_TASK\"}"
    dict = JSON.parse(message)
    detail = DaemonTask.Details(dict)
    @test detail.id == 2
    @test detail.name == "FAKE_TASK"
    dict["id"] = "2"
    detial = DaemonTask.Details(dict)
    @test detail.id == 2
    @test detail.name == "FAKE_TASK"
end

type NewTask <: DaemonTask.DaemonTaskDetails end

function test_execute_undefined_task()
    task = NewTask()
    @test_throws ErrorException DaemonTask.execute(task)
end

#=
 =function test_parse_empty()
 =    message = ""
 =    @test_throws ErrorException DaemonTask.parse(message)
 =    return true
 =end
 =
 =function test_parse_bad_convert()
 =    message = "test message"
 =    @test_throws ErrorException DaemonTask.parse(message)
 =    return true
 =end
 =
 =function get_test_block_match_task()
 =    return BlockMatchTask.BlockMatchTaskDetails(DaemonTask.Details(
 =        DaemonTask.TASK_TYPE_BLOCK_MATCH,
 =        "/",
 =        ["file1", "file2"],
 =        [(1,1,1,1), (2,2,2,2)]
 =        )
 =    )
 =end
 =
 =function test_to_daemon_task_good()
 =    task = get_test_block_match_task()
 =    json_string = JSON.json(task)
 =    json_dict = JSON.parse(json_string)
 =
 =    parsed_task = DaemonTask.to_daemon_task(json_dict)
 =
 =    @test parsed_task.details.taskType == task.details.taskType
 =    @test parsed_task.details.baseDirectory == task.details.baseDirectory
 =    @test length(parsed_task.details.files) == length(task.details.files)
 =    for i in 1:length(task.details.files)
 =        @test parsed_task.details.files[i] == task.details.files[i]
 =    end
 =    @test length(parsed_task.details.indices) == length(task.details.indices)
 =    for i in 1:length(task.details.indices)
 =        @test parsed_task.details.indices[i] == task.details.indices[i]
 =    end
 =
 =    return true
 =end
 =
 =function test_to_daemon_task_error()
 =    task = get_test_block_match_task()
 =    json_string = JSON.json(task)
 =    json_dict = JSON.parse(json_string)
 =
 =    details_dict = json_dict["details"]
 =
 =    delete!(details_dict, "taskType")
 =    @test_throws KeyError DaemonTask.to_daemon_task(json_dict)
 =
 =    details_dict["taskType"] = "BAD TYPE"
 =    @test_throws ErrorException DaemonTask.to_daemon_task(json_dict)
 =    details_dict["taskType"] = DaemonTask.TASK_TYPE_BLOCK_MATCH
 =
 =    delete!(details_dict, "files")
 =    @test_throws KeyError DaemonTask.to_daemon_task(json_dict)
 =    return true
 =end
 =
 =#
function __init__()
    test_daemon_task_details_no_id()
    test_daemon_task_details_empty_id()
    test_daemon_task_details_bad_id()
    test_daemon_task_details_no_name()
    test_daemon_task_details_empty_name()
    test_daemon_task_details_good()

    test_execute_undefined_task()
    #=@test test_parse_empty()=#
    #=@test test_parse_bad_convert()=#
    #=@test test_to_daemon_task_good()=#
    #=@test test_to_daemon_task_error()=#
end

end # module TestDaemonTask
