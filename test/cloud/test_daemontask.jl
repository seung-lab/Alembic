
module TestDaemonTask

using Base.Test

import JSON
import Julimaps.Cloud.DaemonTask

function test_parse_empty()
    message = ""
    @test_throws ErrorException DaemonTask.parse(message)
    return true
end

function test_parse_bad_convert()
    message = "test message"
    @test_throws ErrorException DaemonTask.parse(message)
    return true
end

function get_test_block_match_task()
    return DaemonTask.BlockMatchTask(DaemonTask.Details(
        DaemonTask.TASK_TYPE_BLOCK_MATCH,
        "/",
        ["file1", "file2"],
        [(1,1,1,1), (2,2,2,2)]
        )
    )
end

function test_to_daemon_task_good()
    task = get_test_block_match_task()
    json_string = JSON.json(task)
    json_dict = JSON.parse(json_string)

    parsed_task = DaemonTask.to_daemon_task(json_dict)

    @test parsed_task.details.taskType == task.details.taskType
    @test parsed_task.details.baseDirectory == task.details.baseDirectory
    @test length(parsed_task.details.files) == length(task.details.files)
    for i in 1:length(task.details.files)
        @test parsed_task.details.files[i] == task.details.files[i]
    end
    @test length(parsed_task.details.indices) == length(task.details.indices)
    for i in 1:length(task.details.indices)
        @test parsed_task.details.indices[i] == task.details.indices[i]
    end

    return true
end

function test_to_daemon_task_error()
    task = get_test_block_match_task()
    json_string = JSON.json(task)
    json_dict = JSON.parse(json_string)

    details_dict = json_dict["details"]

    delete!(details_dict, "taskType")
    @test_throws KeyError DaemonTask.to_daemon_task(json_dict)

    details_dict["taskType"] = "BAD TYPE"
    @test_throws ErrorException DaemonTask.to_daemon_task(json_dict)
    details_dict["taskType"] = DaemonTask.TASK_TYPE_BLOCK_MATCH

    delete!(details_dict, "files")
    @test_throws KeyError DaemonTask.to_daemon_task(json_dict)
    return true
end

function __init__()
    @test test_parse_empty()
    @test test_parse_bad_convert()
    @test test_to_daemon_task_good()
    @test test_to_daemon_task_error()
end

end # module TestDaemonTask
