module TestBasicTaskInfo

using Base.Test
using CloudTest.MockTasks

import Julimaps.Cloud.Tasks.BasicTask
import Julimaps.Cloud.Tasks.BlockMatchTask
import JSON

function test_basic_info_info_no_id()
    task = make_valid_basic_info()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "id")
    @test_throws KeyError BasicTask.Info(dict)
end

function test_basic_info_info_empty_id()
    task = make_valid_basic_info()
    dict = JSON.parse(JSON.json(task))
    dict["id"] = ""
    # fail ArgumentError: premature end of integer
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_info_bad_id()
    task = make_valid_basic_info()
    dict = JSON.parse(JSON.json(task))
    dict["id"] = "a"
    #@test_throws ArgumentError: trying to convert "a"
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_info_no_name()
    task = make_valid_basic_info()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "name")
    @test_throws KeyError BasicTask.Info(dict)
end

function test_basic_info_info_empty_name()
    task = make_valid_basic_info()
    dict = JSON.parse(JSON.json(task))
    dict["name"] = ""
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_info_good()
    task = make_valid_basic_info()
    dict = JSON.parse(JSON.json(task))
    info = BasicTask.Info(dict)
    @test info != nothing
    @test info.id == TEST_ID
    @test info.name == TEST_TASK_NAME
end

function __init__()
    test_basic_info_info_no_id()
    test_basic_info_info_empty_id()
    test_basic_info_info_bad_id()
    test_basic_info_info_no_name()
    test_basic_info_info_empty_name()
    test_basic_info_info_good()
end

end # module TestBasicTaskInfo
