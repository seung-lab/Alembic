module TestBasicTaskInfo

using Base.Test
using CloudTest.JulitasksTests.Utils.TestTasks

import Julimaps.Cloud.Julitasks.Tasks.BasicTask
import JSON

function test_basic_info_no_id()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "id")
    BasicTask.Info(dict)
end

function test_basic_info_empty_id()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    dict["id"] = ""
    # fail ArgumentError: premature end of integer
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_bad_id()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    dict["id"] = "a"
    #@test_throws ArgumentError: trying to convert "a"
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_no_name()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "name")
    @test_throws KeyError BasicTask.Info(dict)
end

function test_basic_info_empty_name()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    dict["name"] = ""
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_no_base_directory()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "baseDirectory")
    @test_throws KeyError BasicTask.Info(dict)
end

function test_basic_info_empty_base_directory()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    dict["baseDirectory"] = "   "
    @test_throws ArgumentError BasicTask.Info(dict)
end

function test_basic_info_no_inputs()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "inputs")
    @test_throws KeyError BasicTask.Info(dict)
end

function test_basic_info_empty_inputs()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    dict["inputs"] = []
    @test_throws ArgumentError BasicTask.Info(dict)
end


function test_basic_info_info_good()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    new_info = BasicTask.Info(dict)
    @test new_info != nothing
    @test new_info.id == TEST_ID
    @test new_info.name == TEST_TASK_NAME
    @test new_info.baseDirectory == TEST_BASE_DIRECTORY
    @test length(new_info.inputs) == length(TEST_INPUTS)
    for index in 1:length(TEST_INPUTS)
        @test new_info.inputs[index] == TEST_INPUTS[index]
    end
end

function __init__()
    test_basic_info_no_id()
    test_basic_info_empty_id()
    test_basic_info_bad_id()
    test_basic_info_no_name()
    test_basic_info_empty_name()
    test_basic_info_no_base_directory()
    test_basic_info_empty_base_directory()
    test_basic_info_no_inputs()
    test_basic_info_empty_inputs()

    test_basic_info_info_good()
end

end # module TestBasicTaskInfo
