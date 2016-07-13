module TestBasicTaskInfo

using Base.Test
using CloudTest.TestTasks

import Julimaps.Cloud.Tasks.BasicTask
import Julimaps.Cloud.Tasks.BlockMatchTask
import JSON

function test_basic_info_no_id()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "id")
    @test_throws KeyError BasicTask.Info(dict)
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

function test_basic_info_no_files()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "files")
    @test_throws KeyError BasicTask.Info(dict)
end

function test_basic_info_empty_files()
    info = make_valid_basic_info()
    dict = JSON.parse(JSON.json(info))
    dict["files"] = []
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
    @test length(new_info.files) == length(TEST_FILES)
    for index in 1:length(TEST_FILES)
        @test new_info.files[index] == TEST_FILES[index]
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
    test_basic_info_no_files()
    test_basic_info_empty_files()

    test_basic_info_info_good()
end

end # module TestBasicTaskInfo
