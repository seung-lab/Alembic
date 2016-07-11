module TestAlignmentTaskInfo

import Julimaps.Cloud.Tasks.AlignmentTask
import JSON

using Base.Test
using CloudTest.MockTasks

function test_no_base_directory()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "baseDirectory")
    @test_throws KeyError AlignmentTask.Info(dict)
end

function test_empty_base_directory()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    dict["baseDirectory"] = "   "
    @test_throws ArgumentError AlignmentTask.Info(dict)
end

function test_no_files()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "files")
    @test_throws KeyError AlignmentTask.Info(dict)
end

function test_empty_files()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    dict["files"] = []
    @test_throws ArgumentError AlignmentTask.Info(dict)
end

function test_no_indices()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "indices")
    @test_throws KeyError AlignmentTask.Info(dict)
end

function test_empty_indices()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    dict["indices"] = []
    @test_throws ArgumentError AlignmentTask.Info(dict)
end

function test_malformed_indices()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    dict["indices"] = [(1, 3, 4, 2), (1)]
    @test_throws ArgumentError AlignmentTask.Info(dict)
end

function test_good_alignment_task_info()
    task = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(task))
    new_task = AlignmentTask.Info(dict)
    @test nothing != new_task
    @test new_task.baseDirectory == TEST_BASE_DIRECTORY
    @test length(new_task.files) == length(TEST_FILES)
    for index in 1:length(TEST_FILES)
        @test new_task.files[index] == TEST_FILES[index]
    end
    @test length(new_task.indices) == length(TEST_INDICES)
    for index in 1:length(TEST_INDICES)
        @test new_task.indices[index] == TEST_INDICES[index]
    end
end

function __init__()
    test_no_base_directory()
    test_empty_base_directory()
    test_no_files()
    test_empty_files()
    test_no_indices()
    test_empty_indices()
    test_malformed_indices()

    test_good_alignment_task_info()
end

end # module TestAlignmentTaskInfo
