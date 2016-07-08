module TestAlignmentTask

import Julimaps.Cloud.Tasks.AlignmentTask
import JSON

using Base.Test

const TEST_BASE_DIRECTORY = "base_directory"
const TEST_FILES = ["file_1", "file_2"]
const TEST_INDICES = [(1, 2, 3, 4), ( 5, 6, 7, 8)]

function make_valid_alignment_task()
    return AlignmentTask.Details(TEST_BASE_DIRECTORY,
        TEST_FILES, TEST_INDICES)
end

function test_no_base_directory()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "baseDirectory")
    @test_throws KeyError AlignmentTask.Details(dict)
end

function test_empty_base_directory()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    dict["baseDirectory"] = "   "
    @test_throws ArgumentError AlignmentTask.Details(dict)
end

function test_no_files()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "files")
    @test_throws KeyError AlignmentTask.Details(dict)
end

function test_empty_files()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    dict["files"] = []
    @test_throws ArgumentError AlignmentTask.Details(dict)
end

function test_no_indices()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    delete!(dict, "indices")
    @test_throws KeyError AlignmentTask.Details(dict)
end

function test_empty_indices()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    dict["indices"] = []
    @test_throws ArgumentError AlignmentTask.Details(dict)
end

function test_malformed_indices()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    dict["indices"] = [(1, 3, 4, 2), (1)]
    @test_throws ArgumentError AlignmentTask.Details(dict)
end

function test_good_alignment_task()
    task = make_valid_alignment_task()
    dict = JSON.parse(JSON.json(task))
    @test nothing != AlignmentTask.Details(dict)
end

function __init__()
    test_no_base_directory()
    test_empty_base_directory()
    test_no_files()
    test_empty_files()
    test_no_indices()
    test_empty_indices()
    test_malformed_indices()

    test_good_alignment_task()
end
end # module TestAlignmentTask
