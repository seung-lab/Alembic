module TestAlignmentTaskInfo

using Base.Test
using CloudTest.Utils.TestTasks

import Julimaps.Cloud.Tasks.AlignmentTask
import JSON

function test_alignment_task_info_no_indices()
    info = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(info))
    delete!(dict, "indices")
    @test_throws KeyError AlignmentTask.Info(dict)
end

function test_alignment_task_info_empty_indices()
    info = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(info))
    dict["indices"] = []
    @test_throws ArgumentError AlignmentTask.Info(dict)
end

function test_alignment_task_info_malformed_indices()
    info = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(info))
    dict["indices"] = [(1, 3, 4, 2), (1)]
    @test_throws ArgumentError AlignmentTask.Info(dict)
end

function test_alignment_task_info_good()
    info = make_valid_alignment_task_info()
    dict = JSON.parse(JSON.json(info))
    new_info = AlignmentTask.Info(dict)
    @test nothing != new_info
    @test length(new_info.indices) == length(TEST_INDICES)
    for index in 1:length(TEST_INDICES)
        @test new_info.indices[index] == TEST_INDICES[index]
    end
end

function __init__()
    test_alignment_task_info_no_indices()
    test_alignment_task_info_empty_indices()
    test_alignment_task_info_malformed_indices()

    test_alignment_task_info_good()
end

end # module TestAlignmentTaskInfo
