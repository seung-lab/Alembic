module TestTasks

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask

export TEST_ID, TEST_TASK_NAME, TEST_BASE_DIRECTORY, TEST_FILES
export make_valid_basic_info

export TEST_INDICES
export make_valid_alignment_task_info

export MockTaskNoExecute, MockTaskExecute
export make_valid_task_execute, make_valid_task_no_execute

const TEST_ID = 1
const TEST_TASK_NAME = "test_name"
const TEST_BASE_DIRECTORY = "base_directory"
const TEST_FILES = ["file_1", "file_2"]
function make_valid_basic_info()
    return BasicTask.Info(TEST_ID, TEST_TASK_NAME, TEST_BASE_DIRECTORY,
        TEST_FILES)
end

const TEST_INDICES = [(1, 2, 3, 4), ( 5, 6, 7, 8)]
function make_valid_alignment_task_info()
    return AlignmentTask.Info(TEST_INDICES)
end

type MockTaskNoExecute <: DaemonTask.Details
    basicInfo::BasicTask.Info
    payloadInfo::AbstractString
end

type MockTaskExecute <: DaemonTask.Details
    basicInfo::BasicTask.Info
    payloadInfo::AbstractString
end
# register a registered test task
function DaemonTask.execute(task::MockTaskExecute)
    println("executing")
end

function make_valid_task_no_execute()
    return MockTaskNoExecute(make_valid_basic_info(),
        "Mock Task No Execute")
end

function make_valid_task_execute()
    return MockTaskExecute(make_valid_basic_info(),
        "Mock Task Execute")
end

end # module MockTasks
