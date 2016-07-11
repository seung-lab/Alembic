module MockTasks

import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BasicTask
import Julimaps.Cloud.Tasks.AlignmentTask
import Julimaps.Cloud.Tasks.BlockMatchTask

export TEST_ID, TEST_TASK_NAME, TEST_BASE_DIRECTORY, TEST_FILES
export make_valid_basic_info

export TEST_INDICES
export make_valid_alignment_task_info

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

end # MockTasks
