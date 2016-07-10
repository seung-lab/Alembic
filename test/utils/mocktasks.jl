module MockTasks

import Julimaps.Cloud.Tasks.AlignmentTask
import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BlockMatchTask

export TEST_BASE_DIRECTORY, TEST_FILES, TEST_INDICES
export make_valid_alignment_task

export TEST_ID, TEST_TASK_NAME
export make_valid_daemon_task
const TEST_BASE_DIRECTORY = "base_directory"
const TEST_FILES = ["file_1", "file_2"]
const TEST_INDICES = [(1, 2, 3, 4), ( 5, 6, 7, 8)]

function make_valid_alignment_task()
    return AlignmentTask.Details(TEST_BASE_DIRECTORY,
        TEST_FILES, TEST_INDICES)
end

const TEST_ID = 1
const TEST_TASK_NAME = "test_name"

function make_valid_daemon_task()
    return DaemonTask.Details(TEST_ID, TEST_TASK_NAME)
end

end # MockTasks
