module TestTasks

import Julimaps.Cloud.Run.Tasks.AlignmentTask
import Julimaps.Cloud.Run.Tasks.BlockMatchTask

export TEST_INDICES
export make_valid_alignment_task_info

const TEST_INDICES = [(1, 2, 3, 4), ( 5, 6, 7, 8)]
function make_valid_alignment_task_info()
    return AlignmentTask.Info(TEST_INDICES)
end

end # module TestTasks
