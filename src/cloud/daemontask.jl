module DaemonTask

import JSON

export DaemonTaskDetails
export BlockMatchTask
export RenderTask
export parse_task

# These types identify the task that we like to perform
abstract DaemonTaskDetails

type Details
    id::Int
    taskType::AbstractString
    payload::AbstractString
end

type RenderTask <: DaemonTaskDetails
    details::Details
end

#=
 =Set the constant identifiers for each task.
 =Use these names to create a lookup table for each constructor of the
 =corresponding task
 =#
const TASK_TYPE_BLOCK_MATCH = "BLOCK_MATCH"
const TASK_TYPE_RENDER = "RENDER"
const TASKS = Dict()
TASKS[TASK_TYPE_BLOCK_MATCH] = BlockMatchTask
TASKS[TASK_TYPE_RENDER] = RenderTask


function execute(render_task_details::RenderTask)
    print("Running Rendering with indicies $(render_task_details.indices) for
        $(render_task_details.base_directory)")
end


end # module Task
