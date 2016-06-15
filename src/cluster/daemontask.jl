module DaemonTask

import JSON

export DaemonTask
export BlockMatchTask
export RenderTask
export parse_task

# These types identify the task that we like to perform
abstract DaemonTask

type BlockMatchTask <: DaemonTask
    taskId::AbstractString
    indices::AbstractString
    base_directory::AbstractString
end

type RenderTask <: DaemonTask
    taskId::AbstractString
    indices::AbstractString
    base_directory::AbstractString
end

#=
 =Set the constant identifiers for each task.
 =Use these names to create a lookup table for each constructor of the
 =corresponding task
 =#
const TASK_ID_BLOCK_MATCH = "BLOCK_MATCH"
const TASK_ID_RENDER = "RENDER"
const TASKS = Dict()
TASKS[TASK_ID_BLOCK_MATCH] = BlockMatchTask
TASKS[TASK_ID_RENDER] = RenderTask

#=
 = Given a task in json form, convert it into the correct type
 = Returns DaemonTask
 =#
function parse(message::ASCIIString)
    if length(strip(message))
        error("Trying to parse empty string for task")
    end

    if !haskey(message, "taskId")  || !haskey(message, "indices")
        error("Missing task parameters")
    end

    dictionary = JSON.parse(recieve_response.obj.messageSet.body)
    return TASKS[json.taskId](taskId=dictionary["taskId"],
        task=dictionary["name"], indices=dictionary["indicies"])
end

function execute(block_match_task::BlockMatchTask)
    print("Running BlockMatching with indicies $(block_match_task.indices) for
        $(block_match_task.base_directory)")
end

function execute(render_task::RenderTask)
    print("Running Rendering with indicies $(render_task.indices) for
        $(render_task.base_directory)")
end

end # end module Task
