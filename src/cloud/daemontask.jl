module DaemonTask

import JSON

export DaemonTaskDetails
export BlockMatchTask
export RenderTask
export parse_task

# These types identify the task that we like to perform
abstract DaemonTaskDetails

type BlockMatchTaskDetails <: DaemonTaskDetails
    taskId::AbstractString
    indices::AbstractString
    base_directory::AbstractString
end

type RenderTaskDetails <: DaemonTaskDetails
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
TASKS[TASK_ID_BLOCK_MATCH] = BlockMatchTaskDetails
TASKS[TASK_ID_RENDER] = RenderTaskDetails

#=
 = Given a task in json form, convert it into the correct type
 = Returns DaemonTask
 =#
function parse(message::ASCIIString)
    message = strip(message)
    if isempty(message)
        error("Trying to parse empty string for task")
    end

    json = JSON.parse(message)

    if !haskey(json, "taskId")  || !haskey(json, "indices")
        error("Missing task parameters")
    end

    if !haskey(TASKS, json["taskId"])
        error("Unknown task : $(json["taskId"])")
    end

    return TASKS[json["taskId"]](
        json["taskId"],
        json["name"],
        json["indices"]
    )
end

function execute(block_match_task_details::BlockMatchTaskDetails)
    print("Running BlockMatching with indicies
        $(block_match_task_details.indices) for 
        $(block_match_task_details.base_directory)")
end

function execute(render_task_details::RenderTaskDetails)
    print("Running Rendering with indicies $(render_task_details.indices) for
        $(render_task_details.base_directory)")
end

end # end module Task
