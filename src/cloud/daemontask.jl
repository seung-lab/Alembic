module DaemonTask

import JSON

export DaemonTaskDetails
export BlockMatchTask
export RenderTask
export parse_task

# These types identify the task that we like to perform
abstract DaemonTaskDetails

type Details
    taskType::AbstractString
    baseDirectory::AbstractString
    files::Array{AbstractString}
    indices::Array{Tuple{Int64,Int64,Int64,Int64}}
end

type BlockMatchTask <: DaemonTaskDetails
    details::Details
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

    return to_daemon_task(json)
end

function execute(block_match_task_details::BlockMatchTask)
    print("Running BlockMatching with indicies
        $(block_match_task_details.indices) for 
        $(block_match_task_details.base_directory)")
end

function execute(render_task_details::RenderTask)
    print("Running Rendering with indicies $(render_task_details.indices) for
        $(render_task_details.base_directory)")
end

function to_daemon_task(dictionary::Dict)
    if !haskey(dictionary, "details")
        error("Missing details")
    end

    details = dictionary["details"]

    if !haskey(TASKS, details["taskType"])
        error("Unknown task : $(details["taskType"])")
    end

    indices = []
    for index in details["indices"]
        push!(indices, (index[1],index[2],index[3],index[4]))
    end

    return TASKS[details["taskType"]](
        Details(
            details["taskType"],
            details["baseDirectory"],
            details["files"],
            indices
        )
    )
end

end # module Task
