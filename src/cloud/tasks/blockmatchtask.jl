module BlockMatchTask
import DaemonTask

export BlockMatchTask
export name, execute, task_type

name = "BLOCK_MATCH"

type AlignmentDetails
    baseDirectory::AbstractString
    files::Array{AbstractString}
    indices::Array{Tuple{Int64,Int64,Int64,Int64}}
end

typealias task_type BlockMatchTask

type BlockMatchTask <: DaemonTask.DaemonTaskDetails
    details::DaemonTask.TaskDetails
    payload::AlignmentDetails
end

function to_daemon_task(dictionary::Dict, tasks::Dict{AbstractString, Module})
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

end # module BlockMatchTask
