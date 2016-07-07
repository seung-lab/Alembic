module DaemonTask

import JSON

export DaemonTaskDetails, Details, execute

"""
    DaemonTask

This is the base composite abstract class used to compose Details and Payload
i.e. compose a task with
```julia
type YourDaemonTask <: DaemonTaskDetails
    details::Details
    payload::YourTaskDetails
end
```
"""
abstract DaemonTaskDetails

"""
    DaemonTask.Details

Type contains basic information for a daemon task

"""
type Details
    id::Int64
    name::AbstractString
end

function Details(dict::Dict{AbstractString, Any}) 
    id = typeof(dict["id"]) <: Int ?  dict["id"] : parse(Int64, dict["id"])
    if isempty(strip(dict["name"]))
        throw(ArgumentError("Task name can not be empty"))
    end
    return Details(id, dict["name"])
end

"""
    execute(task::DaemonTask)

Executes the given task. Must be overriden for new tasks
"""
function execute(task::DaemonTaskDetails)
    error("Execute is unimplemented for this task $task")
end

end # module DaemonTask
