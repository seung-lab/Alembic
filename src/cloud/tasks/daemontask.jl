module DaemonTask

import JSON

export Details, Info, execute

"""
    DaemonTask.Details

This is the base composite abstract class used to compose Details and Payload
i.e. compose a task with
```julia
type YourDaemonTaskDetails <: DaemonTask.Details
    basicInfo::BasicTask.Info
    taskInfo::YourTaskDetails
end
```
"""
abstract Details

"""
    DaemonTask.Result

This object contains the outcome of performing the task.
"""
type Result
    success::Bool
    filename::AbstractString
end

# Test to see if the execute function exists for this type
function can_execute(task_type::Type)
    if !(task_type <: Details)
        return false
    end

    execute_methods = methods(DaemonTask.execute, Any[task_type])
    if length(execute_methods) == 0
        return false
    end
    for execute_method in execute_methods
        sig_types = execute_method.sig.types
        if  length(sig_types) > 0 && sig_types[1] == task_type
            return true
        end
    end
    return false
end

"""
    execute(task::DaemonTask)

Executes the given task. Must be overriden for new tasks
"""
function execute(task::Details)
    error("Execute is unimplemented for this task $task")
end

end # module DaemonTask
