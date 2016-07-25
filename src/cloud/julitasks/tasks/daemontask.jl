module DaemonTask

using ...Julitasks.Types

export Info, run

"""
    DaemonTask.Result

This object contains the outcome of performing the task.
"""
type Result
    success::Bool
    outputs::Array{AbstractString, 1}
end

# Test to see if the execute function exists for this type
function can_execute(task_type::Type)
    if !(task_type <: DaemonTaskDetails)
        return false
    end

    execute_methods = methods(execute, Any[task_type, DatasourceService])

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
    run(task::DaemonTaskDetails, datasource::DatasourceService)

Run the current task. 1. Prepare, 2. Execute 3. Finalize
"""
function run(task::DaemonTaskDetails, datasource::DatasourceService)
    prepare(task, datasource)
    result = execute(task, datasource)
    finalize(task, datasource, result)
end

"""
    prepare(task::DaemonTaskDetails, datasource::DatasourceService)

prepare what is needed for the task
"""
function prepare(task::DaemonTaskDetails, datasource::DatasourceService)
    error("Prepare is unimplemented for this task $task")
end

"""
    execute(task::DaemonTaskDetails, datasource::DatasourceService)

Executes the given task. Must be overriden for new tasks
"""
function execute(task::DaemonTaskDetails, datasource::DatasourceService)
    error("Execute is unimplemented for this task $task")
end

"""
    finalize(daemon::DaemonService, task::DaemonTaskDetails,

After task has completed, perform this action
"""
function finalize(task::DaemonTaskDetails, datasource::DatasourceService)
    error("finalize is unimplemented for this task $task")
end

end # module DaemonTask
