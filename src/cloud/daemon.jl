module Daemon

import Julimaps.Cloud.Queue
import Julimaps.Cloud.Bucket
import Julimaps.Cloud.DaemonTask

export DaemonService
export run

type DaemonService
    queue::Queue.QueueService
    bucket::Bucket.BucketService
    poll_frequency_seconds::Int64
    tasks::Dict{AbstractString, Module}
end

function run(daemon::DaemonService)
    while true
        try
            message = Queue.pop_message(daemon.queue)

            if isempty(message)
                println("No messages found in $(Queue.string(daemon.queue))")
            else
                #=task = DaemonTask.parse(message)=#
                task = parse(daemon, message)

                println("Task is $(task.taskId)")

                #=DaemonTask.execute(task)=#
            end
        catch e
            showerror(STDERR, e, catch_backtrace(); backtrace = true)
            println(STDERR) #looks like showerror doesn't include a newline
        end

        sleep(daemon.poll_frequency_seconds)
    end
end

function register(daemon::DaemonService, task_module::Module)
    symbols = names(task_module, true)
    # Much rather use in, but doesn't seem to work with array of symbols
    if findfirst(symbols, :execute) <= 0
        error("Module $task_module does not contain execute function")
    end
    if findfirst(symbols, :task_type) <= 0
        error("Module $task_module does not contain a task type function")
    end
    if !task_module.task_type <: DaemonTask.DaemonTaskDetails
        error("Module $task_module task type does not subtype
            DaemonTaskDetails")
    end
    daemon.tasks[task_module.name] = task_module
end

#=
 = Given a task in json form, convert it into the correct type
 = Returns DaemonTask
 =#
function parse(daemon::DaemonService, message::ASCIIString)
    message = strip(message)
    if isempty(message)
        error("Trying to parse empty string for task")
    end

    json = JSON.parse(message)

    if !haskey(dictionary, "details")
        error("Could not find task details from parsing message ($message)")
    end

    return to_daemon_task(json, daemon.tasks)
end


function register(daemon::DaemonService, task_name::AbstractString,
    task_type::DataType)
    if !task_type <: DaemonTask.DaemonTaskDetails
        throw TypeError("Trying to register task that is not a
            DaemonTaskDetail, got $task_name with $task_type")
    end
    daemon.tasks[task_name] = task_type
end

end # end module Daemon
