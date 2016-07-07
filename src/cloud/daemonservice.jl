module Daemon

import Julimaps.Cloud.Queues.Queue
import Julimaps.Cloud.Buckets.Bucket
import Julimaps.Cloud.Tasks.DaemonTask
import JSON

export DaemonService
export run

type DaemonService
    queue::Queue.QueueService
    bucket::Bucket.BucketService
    poll_frequency_seconds::Int64
    tasks::Dict{AbstractString, Type}
    DaemonService(queue::Queue.QueueService, bucket::Bucket.BucketService,
        poll_frequency_seconds::Int64) = 
            new(queue, bucket, poll_frequency_seconds,
                Dict{AbstractString, Module}())
end

function run(daemon::DaemonService)
    while true
        try
            message = Queue.pop_message(daemon.queue)

            print("Message received is$(message)")

            if isempty(message)
                println("No messages found in $(Queue.string(daemon.queue))")
            else
                task = DaemonTask.parse(message)

                println("Task is $(task.details.id), $(task.details.name)")

                success = DaemonTask.execute(task)

                if !success
                    println("Task $(task.details.id), $(task.details.name) was
                    not successful")
                end
            end
        catch e
            showerror(STDERR, e, catch_backtrace(); backtrace = true)
            println(STDERR) #looks like showerror doesn't include a newline
        end

        sleep(daemon.poll_frequency_seconds)
    end
end

"""
    register(daemon::DaemonService, task_module::Module)

Register the task type generation for the given task_name

"""
function register!(daemon::DaemonService, task_name::AbstractString, task_type::Type)
    if length(methods(DaemonTask.execute, Any[task_type])) == 0
        error("Can not register $task_type to execute in DaemonTask!")
    end

    if haskey(daemon.tasks, task_name)
        warn("Daemon has already registered $task_name with
        $(daemon.tasks[task_name])")
    end

    daemon.tasks[task_name] = task_type
end

function parse(daemon::DaemonService, text::ASCIIString)
    text = strip(text)

    if isempty(text)
        error("Trying to parse empty string for task")
    end

    message = JSON.parse(text)

    if !haskey(message, "details")
        error("Could not find task details from parsing message ($text)")
    end

    details = DaemonTask.Details(message["details"])

    if !haskey(message, "payload")
        error("Could not find task payload from parsing message ($text)")
    end

    if !haskey(daemon.tasks, details.name)
        error("Task $(details.name) is not registered with the daemon")
    end

    return daemon.tasks[task_name](details, message["payload"])
end

end # end module Daemon
