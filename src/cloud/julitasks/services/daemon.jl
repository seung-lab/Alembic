module Daemon

using ...Julitasks.Types

import Julimaps.Cloud.Julitasks.Services.Queue
import Julimaps.Cloud.Julitasks.Services.Bucket
import Julimaps.Cloud.Julitasks.Services.Datasource
import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask
import JSON

export DaemonService, register, run

type DaemonService
    queue::QueueService
    bucket::BucketService
    datasource::DatasourceService
    poll_frequency_seconds::Int64
    tasks::Dict{AbstractString, Type}
end

DaemonService(queue::QueueService, bucket::BucketService,
    datasource::DatasourceService, poll_frequency_seconds::Int64) = 
        DaemonService(queue, bucket, datasource, poll_frequency_seconds,
            Dict{AbstractString, Module}())

function run(daemon::DaemonService)
    while true
        try
            message = Queue.pop_message(daemon.queue)

            print("Message received is$(message)")

            if isempty(message)
                println("No messages found in $(Queue.string(daemon.queue))")
            else
                task = parse(message)

                println("Task is $(task.details.id), $(task.details.name)")

                success = DaemonTask.run(task, daemon.datasource)
            end
        catch e
            showerror(STDERR, e, catch_backtrace(); backtrace = true)
            println(STDERR) #looks like showerror doesn't include a newline
        end

        sleep(daemon.poll_frequency_seconds)
    end
end

"""
    register!(daemon::DaemonService, task_module::Module)

Register the task type generation for the given task_name

"""
function register!(daemon::DaemonService, task_name::AbstractString,
        task_type::Type)
    if !DaemonTask.can_execute(task_type)
        error("Can not register $task_type with name $task_name " *
            " in Daemon. Could not find a registered method to execute")
    end

    if haskey(daemon.tasks, task_name)
        warn("Daemon has already registered $task_name with " *
            "$(daemon.tasks[task_name])")
    end

    daemon.tasks[task_name] = task_type
end

"""
    parse(daemon::DaemonService, text::ASCIIString)

Parse input JSON message into a task object
"""
function parse(daemon::DaemonService, text::ASCIIString)
    text = strip(text)

    if isempty(text)
        throw(ArgumentError("Trying to parse empty string for task"))
    end

    message = JSON.parse(text)

    basic_info = BasicTask.Info(message["basicInfo"])

    if !haskey(daemon.tasks, basic_info.name)
        error("Task $(basic_info.name) is not registered with the daemon")
    end

    task_type = daemon.tasks[basic_info.name]
    return task_type(basic_info, message["payloadInfo"])
end

end # end module Daemon
