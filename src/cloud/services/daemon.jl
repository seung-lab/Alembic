module Daemon

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Datasource
import Julimaps.Cloud.Tasks.DaemonTask
import JSON

export Service, register, run

type Service
    queue::Queue.Service
    bucket::Bucket.Service
    dataSource::Datasource.Service
    poll_frequency_seconds::Int64
    tasks::Dict{AbstractString, Type}
    Service(queue::Queue.Service, bucket::Bucket.Service,
        datasource::Datasource.Service, poll_frequency_seconds::Int64) = 
            new(queue, bucket, datasource, poll_frequency_seconds,
                Dict{AbstractString, Module}())
end

function run(daemon::Service)
    while true
        try
            message = Queue.pop_message(daemon.queue)

            print("Message received is$(message)")

            if isempty(message)
                println("No messages found in $(Queue.string(daemon.queue))")
            else
                task = parse(message)

                println("Task is $(task.details.id), $(task.details.name)")

                prepare(daemon, task)

                result = DaemonTask.execute(task)

                finalize(daemon, task, result)

            end
        catch e
            showerror(STDERR, e, catch_backtrace(); backtrace = true)
            println(STDERR) #looks like showerror doesn't include a newline
        end

        sleep(daemon.poll_frequency_seconds)
    end
end

"""
    register(daemon::Service, task_module::Module)

Register the task type generation for the given task_name

"""
function register!(daemon::Service, task_name::AbstractString,
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

function parse(daemon::Service, text::ASCIIString)
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

function prepare(daemon::Service, task::DaemonTask.Details)
    for filename in task.basicInfo.files
        Datasource.pull!(daemon.dataSource, filename)
    end
end

function finalize(daemon::Service, task::DaemonTask.Details,
    result::DaemonTask.Result)
    if !result.success
        println("Task $(task.details.id), $(task.details.name) was
            not successful")
    else
        println("Task $(task.details.id), $(task.details.name) was
            completed successfully")
        Datasource.push!(daemon.datasource, result.filename)
    end
end

end # end module Daemon
