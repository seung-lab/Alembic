module Daemon

import Julimaps.Cloud.Services.Queue
import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Datasource
import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BasicTask
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

function prepare(daemon::Service, task::DaemonTask.Details)
    Datasource.pull!(daemon.dataSource, task.basicInfo.files)
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
