"""
    NoOpTask

This module includes the composite type NoOpDetails which includes both the
        generic DaemonTask.Info and an AbstractString as a payload
"""
module NoOpTask

using ...Julitasks.Types

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask
import Julimaps.Cloud.Julitasks.Services.Datasource

export NoOpTaskDetails, name, execute

type NoOpTaskDetails <: DaemonTaskDetails
    basicInfo::BasicTask.Info
    payloadInfo::AbstractString
end

const name = "NO_OP"

function full_input_path{String <: AbstractString}(task::NoOpTaskDetails,
        keys::Array{String, 1})
    return map((key) -> full_input_path(task, key), keys)
end

function full_input_path(task::NoOpTaskDetails,
        key::AbstractString)
    return "$(task.basicInfo.baseDirectory)/$key"
end

function DaemonTask.prepare(task::NoOpTaskDetails,
        datasource::DatasourceService)
    println("Preparing NoOpTask")
    Datasource.get(datasource, full_input_path(task, task.basicInfo.inputs))
end

function DaemonTask.execute(task::NoOpTaskDetails,
        datasource::DatasourceService)
    println("executing task NOOP $(task.basicInfo.id), inputs contain")
    for input in task.basicInfo.inputs
        data_stream = Datasource.get(datasource, full_input_path(task, input))
        if data_stream != nothing
            data = readall(data_stream)
            println("Input: $input contains $(data[1:min(20, end)]) \nto\n" *
                "$(data[max(1, end-20):end])")
        end
    end
    return DaemonTask.Result(true, ["output1"])
end

function DaemonTask.finalize(task::NoOpTaskDetails,
        datasource::DatasourceService, result::DaemonTask.Result)
    if !result.success
        println("Task $(task.basicInfo.id), $(task.basicInfo.name) was " *
            "not successful")
    else
        println("Task $(task.basicInfo.id), $(task.basicInfo.name) was " *
            "completed successfully")
        Datasource.put!(datasource, result.outputs)
    end
end

end # module BlockMatchTask
