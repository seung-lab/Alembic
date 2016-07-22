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

function DaemonTask.prepare(task::NoOpTaskDetails,
        datasource::DatasourceService)
    println("Preparing NoOpTask")
    Datasource.pull!(datasource, task.basicInfo.files)
end

function DaemonTask.execute(task::NoOpTaskDetails,
        datasource::DatasourceService)
    println("executing task $task, files contain")
    file = Datasource.pull!(datasource, task.basicInfo.files[1])
    if file != nothing
        println(readall(file))
    end
    return DaemonTask.Result(true, ["output"])
end

function DaemonTask.finalize(task::NoOpTaskDetails,
        datasource::DatasourceService, result::DaemonTask.Result)
    if !result.success
        println("Task $(task.basicInfo.id), $(task.basicInfo.name) was
            not successful")
    else
        println("Task $(task.basicInfo.id), $(task.basicInfo.name) was
            completed successfully")
        Datasource.push!(datasource, result.output)
    end
end

end # module BlockMatchTask
