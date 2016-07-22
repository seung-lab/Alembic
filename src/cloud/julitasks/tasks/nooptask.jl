"""
    NoOpTask

This module includes the composite type NoOpDetails which includes both the
        generic DaemonTask.Info and an AbstractString as a payload
"""
module NoOpTask

using ...Julitasks.Types

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask

export NoOpTaskDetails, name, execute

type NoOpTaskDetails <: DaemonTaskDetails
    basicInfo::BasicTask.Info
    payloadInfo::AbstractString
end

const name = "NO_OP"

function DaemonTask.prepare(task::NoOpTaskDetails,
        datasource::DatasourceService)
    Datasource.pull!(daemon.datasource, task.basicInfo.files)
    println("Preparing NoOpTask")
end

function DaemonTask.execute(task::NoOpTaskDetails,
        datasource::DatasourceService)
    file = datasource.get(task.basicInfo.files[1])
    println(readall(file))
    return DaemonTask.Result(true, "output")
end

function DaemonTask.finalize(task::NoOpTaskDetails,
        datasource::DatasourceService, result::DaemonTask.Result)
    if !result.success
        println("Task $(task.details.id), $(task.details.name) was
            not successful")
    else
        println("Task $(task.details.id), $(task.details.name) was
            completed successfully")
        Datasource.push!(daemon.datasource, result.filename)
    end
end

end # module BlockMatchTask
