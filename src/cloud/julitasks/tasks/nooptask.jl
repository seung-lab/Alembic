"""
    NoOpTask

This module includes the composite type NoOpDetails which includes both the
        generic DaemonTask.Info and an AbstractString as a payload
"""
module NoOpTask

using ...Julitasks.Types

import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask
import Julimaps.Cloud.Julitasks.Tasks.AlignmentTask
import Julimaps.Cloud.Julitasks.Services.Datasource

export NoOpTaskDetails, name, execute

type NoOpTaskDetails <: DaemonTaskDetails
    basicInfo::BasicTask.Info
    payloadInfo::AbstractString
end

const name = "NO_OP"

function DaemonTask.prepare(task::BlockMatchTaskDetails,
        datasource::DatasourceService)
    #=Datasource.pull!(daemon.datasource, task.basicInfo.files)=#
    println("Preparing NoOpTask")
end

function DaemonTask.execute(task::BlockMatchTaskDetails,
        datasource::DatasourceService)
    println("Executing NoOpTask with payload $(task.payloadInfo)")
end

function DaemonTask.finalize(task::BlockMatchTaskDetails,
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
