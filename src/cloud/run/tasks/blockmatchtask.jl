"""
    BlockMatchTask

This module includes the composite type BlockMatchDetails which includes both
the generic DaemonTask.Info and AlignmentTask.Info.
"""
module BlockMatchTask

import Julimaps.Cloud.Run.Tasks.AlignmentTask
import Julimaps.Cloud.Julitasks.Tasks.DaemonTask
import Julimaps.Cloud.Julitasks.Tasks.BasicTask
import Julimaps.Cloud.Julitasks.Services.Daemon
import Julimaps.Cloud.Julitasks.Services.Datasource

export BlockMatchTaskDetails, name, execute

type BlockMatchTaskDetails <: DaemonTask.Details
    basicInfo::BasicTask.Info
    payloadInfo::AlignmentTask.Info
end

BlockMatchTaskDetails(info::BasicTask.Info, dict::Dict{AbstractString, Any}) =
    BlockMatchTaskDetails(info, AlignmentTask.Info(dict))

const name = "BLOCK_MATCH"

function DaemonTask.prepare(task::BlockMatchTaskDetails,
        datasource::Datasource.Service)
    #=Datasource.pull!(daemon.dataSource, task.basicInfo.files)=#
    println("Executing BlockMatchTask")
end

function DaemonTask.execute(task::BlockMatchTaskDetails,
        datasource::Datasource.Service)
    println("Executing BlockMatchTask")
end

function DaemonTask.finalize(task::BlockMatchTaskDetails,
        datasource::Datasource.Service, result::DaemonTask.Result)
    #=
     =if !result.success
     =    println("Task $(task.details.id), $(task.details.name) was
     =        not successful")
     =else
     =    println("Task $(task.details.id), $(task.details.name) was
     =        completed successfully")
     =    Datasource.push!(daemon.datasource, result.filename)
     =end
     =#
    println("Executing BlockMatchTask")
end

end # module BlockMatchTask
