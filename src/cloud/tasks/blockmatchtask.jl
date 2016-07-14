"""
    BlockMatchTask

This module includes the composite type BlockMatchDetails which includes both
the generic DaemonTask.Info and AlignmentTask.Info.
"""
module BlockMatchTask

import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.BasicTask
import Julimaps.Cloud.Tasks.AlignmentTask

export BlockMatchTaskDetails, name, execute

type BlockMatchTaskDetails <: DaemonTask.Details
    basicInfo::BasicTask.Info
    payloadInfo::AlignmentTask.Info
end

BlockMatchTaskDetails(info::BasicTask.Info, dict::Dict{AbstractString, Any}) =
    BlockMatchTaskDetails(info, AlignmentTask.Info(dict))

const name = "BLOCK_MATCH"

function DaemonTask.execute(task::BlockMatchTaskDetails)
    println("Executing BlockMatchTask")
end

end # module BlockMatchTask
