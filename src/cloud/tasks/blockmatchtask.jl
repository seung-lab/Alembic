"""
    BlockMatchTask

This module includes the composite type BlockMatchDetails which includes both
the generic DaemonTask.Details and AlignmentTask.Details.
"""
module BlockMatchTask

import Julimaps.Cloud.Tasks.DaemonTask
import Julimaps.Cloud.Tasks.AlignmentTask

export BlockMatchTask, name, execute

type BlockMatchTaskDetails <: DaemonTask.DaemonTaskDetails
    details::DaemonTask.Details
    payload::AlignmentTask.Details
    BlockMatchTask(details::DaemonTask.Details,
        dict::Dict{AbstractString, Any}) =
            new(details, AlignmentTask.Details(dict))
end

const name = "BLOCK_MATCH"

function DaemonTask.execute(task::BlockMatchTaskDetails)
    println("BlockMatchTask")
end

end # module BlockMatchTask
