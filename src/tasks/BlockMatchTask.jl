#using Alembic

module BlockMatchTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export BlockMatchTaskDetails, NAME, execute, full_input_path, full_output_path

type BlockMatchTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

const NAME = "BLOCKMATCH_TASK"
#const OUTPUT_FOLDER = "output"

# truncates path
function truncate_path(path::AbstractString)
	m = Base.match(Regex("$(Main.TASKS_BASE_DIRECTORY)/(\\S+)"), path);
	return m[1]
end


function full_input_path(task::BlockMatchTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$input"
end

function full_output_path(task::BlockMatchTaskDetails,
        input::AbstractString)
#    path_end = rsearch(input, "/").start + 1

    return "$(task.basic_info.base_directory)/$(task.payload_info)";
end

function DaemonTask.prepare(task::BlockMatchTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs))
end

function DaemonTask.execute(task::BlockMatchTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    Main.MeshSet([tuple(index_array...) for index_array in task.payload_info.indices]...);

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::BlockMatchTaskDetails,
        datasource::DatasourceService, result::DaemonTask.Result)
    if !result.success
        error("Task $(task.basic_info.id), $(task.basic_info.name) was " *
            "not successful")
    else
        println("Task $(task.basic_info.id), $(task.basic_info.name) was " *
            "completed successfully, syncing outputs to remote datasource")
        Datasource.put!(datasource,
            map((output) -> full_output_path(task, output), result.outputs))
	Datasource.remove!(datasource, map((output) -> full_output_path(task, output), result.outputs; only_cache = true))
	Datasource.remove!(datasource, map((input) -> full_input_path(task, input), task.basic_info.inputs; only_cache = true))
    #	task_queue = AWSQueueService(AWS.AWSEnv(), registry_queue_name);
    end
end

end # module BlockMatchTask
