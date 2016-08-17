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

BlockMatchTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = BlockMatchTaskDetails(basic_info, AlembicPayloadInfo(dict));

const NAME = "BLOCKMATCH_TASK"

function full_input_path(task::BlockMatchTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$input"
end

function full_output_path(task::BlockMatchTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::BlockMatchTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::BlockMatchTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    ms = Main.MeshSet([tuple(index_array...) for index_array in task.payload_info.indices]...);
    Main.calculate_stats(ms);

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
	Datasource.remove!(datasource, map((output) -> full_output_path(task, output), result.outputs); only_cache = true)
	Datasource.remove!(datasource, map((input) -> full_input_path(task, input), task.basic_info.inputs); only_cache = true)
	Main.push_registry_updates();

    return DaemonTask.Result(true, task.payload_info.outputs)

    end
end

end # module BlockMatchTask
