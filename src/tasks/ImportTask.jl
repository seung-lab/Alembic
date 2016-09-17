#using Alembic
module ImportTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export ImportTaskDetails, NAME, execute, full_input_path, full_output_path

type ImportTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

ImportTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = ImportTaskDetails(basic_info, AlembicPayloadInfo(dict));


#TODO MAKE THIS PROPER AND NOT ALL OVER THE PLACE
# creates a task to make a MeshSet for index
function ImportTaskDetails(z_index)
#	inputs = unique(vcat(inputs_images, inputs_registry, inputs_meshes, inputs_meshset, input_transform))

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, [string(z_index)]) 
	task = ImportTaskDetails(basic_info, AlembicPayloadInfo([], []));
#	return vcat(inputs..., output)
	return task
end

const NAME = "IMPORT_TASK"

function full_input_path(task::ImportTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::ImportTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::ImportTaskDetails,
        datasource::DatasourceService)
end

function DaemonTask.execute(task::ImportTaskDetails,
        datasource::DatasourceService)
    z_index = parseint(task.basic_info.inputs[1])

    Main.gentrify_tiles(z_index);
    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::ImportTaskDetails,
        datasource::DatasourceService, result::DaemonTask.Result)
    if !result.success
        error("Task $(task.basic_info.id), $(task.basic_info.name) was " *
            "not successful")
    else
        println("Task $(task.basic_info.id), $(task.basic_info.name) was " *
            "completed successfully, syncing outputs to remote datasource")
	Main.push_registry_updates();
    return DaemonTask.Result(true, task.payload_info.outputs)
    end
end

end # module ImportTask
