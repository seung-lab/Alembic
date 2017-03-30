#using Alembic
module MaskTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export MaskTaskDetails, NAME, execute, full_input_path, full_output_path

type MaskTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

MaskTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = MaskTaskDetails(basic_info, AlembicPayloadInfo(dict));


#TODO MAKE THIS PROPER AND NOT ALL OVER THE PLACE
# creates a task to make a MeshSet for index
function MaskTaskDetails(index::Main.Index)
	indices = [index];
	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = map(Main.truncate_path, map(Main.get_registry_path, indices));
	inputs_masks = map(Main.truncate_path, map(Main.get_path, repeated("mask"), indices));


	#inputs = unique(vcat(inputs_images, inputs_registry, inputs_meshes, inputs_meshset, input_transform))
	inputs = unique(vcat(inputs_images, inputs_registry, inputs_masks))
	
#	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
#	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	output_image = Main.truncate_path(Main.get_path(index));
	#output_image_thumbnail = Main.truncate_path(Main.get_path("thumbnail", index));
#	outputs = unique([output_meshset, output_stats, output_transform, output_reviews...])
#	output_tform = 

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	#task = MaskTaskDetails(basic_info, AlembicPayloadInfo([index], vcat([output_image, output_image_thumbnail],inputs_meshset)));
	task = MaskTaskDetails(basic_info, AlembicPayloadInfo([index], [output_image]));
#	return vcat(inputs..., output)
	return task
end

const NAME = "MASK_TASK"

function full_input_path(task::MaskTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::MaskTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::MaskTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::MaskTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    Main.apply_mask(task.payload_info.indices[1]);
    #Main.calculate_stats(ms);

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::MaskTaskDetails,
        datasource::DatasourceService, result::DaemonTask.Result)
    if !result.success
        error("Task $(task.basic_info.id), $(task.basic_info.name) was " *
            "not successful")
    else
        println("Task $(task.basic_info.id), $(task.basic_info.name) was " *
            "completed successfully, syncing outputs to remote datasource")

        Datasource.put!(datasource,
            map((output) -> full_output_path(task, output), result.outputs))
	Datasource.delete!(datasource, map((output) -> full_output_path(task, output), result.outputs); only_cache = true)
	Datasource.delete!(datasource, map((input) -> full_input_path(task, input), task.basic_info.inputs); only_cache = true)
	Main.push_registry_updates();

    return DaemonTask.Result(true, task.payload_info.outputs)
    end
end

end # module MaskTask
