#using Alembic
module RenderTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export RenderTaskDetails, NAME, execute, full_input_path, full_output_path

type RenderTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

RenderTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = RenderTaskDetails(basic_info, AlembicPayloadInfo(dict));


#TODO MAKE THIS PROPER AND NOT ALL OVER THE PLACE
# creates a task to make a MeshSet for index
function RenderTaskDetails(index::Main.Index)
	indices = Main.ancestors(index);

	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = map(Main.truncate_path, map(Main.get_registry_path, indices));

      if Main.is_prealigned(index) 
	input_transform = [Main.truncate_path(Main.get_path("cumulative_transform", index))]
	inputs_meshes = [];
	inputs_meshset = [];
      elseif Main.is_montaged(index)
#	inputs_meshes = map(Main.truncate_path, map(Main.get_path, repeated("Mesh"), indices));
	inputs_meshes = [];
	inputs_meshset = [Main.truncate_path(Main.get_path("MeshSet", index))];
	input_transform = [];
      elseif Main.is_aligned(index)
	inputs_meshes = map(Main.truncate_path, map(Main.get_path, repeated("Mesh"), indices));
	input_transform = [];
	inputs_meshset = [];
      end

	inputs = unique(vcat(inputs_images, inputs_registry, inputs_meshes, inputs_meshset, input_transform))
	
#	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
#	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	output_image = Main.truncate_path(Main.get_path(index));
#	outputs = unique([output_meshset, output_stats, output_transform, output_reviews...])
#	output_tform = 

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = RenderTaskDetails(basic_info, AlembicPayloadInfo([index], [output_image]));
#	return vcat(inputs..., output)
	return task
end

const NAME = "RENDER_TASK"

function full_input_path(task::RenderTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::RenderTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::RenderTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::RenderTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

#=    if length(task.payload_info.indices) == 2
    ms = Main.MeshSet([tuple(index_array...) for index_array in task.payload_info.indices]...; solve=false);
    split_meshset(ms);
    else
    ms = Main.MeshSet([tuple(index_array...) for index_array in task.payload_info.indices]...);
    Main.calculate_stats(ms);
    end=#
    #actually only a tuple
    Main.render(task.payload_info.indices...);
    #Main.calculate_stats(ms);

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::RenderTaskDetails,
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

end # module RenderTask
