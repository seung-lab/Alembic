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


# creates a task to make a MeshSet for index
function BlockMatchTaskDetails(index::Main.Index)
  #=
	if Main.is_montaged(index)	indices = Main.get_index_range(Main.prevstage(index),Main.prevstage(index))
	elseif Main.is_prealigned(index) indices = [Main.prevstage(index), Main.get_preceding(Main.prevstage(index))]
	end
	=#
	indices = Main.ancestors(index);
	possible_pairs = collect([(indexA, indexB) for indexA in indices, indexB in indices])
	#note the flipped order for the prealignment
	possible_pairs = possible_pairs[map(pair -> (Main.is_adjacent(pair[1], pair[2])) || (Main.is_diagonal(pair[1], pair[2])) || (Main.is_preceding(pair[2], pair[1])), possible_pairs)]

	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = map(Main.truncate_path, map(Main.get_registry_path, indices));
	inputs = unique(vcat(inputs_images, inputs_registry))
	
	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
	output_stats = Main.truncate_path(Main.get_path("stats", index))
	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	outputs = unique([output_meshset, output_stats, output_transform, output_reviews...])
#	output_tform = 

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = BlockMatchTaskDetails(basic_info, AlembicPayloadInfo([index], outputs));
#	return vcat(inputs..., output)
	return task
end

function BlockMatchTaskDetails(first_index::Main.Index, last_index::Main.Index)
  	#= indices = Main.get_index_range(first_index, last_index);
	possible_pairs = collect([(indexA, indexB) for indexA in indices, indexB in indices])
	possible_pairs = possible_pairs[map(pair -> Main.is_preceding(pair[1], pair[2], Main.get_params(last_index)["match"]["depth"]), possible_pairs)] =#

	possible_pairs = [(first_index, last_index)]
	indices = Main.get_index_range(first_index, last_index);

	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = map(Main.truncate_path, map(Main.get_registry_path, indices));
	inputs = unique(vcat(inputs_images, inputs_registry))
	
	output_meshset = Main.truncate_path(Main.get_path("MeshSet", (first_index, last_index)))
	output_meshes = []
	#output_meshes = map(Main.truncate_path, map(Main.get_path, repeated("Mesh"), indices))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	output_matches = map(Main.truncate_path, map((pair) -> Main.get_path("Match", pair), possible_pairs))
	outputs = unique([output_meshset, output_meshes..., output_matches..., output_reviews...])

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = BlockMatchTaskDetails(basic_info, AlembicPayloadInfo([first_index, last_index], outputs));
	return task
end

const NAME = "BLOCKMATCH_TASK"

function full_input_path(task::BlockMatchTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
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

    # alignment case
    if length(task.payload_info.indices) == 2
    ms = Main.MeshSet(task.payload_info.indices...; solve=false);
    Main.split_meshset(ms);
    # all other case
    else
    ms = Main.MeshSet(task.payload_info.indices...);
    Main.calculate_stats(ms);
    #Main.render(ms; review=true);
    end
    Main.render(ms; review=true);
    #Main.render(ms; review=true);
    #Main.calculate_stats(ms);

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
	Datasource.delete!(datasource, map((output) -> full_output_path(task, output), result.outputs); only_cache = true)
	Datasource.delete!(datasource, map((input) -> full_input_path(task, input), task.basic_info.inputs); only_cache = true)
	Main.push_registry_updates();

    return DaemonTask.Result(true, task.payload_info.outputs)
    end
end

end # module BlockMatchTask
