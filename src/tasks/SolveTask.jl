#using Alembic
module SolveTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export SolveTaskDetails, NAME, execute, full_input_path, full_output_path

type SolveTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

SolveTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = SolveTaskDetails(basic_info, AlembicPayloadInfo(dict));

#=
# creates a task to make a MeshSet for index
function SolveTaskDetails(index::Main.Index)
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
	task = SolveTaskDetails(basic_info, AlembicPayloadInfo([index], outputs));
#	return vcat(inputs..., output)
	return task
end
=#

function SolveTaskDetails(first_index::Main.Index, last_index::Main.Index)
	indices = Main.get_index_range(first_index, last_index);
	possible_pairs = collect([(indexA, indexB) for indexA in indices, indexB in indices])
	possible_pairs = possible_pairs[map(pair -> Main.is_preceding(pair[1], pair[2]), possible_pairs)]

#	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = map(Main.truncate_path, map(Main.get_registry_path, indices));
	input_meshes = map(Main.truncate_path, map(Main.get_path, repeated("Mesh"), indices))
	input_matches = map(Main.truncate_path, map((pair) -> Main.get_path("Match", pair), possible_pairs))
	inputs = unique(vcat(input_meshes, input_matches, inputs_registry))
	
#	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
	output_meshes = map(Main.truncate_path, map(Main.get_path, repeated("Mesh"), indices))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
#	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	output_matches = map(Main.truncate_path, map((pair) -> Main.get_path("Match", pair), possible_pairs))
	outputs = unique([output_meshes..., output_matches...])

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = SolveTaskDetails(basic_info, AlembicPayloadInfo([first_index, last_index], outputs));
	return task
end

const NAME = "SOLVE_TASK"

function full_input_path(task::SolveTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::SolveTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::SolveTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::SolveTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    ms = Main.compile_meshset([tuple(index_array...) for index_array in task.payload_info.indices]...);
    #temporary
    Main.clear_filters!(ms);
    ms.properties["params"] = Main.PARAMS_ALIGNMENT
    for match in ms.matches
    match.properties["params"] = Main.PARAMS_ALIGNMENT
  end
    for mesh in ms.meshes
    mesh.properties["params"] = Main.PARAMS_ALIGNMENT
  end
  Main.filter!(ms);


    Main.fix_ends!(ms)
    Main.solve!(ms);
    Main.split_meshset(ms);

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::SolveTaskDetails,
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

end # module SolveTask
