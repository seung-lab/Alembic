#using Alembic
module SaveStackTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export SaveStackTaskDetails, NAME, execute, full_input_path, full_output_path

type SaveStackTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

SaveStackTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = SaveStackTaskDetails(basic_info, AlembicPayloadInfo(dict));


#TODO MAKE THIS PROPER AND NOT ALL OVER THE PLACE
# creates a task to make a MeshSet for index
function SaveStackTaskDetails(index::Main.Index)

	inputs_images = [Main.truncate_path(Main.get_path(index))];
	inputs_registry = [Main.truncate_path(Main.get_registry_path(index))];

	inputs = unique(vcat(inputs_images, inputs_registry))
	
#	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
#	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	# output_image = Main.truncate_path(Main.get_path(index));
	output_image = Main.truncate_path(Main.get_path(finished(index)));
#	outputs = unique([output_meshset, output_stats, output_transform, output_reviews...])
#	output_tform = 

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = SaveStackTaskDetails(basic_info, AlembicPayloadInfo([index], [output_image]));
#	return vcat(inputs..., output)
	return task
end

const NAME = "SAVE_STACK_TASK"

function full_input_path(task::SaveStackTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::SaveStackTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::SaveStackTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::SaveStackTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end
    
    slice = (45000:54999, 23000:23049)
    origin = [0,0]
    x_slice = [slice[1][1], slice[1][end]] + origin
    y_slice = [slice[2][1], slice[2][end]] + origin
    index = task.payload_info.indices[1];
    img = Main.make_stack(index, index, slice)
    f = Main.h5open(Main.get_path(Main.finished(index)), "w")
    chunksize = min(1000, min(size(img)...))
    f["img", "chunk", (chunksize, chunksize)] = img
    f["origin"] = origin
    f["x_slice"] = x_slice
    f["y_slice"] = y_slice
    f["z_slice"] = z_slice
    close(f)

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::SaveStackTaskDetails,
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

end # module SaveStackTask
