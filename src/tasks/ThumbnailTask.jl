#using Alembic
module ThumbnailTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export ThumbnailTaskDetails, NAME, execute, full_input_path, full_output_path

type ThumbnailTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

ThumbnailTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = ThumbnailTaskDetails(basic_info, AlembicPayloadInfo(dict));


#TODO MAKE THIS PROPER AND NOT ALL OVER THE PLACE
# creates a task to make a MeshSet for index
function ThumbnailTaskDetails(index::Main.Index)

	inputs_images = [Main.truncate_path(Main.get_path(index))];
	inputs_registry = [Main.truncate_path(Main.get_registry_path(index))];

	inputs = unique(vcat(inputs_images, inputs_registry))
	
#	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
#	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	# output_image = Main.truncate_path(Main.get_path(index));
	output_image_thumbnail = Main.truncate_path(Main.get_path("thumbnail", index));
#	outputs = unique([output_meshset, output_stats, output_transform, output_reviews...])
#	output_tform = 

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = ThumbnailTaskDetails(basic_info, AlembicPayloadInfo([index], [output_image_thumbnail]));
#	return vcat(inputs..., output)
	return task
end

const NAME = "THUMBNAIL_TASK"

function full_input_path(task::ThumbnailTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::ThumbnailTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::ThumbnailTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::ThumbnailTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    dir_path = joinpath(Main.ALIGNED_DIR_PATH, Main.THUMBNAIL_DIR)
    if !isdir(dir_path)    
        mkdir(dir_path)
    end
    
    thumbnail_scale = 0.02
    index = task.payload_info.indices[1];
    img = Main.load(index)
    thumbnail, _ = Main.imscale(img, thumbnail_scale)
    Main.write_thumbnail(thumbnail, index, thumbnail_scale)

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::ThumbnailTaskDetails,
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

end # module ThumbnailTask
