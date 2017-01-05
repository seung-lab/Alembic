#using Alembic
module RenderReviewTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export RenderReviewTaskDetails, NAME, execute, full_input_path, full_output_path

type RenderReviewTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

RenderReviewTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = RenderReviewTaskDetails(basic_info, AlembicPayloadInfo(dict));


#TODO MAKE THIS PROPER AND NOT ALL OVER THE PLACE
# creates a task to make a MeshSet for index
function RenderReviewTaskDetails(src_index::Main.Index, dst_index::Main.Index)

    indices = [src_index, dst_index]
	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = [Main.truncate_path(Main.get_registry_path(src_index))];

	inputs = unique(vcat(inputs_images, inputs_registry))
	
#	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
#	output_stats = Main.truncate_path(Main.get_path("stats", index))
#	output_transform = Main.truncate_path(Main.get_path("relative_transform", index))
#	output_reviews = map(Main.truncate_path, map((pair) -> Main.get_path("review", pair), possible_pairs))
	# output_image = Main.truncate_path(Main.get_path(index));
	output_image_thumbnail = Main.truncate_path(Main.get_path("review", (indices...)));
#	outputs = unique([output_meshset, output_stats, output_transform, output_reviews...])
#	output_tform = 

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = RenderReviewTaskDetails(basic_info, AlembicPayloadInfo(indices, [output_image_thumbnail]));
#	return vcat(inputs..., output)
	return task
end

const NAME = "RENDER_REVIEW_TASK"

function full_input_path(task::RenderReviewTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::RenderReviewTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::RenderReviewTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
	Main.reload_registries();
end

function DaemonTask.execute(task::RenderReviewTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    s = 0.10
    tform = Main.make_scale_matrix(s)

    src_index = task.payload_info.indices[1];
    dst_index = task.payload_info.indices[2];

    src_img = Main.get_image(src_index)
    dst_img = Main.get_image(dst_index)
    src_offset = Main.get_offset(src_index)
    dst_offset = Main.get_offset(dst_index)
    src_img, src_offset = Main.imwarp(src_img, tform, src_offset)
    dst_img, dst_offset = Main.imwarp(dst_img, tform, dst_offset)
    
    path = Main.get_path("review", (src_index, dst_index))
    Main.write_review_image(path, src_img, src_offset, dst_img, dst_offset, s, tform)

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::RenderReviewTaskDetails,
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

end # module RenderReviewTask
