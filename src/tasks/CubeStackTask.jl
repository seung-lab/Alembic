#using Alembic
module CubeStackTask

using ...SimpleTasks.Types

import Main.AlembicPayloadInfo
import SimpleTasks.Tasks.DaemonTask
import SimpleTasks.Tasks.BasicTask
import SimpleTasks.Services.Datasource

export CubeStackTaskDetails, NAME, execute, full_input_path, full_output_path

type CubeStackTaskDetails <: DaemonTaskDetails
    basic_info::BasicTask.Info
    payload_info::AlembicPayloadInfo
end

CubeStackTaskDetails{String <: AbstractString}(basic_info::BasicTask.Info, dict::Dict{String, Any}) = CubeStackTaskDetails(basic_info, AlembicPayloadInfo(dict));


function CubeStackTaskDetails(index::Main.Index)

    inputs_images = [Main.truncate_path(Main.get_path(index))];
    inputs_registry = [Main.truncate_path(Main.get_registry_path(index))];

    inputs = unique(vcat(inputs_images, inputs_registry, cube_origin, cube_dims, overlap, grid_dims))
    
    # NO OUTPUT IMAGE IN THE PAYLOAD - too many files
    # The execute function will run a Base.run cp or rsync for all
    # output_image = Main.truncate_path(Main.get_path(Main.finished(index)));

    basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
    task = CubeStackTaskDetails(basic_info, AlembicPayloadInfo([index], []));
#   return vcat(inputs..., output)
    return task
end

const NAME = "CUBE_STACK_TASK"

function full_input_path(task::CubeStackTaskDetails,
        input::AbstractString)
    return "$(task.basic_info.base_directory)/$(input)"
end

function full_output_path(task::CubeStackTaskDetails,
        output::AbstractString)
    return "$(task.basic_info.base_directory)/$(output)";
end

function DaemonTask.prepare(task::CubeStackTaskDetails,
        datasource::DatasourceService)
    Datasource.get(datasource,
        map((input) -> full_input_path(task, input), task.basic_info.inputs); override_cache = true)
    Main.reload_registries();
end

function DaemonTask.execute(task::CubeStackTaskDetails,
        datasource::DatasourceService)
    inputs = task.basic_info.inputs

    if length(inputs) == 0
        return DaemonTask.Result(true, [])
    end

    index = task.payload_info.indices[1];
    indices = Main.get_indices(Main.aligned(0,0))
    z = findfirst(indices, index)
    cube_origin = (-16818, -5872, z)
    cube_extent = (87608, 56257, z)
    cube_dims = (512,512,1)
    overlap = (0,0,0)
    grid_dims = [0,0,0]
    for (k, (o, e, d)) = enumerate(zip(cube_origin, cube_extent, cube_dims))
        grid_dims[k] = ceil(Int64, (e-o+.indices1)/d)
    end
    grid_dims = (grid_dims...)
    grid_dims = (1,1,1)

    Main.save_cubes(cube_origin, cube_dims, overlap, grid_dims);
    src = joinpath(Main.FINISHED_DIR_PATH, "*.h5")
    dst = joinpath("gs://", Main.TASKS_BUCKET_NAME, Main.DATASET, Main.FINISHED_DIR)
    cmd = `gsutil -m cp $src $dst`
    Base.run(cmd)

    return DaemonTask.Result(true, task.payload_info.outputs)
end

function DaemonTask.finalize(task::CubeStackTaskDetails,
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

end # module CubeStackTask
