#include("../SimpleTasks.jl")

module AWSScheduler

using SimpleTasks.Types
using SimpleTasks.Services.AWSQueue
using SimpleTasks.Services.CLIBucket
using SimpleTasks.Services.AWSCLIProvider
using BlockMatchTask

import AWS
import JSON
import SimpleTasks.Services.Bucket
import SimpleTasks.Services.Queue
import SimpleTasks.Tasks.BasicTask
import Main.AlembicPayloadInfo

# creates a task to make a MeshSet for index
function create_task(index::Main.Index)
	if Main.is_montaged(index)	indices = Main.get_index_range(Main.prevstage(index),Main.prevstage(index))
	elseif Main.is_prealigned(index) indices = [Main.prevstage(index), Main.get_preceding(Main.prevstage(index))]
	end

	inputs_images = map(Main.truncate_path, map(Main.get_path, indices));
	inputs_registry = map(Main.truncate_path, map(Main.get_registry_path, indices));
	inputs = unique(vcat(inputs_images, inputs_registry))
	
	output_meshset = Main.truncate_path(Main.get_path("MeshSet", index))
	output_stats = Main.truncate_path(Main.get_path("stats", index))

	basic_info = BasicTask.Info(0, NAME, Main.TASKS_BASE_DIRECTORY, inputs) 
	task = BlockMatchTaskDetails(basic_info, AlembicPayloadInfo([index], [output, output_stats]));
#	return vcat(inputs..., output)
	return task
end

function schedule_blockmatch(index; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = create_task(index);
    Queue.push_message(queue; message_body = JSON.json(task));
#    tasks = map(create_task, indices[1:end-1])
#    map((task) -> Queue.push_message(queue; message_body = JSON.json(task)),
#        tasks)
end
#=
function schedule(queue_name, bucket_name)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create data and upload it to the bucket service
    indices = 0:9
    inouts = create_input_files(indices)
    map((inout) -> Bucket.upload(bucket, inout[1], inout[2]), inouts)

    # create tasks from the inputs and add them to the queue
    tasks = map(create_task, indices[1:end-1])
    map((task) -> Queue.push_message(queue; message_body = JSON.json(task)),
        tasks)
end

function __init__()
    schedule("task-queue-TEST", "seunglab-test")
end
=#

end # module AWSScheduler
