#include("../SimpleTasks.jl")

module AWSScheduler

using SimpleTasks.Types
using SimpleTasks.Services.AWSQueue
using SimpleTasks.Services.CLIBucket
#using SimpleTasks.Services.AWSCLIProvider
using ImportTask
using BlockMatchTask
using RenderTask
using SolveTask
using ThumbnailTask
using RenderReviewTask
using SaveStackTask

import AWS
import JSON
import SimpleTasks.Services.Bucket
import SimpleTasks.Services.Queue
import SimpleTasks.Tasks.BasicTask
import Main.AlembicPayloadInfo

function schedule_import(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = ImportTask.ImportTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
end


function schedule_blockmatch(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = BlockMatchTask.BlockMatchTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
end

function schedule_render(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = RenderTask.RenderTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
end

function schedule_solve(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = SolveTask.SolveTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
end

function schedule_thumbnail(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = ThumbnailTask.ThumbnailTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
end

function schedule_render_review(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = RenderReviewTask.RenderReviewTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
end

function schedule_save_stack(args...; queue_name = Main.TASKS_TASK_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
#    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    # create tasks from the inputs and add them to the queue
    task = RenderReviewTask.SaveStackTaskDetails(args...);
    Queue.push_message(queue; message_body = JSON.json(task));
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

function collect_errors(queue_name = Main.TASKS_ERROR_QUEUE_NAME, bucket_name = Main.TASKS_BUCKET_NAME)
    env = AWS.AWSEnv()
    queue = AWSQueueService(env, queue_name)
    m = Queue.pop_message(queue)
    while m != ""
        d = JSON.parse(m)
        indices = [(ind...) for ind in d["task"]["payload_info"]["indices"]]
        index = indices[1]
        fn = joinpath(PREMONTAGED_DIR_PATH, "error", string(index, ".txt"))
        writedlm(fn, d)
        m = Queue.pop_message(queue)
    end
end

end # module AWSScheduler
  