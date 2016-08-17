global TASKS_TASK_QUEUE_NAME = "task-queue-TEST";
global TASKS_ERROR_QUEUE_NAME = "error-queue-TEST";
global TASKS_REGISTRY_QUEUE_NAME = "registry-queue-TEST";
global TASKS_BUCKET_NAME = "seunglab";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;

type AlembicPayloadInfo
  	indices::Array{Index, 1} # array of input indices
	outputs::Array{Any, 1} # array of outputs
#	registry_updates::Array{Any, 1} # array of updates to the registry
end

function AlembicPayloadInfo{String <: AbstractString}(dict::Dict{String, Any})
  # do your parsing from JSON parsed dictionary here
  return AlembicPayloadInfo([tuple(index_array...) for index_array in dict["indices"]], dict["outputs"])
end

function push_registry_updates(queue_name = TASKS_REGISTRY_QUEUE_NAME)
    env = AWS.AWSEnv()
    queue = SimpleTasks.Services.AWSQueue.AWSQueueService(env, queue_name)

    # create tasks from the inputs and add them to the queue
    while length(REGISTRY_UPDATES) != 0
    SimpleTasks.Services.AWSQueue.Queue.push_message(queue; message_body = JSON.json(shift!(REGISTRY_UPDATES)));
    end
end

function pull_registry_updates(queue_name = TASKS_REGISTRY_QUEUE_NAME)
    env = AWS.AWSEnv()
    queue = SimpleTasks.Services.AWSQueue.AWSQueueService(env, queue_name)

    while true
      update_message = SimpleTasks.Services.AWSQueue.Queue.pop_message(queue)
      if length(update_message) == 0 break end
    # create tasks from the inputs and add them to the queue
      dict = JSON.parse(update_message);
      update_registry(tuple(dict["index"]...); rotation = Float64(dict["rotation"]), offset = Array{Int64, 1}(dict["offset"]), image_size = Array{Int64, 1}(dict["image_size"]), rendered = Bool(dict["rendered"])); 
    end
end
