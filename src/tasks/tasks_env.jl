global TASKS_TASK_QUEUE_NAME = "task-queue-TEST";
global TASKS_ERROR_QUEUE_NAME = "error-queue-TEST";
global TASKS_ERROR_QUEUE_NAME = "registry-queue-TEST";
global TASKS_BUCKET_NAME = "seunglab";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;

type AlembicPayloadInfo
  	indices::Array{Index, 1} # array of input indices
	outputs::Array{Any, 1} # array of outputs
#	registry_updates::Array{Any, 1} # array of updates to the registry
end
