global BUCKET = joinpath(homedir(), "seunglab")
global DATASET = "datasets/piriform"
global ROI_FIRST = (1,2,0,0);
global ROI_LAST = (8,173,0,0);
global DATASET_RESOLUTION = [7,7,40]

global TASKS_LOCALE = "gcs"
#global TASKS_LOCALE = "aws"

if TASKS_LOCALE == "aws"
global TASKS_TASK_QUEUE_NAME = "task-queue-TEST";
global TASKS_ERROR_QUEUE_NAME = "error-queue-TEST";
global TASKS_DONE_QUEUE_NAME = "done-queue-TEST";
global TASKS_REGISTRY_QUEUE_NAME = "registry-queue-TEST";
global TASKS_BUCKET_NAME = "seunglab";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;
end

if TASKS_LOCALE == "gcs"
global TASKS_TASK_QUEUE_NAME = "task-queue-pinky";
global TASKS_ERROR_QUEUE_NAME = "error-queue-pinky";
global TASKS_DONE_QUEUE_NAME = "done-queue-pinky";
global TASKS_REGISTRY_QUEUE_NAME = "registry-queue-pinky";
#=global TASKS_TASK_QUEUE_NAME = "task-queue-GCS";
global TASKS_ERROR_QUEUE_NAME = "error-queue-GCS";
global TASKS_DONE_QUEUE_NAME = "done-queue-GCS";
global TASKS_REGISTRY_QUEUE_NAME = "registry-queue-GCS"; =#
global TASKS_BUCKET_NAME = "image_assembly";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;
end