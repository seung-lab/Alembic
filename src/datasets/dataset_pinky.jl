global BUCKET = "/home/ubuntu"
global DATASET = "datasets/pinky_40percent"
#global DATASET = "datasets/pinky_cropped"
#global DATASET = "datasets/pinky_demo"
global ROI_FIRST = (1,3484,0,0);
global ROI_LAST = (1,4491,0,0);
global DATASET_RESOLUTION = [4,4,40]

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
global TASKS_BUCKET_NAME = "seunglab_alembic";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;
end
