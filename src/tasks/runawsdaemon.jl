#include("../SimpleTasks.jl")

using Alembic

module RunAWSDaemon

using SimpleTasks.Types
using SimpleTasks.Services.AWSQueue
using SimpleTasks.Services.CLIBucket
using SimpleTasks.Services.AWSCLIProvider
using SimpleTasks.Services.FileSystemCache
using SimpleTasks.Services.BucketCacheDatasource
using SimpleTasks.Services.Daemon
using BlockMatchTask
using RenderTask
using SolveTask
using ImportTask

#using BlockMatchTasks

import AWS

type RunConfig
    task_queue_name::ASCIIString
    error_queue_name::ASCIIString
    done_queue_name::ASCIIString
    bucket_name::ASCIIString
    cache_directory::ASCIIString
    poll_frequency_seconds::Int64
end

#=
 = Create the queue and bucket service and start the daemon
 =#
function run(task_queue_name, error_queue_name, done_queue_name, bucket_name,
        cache_directory, poll_frequency_seconds)
    # Load AWS credentials via AWS library (either through environment
    # variables or ~/.awssecret or query permissions server)
    env = AWS.AWSEnv()

    task_queue = AWSQueueService(env, task_queue_name)

    error_queue = AWSQueueService(env, error_queue_name)

    done_queue = AWSQueueService(env, done_queue_name)

    bucket = CLIBucketService(AWSCLIProvider.Details(env), bucket_name)

    cache = FileSystemCacheService(cache_directory)

    datasource = BucketCacheDatasourceService(bucket, cache)

    daemon = DaemonService(task_queue, error_queue, done_queue, bucket, datasource,
        poll_frequency_seconds)

    register!(daemon, ImportTask.NAME, ImportTaskDetails)
    register!(daemon, BlockMatchTask.NAME, BlockMatchTaskDetails)
    register!(daemon, RenderTask.NAME, RenderTaskDetails)
    register!(daemon, SolveTask.NAME, SolveTaskDetails)

    Daemon.run(daemon)
end

#=
 =Parse ARG into run daemon configuration.
 =Return RunConfig
 =#
function parse_args()
    if length(ARGS) < 3
        error("Not enough arguments given, (given $ARGS) sample usage:
            -- julia daemon.jl task_queue_name error_queue_name bucket_name " *
            " cache_directory poll_frequency_seconds")
    end

    run_config = RunConfig(
        ASCIIString(ARGS[1]),
        ASCIIString(ARGS[2]),
        ASCIIString(ARGS[3]),
        ASCIIString(ARGS[4]),
        parse(Int64, ARGS[5])
    )
    return run_config
end

#Run main!
function __init__()
#    run_config = parse_args()
    #=
     =run_config = RunConfig("task-queue-TEST", "error-queue-TEST", "seunglab-test",
     =    "/var/tmp/", 5)
     =#
    run_config = RunConfig(
#=    	"task-queue-TEST",
    	"error-queue-TEST",
	"seunglab",
	joinpath(homedir(), "cache"),
	10=#
	Main.TASKS_TASK_QUEUE_NAME,
	Main.TASKS_ERROR_QUEUE_NAME,
	Main.TASKS_DONE_QUEUE_NAME,
	Main.TASKS_BUCKET_NAME,
	Main.TASKS_CACHE_DIRECTORY,
	Main.TASKS_POLL_FREQUENCY
    )
    run(run_config.task_queue_name, run_config.error_queue_name, run_config.done_queue_name,
        run_config.bucket_name, run_config.cache_directory,
        run_config.poll_frequency_seconds)
end

end # end module RunDaemon
