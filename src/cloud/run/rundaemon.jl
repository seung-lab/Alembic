include("../../Julimaps.jl")

module RunDaemon

using Julimaps.Cloud.Julitasks.Types
using Julimaps.Cloud.Julitasks.Services.AWSQueue
using Julimaps.Cloud.Julitasks.Services.AWSCLIBucket
using Julimaps.Cloud.Julitasks.Services.FileSystemCache
using Julimaps.Cloud.Julitasks.Services.BucketCacheDatasource
using Julimaps.Cloud.Julitasks.Services.Daemon
using Julimaps.Cloud.Julitasks.Examples.NoOpTask

import Julimaps
import AWS

type RunConfig
    queue_name::ASCIIString
    bucket_name::ASCIIString
    cache_directory::ASCIIString
    poll_frequency_seconds::Int64
end

#=
 = Create the queue and bucket service and start the daemon
 =#
function run(queue_name, bucket_name, cache_directory, poll_frequency_seconds)
    # Load AWS credentials via AWS library (either through environment
    # variables or ~/.awssecret or query permissions server)
    env = AWS.AWSEnv()

    queue = AWSQueueService(env, queue_name)

    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey, bucket_name)

    cache = FileSystemCacheService(cache_directory)

    datasource = BucketCacheDatasourceService(bucket, cache)

    daemon = DaemonService(queue, bucket, datasource, poll_frequency_seconds)

    register!(daemon, NoOpTask.NAME, NoOpTaskDetails)

    Daemon.run(daemon)
end

#=
 =Parse ARG into run daemon configuration.
 =Return RunConfig
 =#
function parse_args()
    if length(ARGS) < 3
        error("Not enough arguments given, (given $ARGS) sample usage:
            -- julia daemon.jl queue_name bucket_name cache_directory
            poll_frequency_seconds")
    end

    run_config = RunConfig(
        ASCIIString(ARGS[1]),
        ASCIIString(ARGS[2]),
        ASCIIString(ARGS[3]),
        parse(Int64, ARGS[4])
    )
    return run_config
end

#Run main!
function __init__()
    #=run_config = parse_args()=#
    run_config = RunConfig("test-queue3", "seunglab-alignment", "/var/tmp/", 20)
    run(run_config.queue_name, run_config.bucket_name,
        run_config.cache_directory, run_config.poll_frequency_seconds)
end

end # end module RunDaemon
