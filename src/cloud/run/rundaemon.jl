include("../Julimaps.jl")

module RunDaemon

using Julimaps.Cloud.Services.AWSQueue
using Julimaps.Cloud.Services.AWSBucket
using Julimaps.Cloud.Services.FileSystemCache
using Julimaps.Cloud.Services.BucketCacheDatasource
using Julimaps.Cloud.Services.Daemon

import Julimaps
import AWS

type RunConfig
    queue_name::ASCIIString
    bucket_name::ASCIIString
    poll_frequency_seconds::Int64
end

#=
 = Create the queue and bucket service and start the daemon
 =#
function main()
    run_config = parse_args()

    # Load AWS credentials via AWS library (either through environment
    # variables or ~/.awssecret or query permissions server)
    env = AWS.AWSEnv()

    queue =AWSQueueService(env, run_config.queue_name)

    bucket = AWSBucketService(env, run_config.bucket_name)

    cache = FileSystemCacheService("

    daemon = Daemon.DaemonService(queue, bucket, run_config.poll_frequency_seconds)

    Daemon.run(daemon)
end

#=
 =Parse ARG into run daemon configuration.
 =Return RunConfig
 =#
function parse_args()
    if length(ARGS) < 3
        error("Not enough arguments given, (given $ARGS) sample usage:
            -- julia daemon.jl queue_name bucket_name")
    end

    run_config = RunConfig(
        ASCIIString(ARGS[1]),
        ASCIIString(ARGS[2]),
        parse(Int64, ARGS[3])
    )
    return run_config
end

#Run main!
function __init__()
    main()
end

end # end module RunDaemon
