#=Pkg.clone("git@github.com:JuliaParallel/AWS.jl.git")=#
#=Pkg.clone("git@github.com:JuliaLang/JSON.jl.git")=#
module RunDaemon

import Daemon
import QueueService

type RunConfig
    queue_name::ASCIIString
    bucket_name::ASCIIString
    poll_frequency_seconds::Int
end

#=
 = Create the queue and bucket service and start the daemon
 =#
function main()
    run_config = parse_args()

    # Load AWS credentials via AWS library (either through environment
    # variables or ~/.awssecret or query permissions server)
    env = AWS.AWSEnv();

    queue = QueueService.AWSQueueService(env, run_config.queue_name)

    bucket = BucketService.AWSBucketService(env, run_config.bucket_name)

    daemon = Daemon.Daemon(queue, bucket, run_config.poll_frequency_seconds)

    daemon.run()
end

#=
 =Parse ARG into run daemon configuration.
 =Return RunConfig
 =#
function parse_args()
    if length(ARGS) < 3
        error("Not enough arguments given, sample usage:
            -- julia daemon.jl queue_name bucket_name")
    end

    return RunConfig(
        queue_name=ASCIIString(ARGS[1]),
        bucket_name = ASCIIString(ARGS[2])
    )
end

#Run main!
function __init__()
    main()
end

end # end module RunDaemon
