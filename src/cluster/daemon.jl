import AWS

function main()
    # Load AWS credentials via AWS library (either through environment
    # variables or ~/.awssecret or query permissions server)
    env = AWS.AWSEnv();

    if length(ARGS) < 3
        error("Not enough arguments given, sample usage:
            -- julia daemon.jl queue_name bucket_name base_directory")
    end

    eventsAndFunctions
    queue_name = ASCIIString(ARGS[1])
    bucket_name = ASCIIString(ARGS[2])
    base_directory = ARGS[3]

    # detect queue
    queue_url = get_queue_url(env, queue_name)

    # detect base directory
    verify_directory(env, bucket_name, base_directory)

end

function get_queue_url(env::AWS.AWSEnv, queue_name::AbstractString)
    # try creating it anyway, does not overwrite existing queue
    response = AWS.SQS.CreateQueue(env, queueName=queue_name)
    if response.http_code != 200
        error("Unable to find queue: $queue_name, response:
            ($(response.http_code)), creating a new one")
    end
    return response.obj.queueUrl
end

function verify_directory(env::AWS.AWSEnv, bucket_name::AbstractString,
    base_directory::AbstractString)
    bucket_response = AWS.S3.get_bkt(env, bucket_name)

    if bucket_response.http_code != 200
        error("Unable to access bucket: $bucket_name, response:
            ($(bucket_response.http_code))")
    end

    base_directory_response = AWS.S3.get_object(env, bucket_name, "$base_directory/")
    if base_directory_response.http_code != 200
        error("Unable to locate base directory: $base_directory, response:
            ($(base_directory_response.http_code))")
    end
end

main()
