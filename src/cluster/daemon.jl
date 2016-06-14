module Daemon
import AWS
using QueueService


function run(queue::QueueService, bucket::BucketService,
    poll_frequency_seconds::Int)

    while true
        message = pop_message(queue)
        task = parse_task(message)
        print("Task is $(task.name) with indicies $(task.indices)")

        sleep(poll_frequency_seconds)
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

end # end module Daemon
