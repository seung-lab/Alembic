module Queue

import AWS
import AWS.SQS

export QueueService
export AWSQueueService
export pop_message

abstract QueueService

type AWSQueueService <: QueueService
    name::ASCIIString
    url::ASCIIString
    env::AWS.AWSEnv

    AWSQueueService(env::AWS.AWSEnv, name::ASCIIString) =
        new(name, get_queue_url(env, name), env)
end

# Given our aws environment and bucket name, find the correct url
function get_queue_url(env::AWS.AWSEnv, queue_name::AbstractString)
    response = SQS.CreateQueue(env, queueName=queue_name)

    # try creating a queue, does not overwrite existing queue so it's safe
    if response.http_code != 200
        error("Unable to create/get queue: $queue_name, response:
            ($(response.http_code))")
    end

    return response.obj.queueUrl
end

# Pops message off the given queue service and returns a task
function pop_message(queue::AWSQueueService)
    receive_response = SQS.ReceiveMessage(queue.env, queue.url)

    if receive_response.http_code != 200
        error("Unable to retrieve task from $(queue.name) from $(queue.url)")
    end

    if length(receive_response.obj.messageSet) < 1
        error("Invalid message received from $(queue.name) from $(queue.url)")
    end

    delete_response = SQS.DeleteMessage(env, queue.url,
    receiptHandle=receive_response.obj.messageSet[1].receiptHandle)

    return strip(receive_response.obj.messageSet[1].body)
end

end # end module QueueService
