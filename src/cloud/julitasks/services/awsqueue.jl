module AWSQueue

using ...Julitasks.Types

import AWS, AWS.SQS
import Julimaps.Cloud.Julitasks.Services.Queue

export AWSQueueService

type AWSQueueService <: QueueService
    env::AWS.AWSEnv
    name::ASCIIString
    url::ASCIIString
end

AWSQueueService(env::AWS.AWSEnv, name::ASCIIString) =
        AWSQueueService(env, name, get_queue_url(env, name))

function Base.string(queue::AWSQueueService)
    return "$(queue.name) from $(queue.url)"
end

"""
    get_queue_url(env::AWS.AWSEnv, queue_name::AbstractString)

Find the correct url with our aws environment and bucket name
"""
function get_queue_url(env::AWS.AWSEnv, queue_name::AbstractString)
    response = SQS.CreateQueue(env, queueName=queue_name)

    # try creating a queue, does not overwrite existing queue so it's safe
    if response.http_code != 200
        error("Unable to create/get queue: $queue_name, response:
            ($(response.http_code))")
    end

    return response.obj.queueUrl
end

"""
    Queue.pop_message(queue::AWSQueueService)

Pop a message of the aws queue. Supplied to Queue for multi dispatch
"""
function Queue.pop_message(queue::AWSQueueService)
    receive_response = SQS.ReceiveMessage(queue.env; queueUrl = queue.url)

    if receive_response.http_code != 200
        error("Unable to retrieve task from $queue")
    end

    if length(receive_response.obj.messageSet) < 1
        return ""
    end

    receiptHandle=receive_response.obj.messageSet[1].receiptHandle

    delete_response = SQS.DeleteMessage(queue.env; queueUrl = queue.url,
        receiptHandle = receiptHandle)

    return strip(receive_response.obj.messageSet[1].body)
end

end # module AWSQueue
