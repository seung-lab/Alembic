import HttpCommon
import AWS

type AWSQueueService <: QueueService
    env::AWS.AWSEnv
end

function pop_task(queue::AWSQueueService)
    sqs_response = AWS.SQS.ReceiveMessage(queue.env, queue.url)
    if length(sqs_response.obj.messageSet) < 1
        error("Invalid message received from $(queue.name) from $(queue.url)")
    end
    return HttpResponse.Response(status = sqs_response.http_code, kkj
end

