import JSON
import AWS
import task

type QueueService
    name::ASCIIString
    url::ASCIIString
end

function receive_message(queue::GCSQueueService)
    error("NOT IMPLEMENTED")
end

# Pops message off the given queue service and returns a task
function pop_task(queue::AWSQueueService)
    receive_response = receive_message(queue)

    if receive_response.http_code != 200
        error("Unable to retrieve task from $(queue.name) from $(queue.url)")
    end

    if length(receive_response.obj.messageSet) < 1
        || length(strip(receive_response.obj.messageSet.body))
        error("Invalid message received from $(queue.name) from $(queue.url)")
    end

    delete_response = AWS.SQS.DeleteMessage(env, queue.url,
    receiptHandle=receive_response.obj.messageSet[1].receiptHandle)
    task = parse_task(strip(receive_response.obj.messageSet[1].body))
    return task
end

function parse_task(message::ASCIIString)
    if length(strip(sqs_response.obj.messageSet.body))
        || !haskey(message, "name")  || !haskey(message, "indices")
        error("Invalid message body retrieved")
    end


    dictionary = JSON.parse(recieve_response.obj.messageSet.body)
    return task(dictionary["name"], dictionary["indicies"])
end


function get_task(queue::QueueService)
    return task;
