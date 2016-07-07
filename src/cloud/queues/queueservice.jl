module Queue

export QueueService, pop_message

abstract QueueService

"""
    pop_message(queue::QueueService)

Pop a message of the queue. Override this function in implementation
"""
function pop_message(queue::QueueService)
    error("pop_message for $queue is not implemented")
end

end # module Queue
