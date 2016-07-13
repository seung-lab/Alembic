module Queue

export Service, pop_message

abstract Service

"""
    pop_message(queue::QueueService)

Pop a message of the queue.
"""
function pop_message(queue::Service)
    error("pop_message for $queue is not implemented")
end

end # module Queue
