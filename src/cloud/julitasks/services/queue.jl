module Queue

using ...Julitasks.Types

export pop_message

"""
    pop_message(queue::QueueService)

Pop a message of the queue.
"""
function pop_message(queue::QueueService)
    error("pop_message for $queue is not implemented")
end

end # module Queue
