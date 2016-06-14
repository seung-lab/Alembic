module Task

import JSON
export Task
export parse_task

type Task
    name::AbstractString
    indices::AbstractString
    base_directory::AbstractString
end

function parse_task(message::ASCIIString)
    if length(strip(message))
        error("Trying to parse empty string for task")
    end

    if !haskey(message, "name")  || !haskey(message, "indices")
        error("Missing task parameters")
    end

    dictionary = JSON.parse(recieve_response.obj.messageSet.body)
    return task(task=dictionary["name"], indices=dictionary["indicies"])
end

end # end module Task
