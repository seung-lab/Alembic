module BasicTask

export Info

"""
    BasicTask.Info

Type contains basic information for a daemon task

"""
type Info
    id::Int64
    name::AbstractString
    baseDirectory::AbstractString
    inputs::Array{AbstractString, 1}
end

function Info(dict::Dict{AbstractString, Any})
    if haskey(dict, "id")
        id = typeof(dict["id"]) <: Int ?  dict["id"] : parse(Int64, dict["id"])
    else
        id = -1
    end

    if isempty(strip(dict["name"]))
        throw(ArgumentError("Task name can not be empty"))
    end

    # parse base directory
    if isempty(strip(dict["baseDirectory"]))
        throw(ArgumentError("Payload does not include a baseDirectory"))
    end

    # parse input list
    if typeof(dict["inputs"]) != Array{Any, 1} ||
            length(dict["inputs"]) == 0
        throw(ArgumentError("Payload does not include a input list"))
    end

    return Info(id, dict["name"], dict["baseDirectory"], dict["inputs"])
end

end # module BasicTask
