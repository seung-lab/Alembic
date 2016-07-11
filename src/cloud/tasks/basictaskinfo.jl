module BasicTask

export Info

"""
    DaemonTask.Info

Type contains basic information for a daemon task

"""
type Info
    id::Int64
    name::AbstractString
    baseDirectory::AbstractString
    files::Array{AbstractString, 1}
end

function Info(dict::Dict{AbstractString, Any}) 
    id = typeof(dict["id"]) <: Int ?  dict["id"] : parse(Int64, dict["id"])
    if isempty(strip(dict["name"]))
        throw(ArgumentError("Task name can not be empty"))
    end

    # parse base directory
    if isempty(strip(dict["baseDirectory"]))
        throw(ArgumentError("Payload does not include a baseDirectory"))
    end

    # parse file list
    if typeof(dict["files"]) != Array{Any, 1} ||
            length(dict["files"]) == 0
        throw(ArgumentError("Payload does not include a file list"))
    end

    return Info(id, dict["name"], dict["baseDirectory"], dict["files"])
end

end # module BasicTask
