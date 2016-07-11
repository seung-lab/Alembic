module BasicTask

export Info

"""
    DaemonTask.Info

Type contains basic information for a daemon task

"""
type Info
    id::Int64
    name::AbstractString
end

function Info(dict::Dict{AbstractString, Any}) 
    id = typeof(dict["id"]) <: Int ?  dict["id"] : parse(Int64, dict["id"])
    if isempty(strip(dict["name"]))
        throw(ArgumentError("Task name can not be empty"))
    end
    return Info(id, dict["name"])
end

end # module BasicTask
