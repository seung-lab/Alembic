module AlignmentTask

export Info

"""
    AlignmentTask.Info

This module is meant as a composable type for Alignment tasks such
BlockMatchTask and RenderTask

"""
type Info
    indices::Array{Tuple{Int64,Int64,Int64,Int64}}
end

function Info(dict::Dict{AbstractString, Any})
    # parse indices list
    if typeof(dict["indices"]) != Array{Any, 1} ||
            length(dict["indices"]) == 0
        throw(ArgumentError("Payload does not include indices"))
    end
    indices = Array{Tuple{Int64, Int64, Int64, Int64}, 1}()
    for index = dict["indices"]
        if length(index) != 4 # tuple size
            throw(ArgumentError("Index needs to be length of 4, received
                $index"))
        end
        push!(indices, (index[1], index[2], index[3], index[4]))
    end

    return Info(indices)
end

end # module AlignmentTask
