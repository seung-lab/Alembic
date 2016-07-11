module CloudTest
include("utils/mockservices.jl")
include("utils/mocktasks.jl")
include("tasks/test_daemontask.jl")
include("tasks/test_blockmatchtask.jl")
include("tasks/test_basictaskinfo.jl")
include("tasks/test_alignmenttaskinfo.jl")
include("test_daemonservice.jl")
end # module CloudTest
