module Services
include("bucket.jl")
include("awsbucket.jl")
include("awsclibucket.jl")

include("queue.jl")
include("awsqueue.jl")

include("cache.jl")
include("filesystemcache.jl")

include("datasource.jl")
include("bucketcachedatasource.jl")

include("daemon.jl")
end # module Services
