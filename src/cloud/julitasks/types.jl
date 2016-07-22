module Types

export BucketService,
    DaemonService,
    DatasourceService,
    QueueService,
    CacheService

abstract BucketService
abstract DatasourceService
abstract QueueService
abstract CacheService

export DaemonTaskDetails

"""
    DaemonTaskDetails

This is the base composite abstract class used to compose Details and Payload
i.e. compose a task with
```julia
type YourDaemonTaskDetails <: DaemonTaskDetails
    basicInfo::BasicTask.Info
    taskInfo::YourTask.Info
end
```
"""
abstract DaemonTaskDetails

end # module Types
