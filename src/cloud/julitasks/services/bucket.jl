module Bucket

using ...Julitasks.Types

export download, upload

"""
    download(bucket::BucketService, remote_file::AbstractString,
local_file::Union{AbstractString, IO})

Download a remote file either to a new location `AbstractString` or a stream `IO`.
If the local_file is not specified, a stream to the download is returned
"""
function download(bucket::BucketService, remote_file::AbstractString,
    local_file::Union{AbstractString, IO, Void}=nothing)
    error("download with $bucket is not implemented")
end

"""
    upload(bucket::BucketService, local_file::Union{AbstractString, IO},
remote_file::AbstractString)

Upload a file `AbstractString` or a stream `IO` to the bucket service.
"""
function upload(bucket::BucketService, local_file::Union{AbstractString, IO},
    remote_file::AbstractString)
    error("upload with $bucket is not implemented")
end

end # module Bucket
