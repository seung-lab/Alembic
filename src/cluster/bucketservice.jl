module BucketService

import AWS
import AWS.S3

export BucketService
export AWSQueueService
export download
export upload

abstract BucketService
type AWSBucketService <: BucketService
    prefix::ASCIIString
    name::ASCIIString
    env::AWS.AWSEnv

    function AWSBucketService(env::AWSBucketService, name::ASCIIString,
        base_directory::ASCIIString)
        this = new()
        this.env = env
        this.name = name
        this.prefix = "s3://"

        # Check to make sure bucket is reachable
        check_reachable(env, name)
    end
end

function check_reachable(env::AWS.AWSEnv, bucket_name::AbstractString)
    bucket_response = S3.get_bkt(env, bucket_name)

    if bucket_response.http_code != 200
        error("Unable to access bucket: $bucket_name, response:
            ($(bucket_response.http_code))")
    end
end

function download(bucket::AWSBucketService, remote_file::ASCIIString,
    local_file::Union{ASCIIString, IO})
    get_response = S3.get_object(bucket.env, bucket.name, remote_file)

    if get_response.http_code != 200
        error("Unable to access file: $remote_file, response:
            ($(bucket_response.http_code))")
    end

    if isa(local_file, ASCIIString)
        local_file = open(local_file, "w")
    end

    write(local_file, get_response.obj)
end

function upload(bucket::AWSBucketService, local_file::Union{ASCIIString, IO},
    remote_file::ASCIIString)
    folder = remote_file[1:rsearch(remote_file, "/").stop]

    #check to make sure folder exists TODO does it matter?!
    if !isempty(folder)
        folder_response = S3.put_object(bucket.env, bucket.name, folder,
            "folder")
        if folder_response != 200
            error("Unable to access folder $folder in bucket $(bucket.name)")
        end
    end

    # TODO FIX THIS
    put_response = S3.put_object(bucket.env, bucket.name, local_file)

    if put_response != 200
        error("Unable to put object $local_file in bucket $(bucket.name),
            response was $(put_response.http_code)")
    end

end
end # end module BucketService
