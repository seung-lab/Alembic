module AWSCLIBucket

using ...Julitasks.Types

import AWS, AWS.S3
import Julimaps.Cloud.Julitasks.Services.Bucket

export AWSCLIBucketService

const AWS_KEY_FILE = "$(homedir())/.aws/config"

type AWSCLIBucketService <: BucketService
    name::ASCIIString
end

AWSCLIBucketService(aws_access_key_id::ASCIIString,
    aws_secret_access_key::ASCIIString, name::ASCIIString) =
    create_access!(aws_access_key_id, aws_secret_access_key) && 
    check_reachable(name) && AWSCLIBucketService(name)

function create_access!(aws_access_key_id::ASCIIString,
        aws_secret_access_key::ASCIIString)
    if stat(AWS_KEY_FILE).size == 0
        println("No access file found for aws cli, creating one!")
        file = open(AWS_KEY_FILE, "w")
        write(file, "[default]\n")
        write(file, "aws_access_key_id = $(aws_access_key_id)\n")
        write(file, "aws_secret_access_key = $(aws_secret_access_key)\n")
        close(file)
    end
    return true
end

function check_reachable(bucket_name::AbstractString)
    try
        println("Checking access to $bucket_name")
        run(pipeline(`aws s3 ls s3://$bucket_name`, DevNull))
    catch
        throw(ArgumentError("Unable to access bucket \"$bucket_name\"")) 
    end
    return true
end

function Bucket.download(bucket::AWSCLIBucketService, remote_file::ASCIIString,
    local_file::Union{ASCIIString, IO})

    if isa(local_file, ASCIIString)
        local_file = open(local_file, "w")
    end

    download_cmd = `aws s3 cp s3://$(bucket.name)/$remote_file -`

    process = nothing
    # ugh julia doesn't support piping directly to IOBuffers yet
    #=https://github.com/JuliaLang/julia/issues/14437=#
    if typeof(local_file) <: IOBuffer
        (s3_output, process) = open(download_cmd, "r")
        # manually read from stream and write to buffer
        write(local_file, readbytes(s3_output))
    else
        # open the cmd in write mode. this automatically takes the 2nd arg
        # (stdio) and uses it as redirection of STDOUT
        (s3_input, process) = open(download_cmd, "w", local_file)
        # for now just make it block until command has completed until i know
        # how to make IOBuffer async/return a process
        wait(process)
    end

    return process
end

function Bucket.upload(bucket::AWSCLIBucketService,
        local_file::Union{ASCIIString, IO}, remote_file::ASCIIString)

    if isa(local_file, ASCIIString)
        local_file = open(local_file, "w")
    end

    upload_cmd = `aws s3 cp - s3://$(bucket.name)/$remote_file`

    process = nothing
    # ugh julia doesn't support piping directly to IOBuffers yet
    #=https://github.com/JuliaLang/julia/issues/14437=#
    if typeof(local_file) <: IOBuffer
        (s3_input, process) = open(upload_cmd, "w")
        # manually read from buffer and write to stream
        write(s3_input, takebuf_array(local_file))
        close(s3_input)
    else
        # open the cmd in write mode. this automatically takes the 2nd arg
        # (stdio) and uses it as redirection of STDOUT.
        (s3_output, process) = open(upload_cmd, "r", local_file)
        # for now just make it block until command has completed until i know
        # how to make IOBuffer async/return a process
        wait(process)
    end

    return process
end

#=
 =    folder = remote_file[1:rsearch(remote_file, "/").stop]
 =
 =    #check to make sure folder exists TODO does it matter?!
 =    if !isempty(folder)
 =        folder_response = S3.put_object(bucket.env, bucket.name, folder,
 =            "folder")
 =        if folder_response != 200
 =            error("Unable to access folder $folder in bucket $(bucket.name)")
 =        end
 =    end
 =
 =    put_response = S3.put_object(bucket.env, bucket.name, local_file)
 =
 =    if put_response != 200
 =        error("Unable to put object $local_file in bucket $(bucket.name),
 =            response was $(put_response.http_code)")
 =    end
 =#
#=end=#

end # module AWSCLIBucket
