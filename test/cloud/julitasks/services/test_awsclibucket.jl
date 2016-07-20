module TestAWSCLIBucket

using Base.Test
using Julimaps.Cloud.Julitasks.Services.AWSCLIBucket
#=using CloudTest.JulitasksTests.Utils.TestTasks=#
#=using CloudTest.JulitasksTests.Utils.MockServices=#

import AWS
import Julimaps.Cloud.Julitasks.Services.Bucket

const TEST_FILE_NAME = "testfile"
const BUCKET_NAME = "seunglab-alignment"
const MAX_INDEX = 100000

function upload_test_file()
    create_test_file()
    try
        run(`aws s3 cp $TEST_FILE_NAME s3://$BUCKET_NAME/$TEST_FILE_NAME`)
    finally
        delete_test_file()
    end
end

function write_numbers(io::IO)
    for i in 1:MAX_INDEX
        write(io, "$i\n")
    end
end

function create_test_file()
    file = open(TEST_FILE_NAME, "w")
    write_numbers(file)
    close(file)
end

function delete_test_file()
    rm(TEST_FILE_NAME)
end

function test_creatable()
    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)
    @test bucket != nothing
end

function test_bad_bucket()
    env = AWS.AWSEnv()
    @test_throws ArgumentError AWSCLIBucketService(env.aws_id, env.aws_seckey,
        "$BUCKET_NAME-bad")
end

function test_download_IO()
    upload_test_file()

    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    buffer = IOBuffer()
    Bucket.download(bucket, TEST_FILE_NAME, buffer)

    index = 0
    while !eof(buffer)
        index = index + 1
        @test parse(Int, readline(buffer)) == index
    end
end

function test_download_file()
    upload_test_file()

    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    new_filename = "new$TEST_FILE_NAME"
    Bucket.download(bucket, TEST_FILE_NAME, new_filename)

    file = open(new_filename, "r")
    index = 0
    try
        while !eof(file)
            index = index + 1
            @test parse(Int, readline(file)) == index
        end
    finally
        close(file)
        rm(new_filename)
    end
    @test index == MAX_INDEX
end

function test_upload_io()
    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    io = IOBuffer()
    write_numbers(io)
    upload_filename = "$(TEST_FILE_NAME)up"

    Bucket.upload(bucket, io, "$upload_filename")

    try
        run(`aws s3 cp s3://$BUCKET_NAME/$upload_filename $upload_filename`)
    catch
        error("Unable to find downloaded file $upload_filename")
    end

    downloaded_file = open(upload_filename, "r")
    index = 0
    try
        while !eof(downloaded_file)
            index = index + 1
            @test parse(Int, readline(downloaded_file)) == index
        end
    finally
        close(downloaded_file)
    end
    @test index == MAX_INDEX
end

function test_upload_file()
end

function __init__()
    test_creatable()
    test_bad_bucket()
    test_download_IO()
    test_download_file()

    test_upload_io()
end

end # module TestAWSCLIBucket
