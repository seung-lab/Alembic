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

function upload_remote_test_file()
    create_local_test_file(TEST_FILE_NAME)
    try
        run(`aws s3 cp $TEST_FILE_NAME s3://$BUCKET_NAME/$TEST_FILE_NAME`)
    finally
        delete_local_test_file()
    end
end

function delete_remote_test_file()
    run(`aws s3 rm s3://$BUCKET_NAME/$TEST_FILE_NAME`)
end

function write_numbers(io::IO)
    for i in 1:MAX_INDEX
        write(io, "$i\n")
    end
end

function create_local_test_file(filename::AbstractString)
    file = open(filename, "w")
    write_numbers(file)
    close(file)
end

function delete_local_test_file()
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
    upload_remote_test_file()

    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    buffer = IOBuffer()
    Bucket.download(bucket, TEST_FILE_NAME, buffer)

    file_indices = split(takebuf_string(buffer), "\n")
    for index in 1:length(file_indices) - 1
        @test parse(Int, file_indices[index]) == index
    end
    # last split from \n is generates an empty string
    @test length(file_indices) - 1 == MAX_INDEX

    delete_remote_test_file()
end

function test_download_file()
    upload_remote_test_file()

    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    download_filename = "new$TEST_FILE_NAME"
    process = Bucket.download(bucket, TEST_FILE_NAME, download_filename)

    file = open(download_filename, "r")
    index = 0
    try
        while !eof(file)
            index = index + 1
            @test parse(Int, readline(file)) == index
        end
    finally
        close(file)
        #=rm(download_filename)=#
    end
    @test index == MAX_INDEX

    delete_remote_test_file()
end

function test_upload_io()
    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    io = IOBuffer()
    write_numbers(io)
    upload_filename = "$(TEST_FILE_NAME)UpIO"
    
    # make sure the file if it exists is removed first
    run(`aws s3 rm s3://$BUCKET_NAME/$upload_filename`)

    Bucket.upload(bucket, io, upload_filename)

    # Verify the uploaded file by manually download the uploaded file
    try
        run(`aws s3 cp s3://$BUCKET_NAME/$upload_filename $upload_filename`)
        run(`aws s3 rm s3://$BUCKET_NAME/$upload_filename`)
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
        rm(upload_filename)
    end
    @test index == MAX_INDEX
end

function test_upload_file()
    env = AWS.AWSEnv()
    bucket = AWSCLIBucketService(env.aws_id, env.aws_seckey,
        BUCKET_NAME)

    upload_filename = "$(TEST_FILE_NAME)UpFile"
    # make sure the file if it exists is removed first
    run(`aws s3 rm s3://$BUCKET_NAME/$upload_filename`)

    create_local_test_file(upload_filename)

    Bucket.upload(bucket, upload_filename, upload_filename)

    # Verify the uploaded file by manually download the uploaded file
    try
        run(`aws s3 cp s3://$BUCKET_NAME/$upload_filename $upload_filename`)
        run(`aws s3 rm s3://$BUCKET_NAME/$upload_filename`)
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
        rm(upload_filename)
    end
    @test index == MAX_INDEX
end

function __init__()
    test_creatable()
    test_bad_bucket()
    test_download_IO()
    test_download_file()

    test_upload_io()
    test_upload_file()
end

end # module TestAWSCLIBucket
