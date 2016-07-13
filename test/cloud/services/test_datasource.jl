module TestDataSource

using Base.Test
using CloudTest.MockServices

import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache
import Julimaps.Cloud.Services.DataSource

function test_pull_empty_cache()
    key = "somekey"

    bucket_text = "mock contents"
    bucket_file = IOBuffer(bucket_text)
    seekstart(bucket_file)
    bucket_files = Dict()
    bucket_files[key] = bucket_file

    cache_values = Dict()

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)

    datasource = DataSource.Service(bucket, cache)

    DataSource.pull!(datasource, key)

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    @test readchomp(new_cache_value) == bucket_text 
end

function  test_pull_with_cache()
    key = "somekey"

    bucket_text = "mock contents"
    bucket_file = IOBuffer(bucket_text)
    seekstart(bucket_file)
    bucket_files = Dict()
    bucket_files[key] = bucket_file

    cache_text = "mock contents cached already"
    cache_io = IOBuffer(cache_text)
    seekstart(cache_io)
    cache_values = Dict()
    cache_values[key] = cache_io

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = DataSource.Service(bucket, cache)

    DataSource.pull!(datasource, key)

    @test haskey(cache.mockValues, key)
    new_cache_io = cache.mockValues[key]
    @test readchomp(new_cache_io) == cache_text
end
function  test_pull_no_key()
end
function  test_pull_force_empty_cache()
end
function  test_pull_force_with_cache()
end
function  test_push()
end
function  test_push_not_exist()
end
function __init__()
    test_pull_empty_cache()
    test_pull_with_cache()
    test_pull_no_key()
    test_pull_force_empty_cache()
    test_pull_force_with_cache()

    test_push()
    test_push_not_exist()
end

end # module TestDataSource
