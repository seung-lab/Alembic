module TestDatasource

using Base.Test
using CloudTest.MockServices
using Julimaps.Cloud.Services.BucketCacheDatasource

import Julimaps.Cloud.Services.Bucket
import Julimaps.Cloud.Services.Cache
import Julimaps.Cloud.Services.Datasource

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

    datasource = BucketCacheDatasourceService(bucket, cache)

    Datasource.pull!(datasource, key)

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    # cache gets updated with the bucket text
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
    datasource = BucketCacheDatasourceService(bucket, cache)

    Datasource.pull!(datasource, key)

    @test haskey(cache.mockValues, key)
    new_cache_io = cache.mockValues[key]
    # cache does not get updated because it is already in there
    @test readchomp(new_cache_io) == cache_text
end

function  test_pull_no_bucket()
    key = "somekey"

    bucket_files = Dict()
    cache_values = Dict()

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = BucketCacheDatasourceService(bucket, cache)

    @test_throws Exception Datasource.pull!(datasource, key)
end

function test_force_pull_empty_cache()
    key = "somekey"

    bucket_text = "mock contents"
    bucket_file = IOBuffer(bucket_text)
    seekstart(bucket_file)
    bucket_files = Dict()
    bucket_files[key] = bucket_file

    cache_values = Dict()

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)

    datasource = BucketCacheDatasourceService(bucket, cache)

    Datasource.pull!(datasource, key; force=true)

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    # cache gets forced to update with the new bucket text
    @test readchomp(new_cache_value) == bucket_text 
end

function  test_force_pull_with_cache()
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
    datasource = BucketCacheDatasourceService(bucket, cache)

    Datasource.pull!(datasource, key; force=true)

    @test haskey(cache.mockValues, key)
    new_cache_io = cache.mockValues[key]
    # cache gets forced to update with the new bucket text
    @test readchomp(new_cache_io) == bucket_text
end

function test_push()
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
    datasource = BucketCacheDatasourceService(bucket, cache)

    result = Datasource.push!(datasource, key)

    @test result == true

    @test haskey(bucket.mockFiles, key)
    bucket_file = bucket.mockFiles[key]
    # pushing back to bucket with new cached file updates bucket
    @test readchomp(bucket_file) == cache_text
end

function  test_push_not_exist()
    key = "somekey"

    bucket_text = "mock contents"
    bucket_file = IOBuffer(bucket_text)
    seekstart(bucket_file)
    bucket_files = Dict()
    bucket_files[key] = bucket_file

    cache_values = Dict()

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = BucketCacheDatasourceService(bucket, cache)

    result = Datasource.push!(datasource, key)

    @test result == false

    @test haskey(bucket.mockFiles, key)
    bucket_file = bucket.mockFiles[key]
    # pushing back to bucket when not existing in cache does not modify bucket
    @test readchomp(bucket_file) == bucket_text
end

function __init__()
    test_pull_empty_cache()
    test_pull_with_cache()
    test_pull_no_bucket()
    test_force_pull_empty_cache()
    test_force_pull_with_cache()

    test_push()
    test_push_not_exist()
end

end # module TestDatasource
