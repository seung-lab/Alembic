module TestDatasource

using Base.Test
using CloudTest.JulitasksTests.Utils.MockServices
using Julimaps.Cloud.Julitasks.Services.BucketCacheDatasource

import Julimaps.Cloud.Julitasks.Services.Bucket
import Julimaps.Cloud.Julitasks.Services.Cache
import Julimaps.Cloud.Julitasks.Services.Datasource

function test_get_empty_cache()
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

    get_value = Datasource.get(datasource, key)

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    # we are returned the correct value
    @test new_cache_value == get_value
    # cache gets updated with the bucket text
    @test readchomp(new_cache_value) == bucket_text 
end

function test_get_multi_empty_cache()
    key = "somekey"
    bucket_text = "mock contents"
    bucket_file = IOBuffer(bucket_text)
    seekstart(bucket_file)

    key2 = "somekey2"
    bucket_text2 = "mock contents2"
    bucket_file2 = IOBuffer(bucket_text2)
    seekstart(bucket_file2)

    bucket_files = Dict()
    bucket_files[key] = bucket_file
    bucket_files[key2] = bucket_file2

    cache_values = Dict()

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)

    datasource = BucketCacheDatasourceService(bucket, cache)

    (get_value1, get_value2) = Datasource.get(datasource, [key, key2])

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    # we are returned the correct value
    @test get_value1 == new_cache_value
    # cache gets updated with the bucket text
    @test readchomp(new_cache_value) == bucket_text

    @test haskey(cache.mockValues, key2)
    new_cache_value2 = cache.mockValues[key2]
    # we are returned the correct value
    @test get_value2 == new_cache_value2
    # cache gets updated with the bucket text
    @test readchomp(new_cache_value2) == bucket_text2
end

function test_get_with_cache()
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

    get_value = Datasource.get(datasource, key)

    @test haskey(cache.mockValues, key)
    new_cache_io = cache.mockValues[key]
    # we are returned the correct value
    @test get_value == new_cache_io
    # cache does not get updated because it is already in there
    @test readchomp(new_cache_io) == cache_text
end

function test_get_no_bucket()
    key = "somekey"

    bucket_files = Dict()
    cache_values = Dict()

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = BucketCacheDatasourceService(bucket, cache)

    @test_throws Exception Datasource.get(datasource, key)
end

function test_override_cache_get_empty_cache()
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

    get_value = Datasource.get(datasource, key; override_cache=true)

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    # we are returned the correct value
    @test get_value == new_cache_value
    # cache gets forced to update with the new bucket text
    @test readchomp(new_cache_value) == bucket_text 
end

function test_override_cache_get_with_cache()
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

    get_value = Datasource.get(datasource, key; override_cache=true)

    @test haskey(cache.mockValues, key)
    new_cache_io = cache.mockValues[key]
    # we are returned the correct value
    @test get_value == new_cache_io
    # cache gets forced to update with the new bucket text
    @test readchomp(new_cache_io) == bucket_text
end

function test_put_no_value()
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

    result = Datasource.put!(datasource, key)

    @test result

    @test haskey(bucket.mockFiles, key)
    new_bucket_file = bucket.mockFiles[key]
    # putting back to bucket with new cached file updates bucket
    @test readchomp(new_bucket_file) == cache_text
end

function test_put_multi_no_value()
    key1 = "somekey1"

    bucket_text1 = "mock contents"
    bucket_file1 = IOBuffer(bucket_text1)
    seekstart(bucket_file1)

    key2 = "somekey2"

    bucket_text2 = "mock contents2"
    bucket_file2 = IOBuffer(bucket_text2)
    seekstart(bucket_file2)

    bucket_files = Dict()
    bucket_files[key1] = bucket_file1
    bucket_files[key2] = bucket_file2

    cache_text1 = "mock contents cached already1"
    cache_io1 = IOBuffer(cache_text1)
    seekstart(cache_io1)

    cache_text2 = "mock contents cached already2"
    cache_io2 = IOBuffer(cache_text2)
    seekstart(cache_io2)

    cache_values = Dict()
    cache_values[key1] = cache_io1
    cache_values[key2] = cache_io2

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = BucketCacheDatasourceService(bucket, cache)

    (result1, result2) = Datasource.put!(datasource, [key1, key2])

    @test result1
    @test result2

    @test haskey(bucket.mockFiles, key1)
    new_bucket_file1 = bucket.mockFiles[key1]
    # putting back to bucket with new cached file updates bucket
    @test readchomp(new_bucket_file1) == cache_text1

    @test haskey(bucket.mockFiles, key2)
    new_bucket_file2 = bucket.mockFiles[key2]
    # putting back to bucket with new cached file updates bucket
    @test readchomp(new_bucket_file2) == cache_text2
end

function test_put_not_exist_no_value()
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

    result = Datasource.put!(datasource, key)

    @test !result

    @test haskey(bucket.mockFiles, key)
    new_bucket_file = bucket.mockFiles[key]
    # putting back to bucket when not existing in cache does not modify bucket
    @test readchomp(new_bucket_file) == bucket_text
end

function test_put_new_value()
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

    new_text = "new value"
    new_value = IOBuffer(new_text)
    seekstart(new_value)
    result = Datasource.put!(datasource, key, new_value)

    @test result == true

    @test haskey(bucket.mockFiles, key)
    new_bucket_file = bucket.mockFiles[key]
    seekstart(new_bucket_file)
    # putting with new value changes both bucket and cache
    @test readchomp(new_bucket_file) == new_text

    @test haskey(cache.mockValues, key)
    new_cache_value = cache.mockValues[key]
    seekstart(new_cache_value)
    # putting with new value changes both bucket and cache
    @test readchomp(new_cache_value) == new_text
end

function test_put_multi_new_value()
    key1 = "somekey1"

    bucket_text1 = "mock contents"
    bucket_file1 = IOBuffer(bucket_text1)
    seekstart(bucket_file1)

    key2 = "somekey2"

    bucket_text2 = "mock contents2"
    bucket_file2 = IOBuffer(bucket_text2)
    seekstart(bucket_file2)

    bucket_files = Dict()
    bucket_files[key1] = bucket_file1
    bucket_files[key2] = bucket_file2

    cache_text1 = "mock contents cached already1"
    cache_io1 = IOBuffer(cache_text1)
    seekstart(cache_io1)

    cache_text2 = "mock contents cached already2"
    cache_io2 = IOBuffer(cache_text2)
    seekstart(cache_io2)

    cache_values = Dict()
    cache_values[key1] = cache_io1
    cache_values[key2] = cache_io2

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = BucketCacheDatasourceService(bucket, cache)

    new_text1 = "new value1"
    new_value1 = IOBuffer(new_text1)
    seekstart(new_value1)

    new_text2 = "new value2"
    new_value2 = IOBuffer(new_text2)
    seekstart(new_value2)

    (result1, result2) = Datasource.put!(datasource, [key1, key2],
        [new_value1, new_value2])

    @test result1
    @test result2

    @test haskey(bucket.mockFiles, key1)
    new_bucket_file1 = bucket.mockFiles[key1]
    seekstart(new_bucket_file1)
    # putting with new value changes both bucket and cache
    @test readchomp(new_bucket_file1) == new_text1

    @test haskey(cache.mockValues, key1)
    new_cache_io1 = cache.mockValues[key1]
    seekstart(new_cache_io1)
    # putting with new value changes both bucket and cache
    @test readchomp(new_cache_io1) == new_text1

    @test haskey(bucket.mockFiles, key2)
    new_bucket_file2 = bucket.mockFiles[key2]
    seekstart(new_bucket_file2)
    # putting with new value changes both bucket and cache
    @test readchomp(new_bucket_file2) == new_text2

    @test haskey(cache.mockValues, key2)
    new_cache_io2 = cache.mockValues[key2]
    seekstart(new_cache_io2)
    # putting with new value changes both bucket and cache
    @test readchomp(new_cache_io2) == new_text2
end

function test_put_not_exist_new_value()
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

    new_text = "new value"
    new_value = IOBuffer(new_text)
    seekstart(new_value)
    result = Datasource.put!(datasource, key, new_value)

    @test result

    @test haskey(bucket.mockFiles, key)
    new_bucket_file = bucket.mockFiles[key]
    # even if the value doesn't exist in cache, the new_text should be used
    @test readchomp(new_bucket_file) == new_text
end

function test_put_multi_new_value_only_cache()
    key1 = "somekey1"

    bucket_text1 = "mock contents"
    bucket_file1 = IOBuffer(bucket_text1)
    seekstart(bucket_file1)

    key2 = "somekey2"

    bucket_text2 = "mock contents2"
    bucket_file2 = IOBuffer(bucket_text2)
    seekstart(bucket_file2)

    bucket_files = Dict()
    bucket_files[key1] = bucket_file1
    bucket_files[key2] = bucket_file2

    cache_text1 = "mock contents cached already1"
    cache_io1 = IOBuffer(cache_text1)
    seekstart(cache_io1)

    cache_text2 = "mock contents cached already2"
    cache_io2 = IOBuffer(cache_text2)
    seekstart(cache_io2)

    cache_values = Dict()
    cache_values[key1] = cache_io1
    cache_values[key2] = cache_io2

    bucket = MockBucketService(bucket_files)
    cache = MockCacheService(cache_values)
    datasource = BucketCacheDatasourceService(bucket, cache)

    new_text1 = "new value1"
    new_value1 = IOBuffer(new_text1)
    seekstart(new_value1)

    new_text2 = "new value2"
    new_value2 = IOBuffer(new_text2)
    seekstart(new_value2)

    (result1, result2) = Datasource.put!(datasource, [key1, key2],
        [new_value1, new_value2]; only_cache=true)

    @test result1
    @test result2

    @test haskey(bucket.mockFiles, key1)
    new_bucket_file1 = bucket.mockFiles[key1]
    new_cache_io1 = cache.mockValues[key1]
    # putting with new value changes both bucket and cache
    @test readchomp(new_bucket_file1) == bucket_text1
    @test readchomp(new_cache_io1) == new_text1

    @test haskey(bucket.mockFiles, key2)
    new_bucket_file2 = bucket.mockFiles[key2]
    new_cache_io2 = cache.mockValues[key2]
    # putting with new value changes both bucket and cache
    @test readchomp(new_bucket_file2) == bucket_text2
    @test readchomp(new_cache_io2) == new_text2

end

function __init__()
    test_get_empty_cache()
    test_get_multi_empty_cache()
    test_get_with_cache()
    test_get_no_bucket()
    test_override_cache_get_empty_cache()
    test_override_cache_get_with_cache()

    test_put_no_value()
    test_put_multi_no_value()
    test_put_not_exist_no_value()

    test_put_new_value()
    test_put_multi_new_value()
    test_put_not_exist_new_value()
end

end # module TestDatasource
