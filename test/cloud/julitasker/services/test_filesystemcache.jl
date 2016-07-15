module TestFileSystemCache

using Base.Test
using CloudTest.JulitasksTests.Utils.TestServices

import Julimaps.Cloud.Julitasks.Services.Cache
import Julimaps.Cloud.Julitasks.Services.FileSystemCache

function test_blank_directory()
    @test_throws ArgumentError FileSystemCache.FileSystemCacheService(" ")
end

function test_good_cache()
    @test FileSystemCache.FileSystemCacheService(TEST_BASE_DIRECTORY) !=
        nothing
end

function test_sanitize_key()
    @test FileSystemCache.sanitize("../../.ssh/id_rsa") == "//.ssh/id_rsa"
end

function test_to_filename()
    cache = make_valid_file_system_cache()
    @test FileSystemCache.to_filename(cache, "file") ==
        "$TEST_BASE_DIRECTORY/file"
end

function test_to_filename_sanitized()
    cache = make_valid_file_system_cache()
    @test FileSystemCache.to_filename(cache, "../../.ssh/id_rsa") ==
        "$TEST_BASE_DIRECTORY///.ssh/id_rsa"
end

function test_exists()
    filename = "testfile"
    full_path_file_name = "$TEST_BASE_DIRECTORY/$filename"
    file = open(full_path_file_name, "w")
    close(file)
    cache = make_valid_file_system_cache()
    @test Cache.exists(cache, filename)
    rm(full_path_file_name)
end

function test_not_exist()
    filename = "testfileDOES NOT EXIST"
    full_path_file_name = "$TEST_BASE_DIRECTORY/$filename"
    cache = make_valid_file_system_cache()
    @test !Cache.exists(cache, filename)
end

function test_put()
    filename = "key"
    test_string = "TEST"
    buffer = IOBuffer(test_string)
    cache = make_valid_file_system_cache()
    Cache.put!(cache, filename, buffer)
    file = open("$TEST_BASE_DIRECTORY/$filename", "r")
    @test isreadable(file)
    @test readall(file) == test_string
    close(file)
end

function test_put_path()
    filename = "a/b/c/key"
    test_string = "TEST"
    buffer = IOBuffer(test_string)
    cache = make_valid_file_system_cache()
    Cache.put!(cache, filename, buffer)
    file = open("$TEST_BASE_DIRECTORY/$filename", "r")
    @test isreadable(file)
    @test readall(file) == test_string
    @test ispath("$TEST_BASE_DIRECTORY/a/b/c/")
    rm("$TEST_BASE_DIRECTORY/$filename")
    rm("$TEST_BASE_DIRECTORY/a/b/c/")
    rm("$TEST_BASE_DIRECTORY/a/b/")
    rm("$TEST_BASE_DIRECTORY/a/")
    close(file)
end

function test_get()
    filename = "a/b/c/key"
    mkpath("$TEST_BASE_DIRECTORY/a/b/c/")
    file = open("$TEST_BASE_DIRECTORY/$filename", "w")
    test_string = "TEST"
    write(file, "TEST")
    close(file)

    cache = make_valid_file_system_cache()
    buffer = Cache.get(cache, filename)
    @test readall(buffer) == test_string
    close(buffer)
end

function test_get_not_exists()
    filename = "testfileDOES NOT EXIST"
    cache = make_valid_file_system_cache()
    buffer = Cache.get(cache, filename)
    @test buffer == nothing
end

function test_delete()
    filename = "key"
    test_string = "TEST"
    buffer = IOBuffer(test_string)
    cache = make_valid_file_system_cache()
    Cache.put!(cache, filename, buffer)
    Cache.delete!(cache, filename)
    @test !isfile("$TEST_BASE_DIRECTORY/$filename")
end

function test_delete_path()
    filename = "a/b/c/key"
    test_string = "TEST"
    buffer = IOBuffer(test_string)
    cache = make_valid_file_system_cache()
    Cache.put!(cache, filename, buffer)
    Cache.delete!(cache, filename)
end

function __init__()
    test_blank_directory()
    test_good_cache()

    test_sanitize_key()
    test_to_filename()
    test_to_filename_sanitized()

    test_exists()
    test_not_exist()

    test_put()
    test_put_path()

    test_get()
    test_get_not_exists()

    test_delete()
    test_delete_path()
end

end # module TestFileSystemCache
