module BucketService

abstract BucketService
type AWSBucketService <: BucketService
    name::ASCIIString
    env::AWS.AWSEnv

    function AWSBucketService(env::AWSBucketService, name::ASCIIString,
        base_directory::ASCIIString)
        this = new()
        this.env = env
        this.name = name

        # Check to make sure bucket is reachable
        check_reachable(env, name)
    end
end

function check_reachable(env::AWS.AWSEnv, bucket_name::AbstractString)
    bucket_response = AWS.S3.get_bkt(env, bucket_name)

    if bucket_response.http_code != 200
        error("Unable to access bucket: $bucket_name, response:
            ($(bucket_response.http_code))")
    end

    #= commentting this out. base directory *should* be specified in task
    #parameters
     =base_directory_response = AWS.S3.get_object(env, bucket_name, "$base_directory/")
     =if base_directory_response.http_code != 200
     =    error("Unable to locate base directory: $base_directory, response:
     =        ($(base_directory_response.http_code))")
     =end
     =#
end

end # end module BucketService
