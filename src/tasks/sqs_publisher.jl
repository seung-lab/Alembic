using Alembic

using AWSCore
using AWSSQS
using JSON

function get_queue(queue_name)
    # secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
    # secret_path = joinpath(homedir(), secrets_subpath)
    # aws_secrets = JSON.parsefile(secret_path)
    # env = aws_config(creds=AWSCredentials(collect(values(aws_secrets))...))
    # return sqs_get_queue(env, queue_name)
    secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
    secret_path = joinpath(homedir(), secrets_subpath)
    aws_secrets = JSON.parsefile(secret_path)
    env = aws_config(creds=AWSCredentials(collect(values(aws_secrets))...))
    # q = sqs_get_queue(env, queue_name)
    # This is a hack, using https://github.com/jingpengw/AWSCore.jl
    # + expanding the sqs_get_queue function. Should fork AWSCore
    # for Alembic.
    r = AWSSQS.sqs(env, "GetQueueUrl", QueueName = queue_name)
    url = r["QueueUrl"]
    return merge(env, Dict(:resource => HTTP.URI(url).path))
end

function publish_message(params_fn::AbstractString, queue_name::AbstractString, 
                                                    pairs_fn::AbstractString)
    q = get_queue(queue_name)
    load_params(params_fn)
    params = JSON.parsefile(params_fn)
    pairs = get_all_overlaps(get_z_range(), PARAMS[:match][:depth], 
                                        symmetric=PARAMS[:match][:symmetric])
    if params_fn != ""
        pairs_arr = readdlm(pairs_fn, ',', Int)
        pairs = [[pairs_arr[i,:]...] for i in 1:size(pairs_arr,1)]
    end

    return publish_messages(q, params, pairs)
end

function publish_messages(queue_name::AbstractString, params::Dict, pairs::Array)
    q = get_queue(queue_name)
    for p in pairs
        name = string(params["task"]["method"], " for ", p)
        println("Publishing $name")
        params_copy = deepcopy(params)
        params_copy["task"]["name"] = name
        params_copy["task"]["pairs"] = p
        sqs_send_message(q, JSON.json(params_copy))
    end
end

if length(ARGS) > 0
    publish_messages(ARGS...)
end
