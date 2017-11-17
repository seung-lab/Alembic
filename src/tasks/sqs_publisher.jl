using Alembic

using AWSCore
using AWSSQS
using JSON

function get_queue(queue_name)
    secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
    secret_path = joinpath(homedir(), secrets_subpath)
    aws_secrets = JSON.parsefile(secret_path)
    env = aws_config(creds=AWSCredentials(collect(values(aws_secrets))...))
    return sqs_get_queue(env, queue_name)
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
    for (z_start, z_stop) in pairs
        println("Publishing $((z_start, z_stop))")
        p = deepcopy(params)
        p["mesh"]["z_start"] = z_start
        p["mesh"]["z_stop"] = z_stop
        sqs_send_message(q, JSON.json(p))
    end
end

if length(ARGS) > 0
    publish_messages(ARGS...)
end
