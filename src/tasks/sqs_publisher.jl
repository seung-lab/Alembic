using Alembic

using AWSCore
using AWSSQS
using JSON

function publish_messages(params_fn, queue_name)
    secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
    secret_path = joinpath(homedir(), secrets_subpath)
    aws_secrets = JSON.parsefile(secret_path)
    env = aws_config(creds=AWSCredentials(collect(values(aws_secrets))...))
    q = sqs_get_queue(env, queue_name)

    load_params_file(params_fn)
    pairs = get_all_overlaps(get_z_range(), PARAMS[:match][:depth], 
                                        symmetric=PARAMS[:match][:symmetric])
    params = JSON.parsefile(params_fn)
    for (z_start, z_stop) in pairs
        println("Publishing $((z_start, z_stop))")
        p = deepcopy(params)
        p["mesh"]["z_start"] = z_start
        p["mesh"]["z_stop"] = z_stop
        sqs_send_message(q, JSON.json(p))
    end
end

function publish_messages(params_fn, queue_name, pairs_fn)
    secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
    secret_path = joinpath(homedir(), secrets_subpath)
    aws_secrets = JSON.parsefile(secret_path)
    env = aws_config(creds=AWSCredentials(collect(values(aws_secrets))...))
    q = sqs_get_queue(env, queue_name)

    pairs_arr = readdlm(pairs_fn, ',', Int)
    pairs = [[pairs_arr[i,:]...] for i in 1:size(pairs_arr,1)]
    params = JSON.parsefile(params_fn)

    for (z_start, z_stop) in pairs
        println("Publishing $((z_start, z_stop))")
        p = deepcopy(params)
        p["mesh"]["z_start"] = z_start
        p["mesh"]["z_stop"] = z_stop
        sqs_send_message(q, JSON.json(p))
    end
end

publish_messages(ARGS...)
