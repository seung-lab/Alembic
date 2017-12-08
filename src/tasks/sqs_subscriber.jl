addprocs()
@time @everywhere using Alembic

using AWSCore
using AWSSQS
using JSON

function receive_messages(queue_name)
    secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
    secret_path = joinpath(homedir(), secrets_subpath)
    aws_secrets = JSON.parsefile(secret_path)
    env = aws_config(creds=AWSCredentials(collect(values(aws_secrets))...))
    q = sqs_get_queue(env, queue_name)
    while true
        m = sqs_receive_message(q)
        if m != nothing
            params = JSON.parse(m[:message])
            name = (params["mesh"]["z_start"], params["mesh"]["z_stop"])
            println("Received $name")
            func = getfield(Main, Symbol(params["task_method"]))
            func(params)
            sqs_delete_message(q, m)
        else
            sleep(5)
        end
    end
end

receive_messages(ARGS[1])
