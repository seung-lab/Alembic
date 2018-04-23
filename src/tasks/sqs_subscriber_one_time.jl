addprocs()
@time using Alembic

using AWSCore
using AWSSQS
using JSON

function receive_message(queue_name)
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
    q = merge(env, Dict(:resource => HTTP.URI(url).path))
    t = 0
    while t == 0
        m = sqs_receive_message(q)
        if m != nothing
            params = JSON.parse(m[:message])
            task = params["task"]
            println("Received $(task["name"])")
            func = getfield(Main, Symbol(task["method"]))
            @time func(params)
            sqs_delete_message(q, m)
            t += 1
        else
            sleep(5)
        end
    end
end

receive_message(ARGS[1])
