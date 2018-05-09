using JSON
include("tasks.jl")

param_file = ARGS[1]

n = 1
if length(ARGS) > 1
    n = parse(ARGS[2])
end

println("Importing Alembic...")
@time begin
    if n > 1
        addprocs(n);
        @everywhere using Alembic
    else
        using Alembic
    end
end

params = JSON.parsefile(param_file)
println("Running match task...")
@time begin
   match_task(params)
end
