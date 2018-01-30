using JSON
include("tasks.jl")

param_file = ARGS[1]

addprocs();
using Alembic

params = JSON.parsefile(param_file)
println("Runnign render task...")
@time begin
   match_task(params)
end
