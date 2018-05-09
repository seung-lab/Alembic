using JSON
include("tasks.jl")

param_file = ARGS[1]

addprocs();
using Alembic

params = JSON.parsefile(param_file)
println("Running render task...")
@time begin
   render_task(params)
end
