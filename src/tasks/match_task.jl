
param_file = ARGS[1]
nprocs = parse(Int64, ARGS[2])

addprocs(nprocs)
@everywhere using Alembic

load_params(param_file)
ms = make_stack()
match!(ms)

obj_name = "meshset"
s = Alembic.StorageWrapper(get_path(obj_name))
s[get_name(ms)] = ms;
