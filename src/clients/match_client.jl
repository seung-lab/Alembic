using Alembic

param_file = ARGS[1]
load_params(param_file)
ms = make_stack()
match!(ms)

obj_name = "meshset"
s = Alembic.StorageWrapper(get_path(obj_name))
s[get_name(ms)] = ms;
