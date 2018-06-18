using Alembic

z_map = Dict()
adj_z = 0
for z in 1:2176
	if z == 465
		z_map[z] = 468
		continue
	end
	if z == 466
		adj_z -= 1
	end
	if z == 469
		adj_z += 1
	end
	if z == 535
		adj_z += 1
	end
	if z == 542
		adj_z += 1
	end
	if z == 549
		adj_z += 1
	end
	if z == 618
		adj_z += 1
	end
	if z == 788
		adj_z += 1
	end
	if z == 792
		adj_z += 1
	end
	if z == 853
		adj_z += 1
	end
	if z == 984
		adj_z += 1
	end
	if z == 1041
		adj_z += 1
	end
	z_map[z] = z+adj_z
end

fn = joinpath(homedir(), ".julia/v0.6/Alembic/src/params/pinky100_v1/180605_z_map.csv")
writedlm(fn, z_map, ",")

cv_path = "gs://neuroglancer/pinky100_v0/image_single_slices"
save(z_map, cv_path, "z_map")