#export montageSection, montageSections, get_mesh_set, set_mesh_set, sync_mesh_set

#using Images
#using HDF5
#using JLD
#include("convolve.jl")
#include("Mesh.jl")
#include("Matches.jl")
#include("MeshSet.jl")
#include("MeshSolve.jl")

function montageSection(n)
	@time Ms, imageArray = loadSection(SESSION, n);
	@time addAllMatches!(Ms, imageArray);
	@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
	printResidualStats(Ms);
	save(Ms);
	imageArray = 0;
	gc();
end


function align_stack(wafer_num, k::UnitRange{Int64})
	@time Ms, imageArray = load_stack(wafer_num, k);
	@time addAllMatches!(Ms, imageArray);
	@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
	printResidualStats(Ms);
	save(Ms);
	return Ms, imageArray
end

@everywhere function get_mesh_set()
	return Ms;
end

@everywhere function sync_mesh_set()
	id = myid()
	@sync @parallel for i in 1:num_procs
		if i == id continue; end
		 remotecall(i, set_mesh_set(id));
	end
end

@everywhere function set_mesh_set(id)
	Ms = remotecall_fetch(id, get_mesh_set);
	return Ms;
end

#=
function montageSections(ran::UnitRange{Int64})
@time for k in minimum(ran)-1:num_concurrent:maximum(ran)
	toFetch = @sync @parallel for l in 1:num_concurrent
	ind = k + l;
	println(ind);
	if ind > maximum(ran) return; end
	println(ind);
	@time montageSection(ind);
	end
	end
end
=#




################################### DIAGNOSTIC #########################################


end #END MESHMODULE
################################### SCRIPT ############################################

