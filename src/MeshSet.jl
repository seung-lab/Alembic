#=using Julimaps
using Params
using HDF5
using JLD
using Images
importall IO
include("Mesh.jl")
include("Matches.jl")
include("MeshSolve.jl")

export save, load

=#

type MeshSet
	N::Int64						# number of meshes in the set
	M::Int64						# number of matches in the set - (a -> b) and (b -> a) are distinct

	indices::Array{Index, 1}				# wafer, section, row, column as a tuple - if tileIndex is 0 then denotes entire section

	n::Int64						# number of nodes in the set across the whole set
	m::Int64						# number of edges in the set across the whole set
	m_i::Int64						# number of internal edges in the set across the whole set
	m_e::Int64						# number of edges between meshes
	
	meshes::Array{Mesh, 1}					# vector of meshes in the set
	nodes_indices::Array{Int64, 1}				# vector of number of nodes before the n-th mesh. To get index at i-th node of n-th mesh, nodes_indices[n] + i.

	matches::Array{Matches, 1}				# vector of matches in the set
	matches_pairs::Pairings				# vector of index (in meshes) - (a, b) means the match is between (meshes[a] -> meshes[b])
end

function isAdjacent(Am, Bm)
	if abs(Am.index[3] - Bm.index[3]) + abs(Am.index[4] - Bm.index[4]) + abs(Am.index[2] - Bm.index[2]) == 1 return true; end
	return false;
end

function isDiagonal(Am, Bm)
	if abs(Am.index[3] - Bm.index[3]) + abs(Am.index[4] - Bm.index[4]) == 2 && Am.index[3] != Bm.index[3] && Am.index[4] != Bm.index[4] return true; end
	return false;
end

function findNode(Ms, mesh_ind, node_ind);
	return Ms.nodes_indices[mesh_ind] + node_ind;
end

function findIndex(Ms, mesh_index_tuple)
	return findfirst(this -> mesh_index_tuple == this.index, Ms.meshes)
end

function makeNewMeshSet()
	N = 0;
	M = 0;

	indices = Array{Index, 1}(0);
	
	n = 0;
	m = 0;
	m_i = 0;
	m_e = 0;

	meshes = Array{Mesh, 1}(0);
	nodes_indices = Array{Int64, 1}(0);

	matches = Array{Matches, 1}(0);
	matches_pairs = Array{Pair, 1}(0);

	return MeshSet(N, M, indices, n, m, m_i, m_e, meshes, nodes_indices, matches, matches_pairs);
end


function addMesh2MeshSet!(Am, Ms)
	push!(Ms.indices, Am.index);
	push!(Ms.meshes, Am);
	if length(Ms.nodes_indices) == 0 push!(Ms.nodes_indices, 0);
	else push!(Ms.nodes_indices, Ms.n);
	end
	Ms.N += 1;
	Ms.m_i += Am.m;
	Ms.m += Am.m;
	Ms.n += Am.n;
	return;
end

function addMatches2MeshSet!(M, Ms)
	if (typeof(M) == Void || M == Void) return; end
	push!(Ms.matches, M);
	push!(Ms.matches_pairs, (findIndex(Ms, M.src_index), findIndex(Ms, M.dst_index)));
	Ms.M += 1;
	Ms.n;
	Ms.m += M.n;
	Ms.m_e += M.n;
	return;
end

function solveMeshSet!(Ms, match_coeff, eta_gradient, eta_newton, ftol_grad, ftol_newton)
	nodes = Points(0);
	nodes_fixed = BinaryProperty(0);
	edges = spzeros(Float64, Ms.n, 0);
	edge_lengths = FloatProperty(0);
	edge_coeffs = FloatProperty(0);
	for i in 1:Ms.N
		cur_mesh = Ms.meshes[i];
		if i == 1 nodes = hcat(cur_mesh.nodes...);
		else nodes = hcat(nodes, hcat(cur_mesh.nodes...)); end
		append!(nodes_fixed, cur_mesh.nodes_fixed);
		append!(edge_lengths, cur_mesh.edge_lengths);
		append!(edge_coeffs, cur_mesh.edge_coeffs);
		if (i == Ms.N)	
			edges = hcat(edges, vcat(spzeros(Float64, Ms.nodes_indices[i], cur_mesh.m), 
						 cur_mesh.edges));
		else
			edges = hcat(edges, vcat(spzeros(Float64, Ms.nodes_indices[i], cur_mesh.m), 
						 cur_mesh.edges, 
						 spzeros(Float64, Ms.n - Ms.nodes_indices[i] - cur_mesh.n, cur_mesh.m)));
		end
	end

	for i in 1:Ms.M
		cur_matches = Ms.matches[i];
		append!(edge_lengths, fill(0.0, cur_matches.n));
		append!(edge_coeffs, fill(convert(Float64, match_coeff), cur_matches.n));
		edges_padded = spzeros(Float64, Ms.n, cur_matches.n);

		for j in 1:Ms.matches[i].n
			edges_padded[findNode(Ms, Ms.matches_pairs[i][1], cur_matches.src_pointIndices[j]), j] = -1;
			edges_padded[findNode(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][1]), j] = cur_matches.dst_weights[j][1];
			edges_padded[findNode(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][2]), j] = cur_matches.dst_weights[j][2];
			edges_padded[findNode(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][3]), j] = cur_matches.dst_weights[j][3];
		end
		edges = hcat(edges, edges_padded);
	end

	SolveMesh!(nodes, nodes_fixed, edges, edge_coeffs, edge_lengths, eta_gradient, eta_newton, ftol_grad, ftol_newton);
	nodes_t = Points(0);
	for i in 1:size(nodes, 2)
       		push!(nodes_t, vec(nodes[:, i]))
       	end
	for i in 1:Ms.N
		cur_mesh = Ms.meshes[i];
		cur_mesh.nodes_t = nodes_t[Ms.nodes_indices[i] + (1:cur_mesh.n)];
	end

end

function save(filename::String, Ms::MeshSet)
	jldopen(filename, "w") do file
		write(file, "MeshSet", Ms);
	end
end

function save(Ms::MeshSet)
	firstindex = Ms.meshes[1].index;
	lastindex = Ms.meshes[Ms.N].index;

	if is_pre_aligned(firstindex)
filename = joinpath(ALIGNED_DIR, string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","), "_aligned.jld"));
	else
	filename = joinpath(MONTAGE_DIR, string(join(firstindex[1:2], ","), "_montaged.jld"));
	end

	jldopen(filename, "w") do file
		write(file, "MeshSet", Ms);
	end
end

function JLD2MeshSet(filename)
	Ms = load(filename, "MeshSet"); 
	return Ms;
end

function getAllOverlaps(Ms)
adjacent_pairs = Pairings(0);
diagonal_pairs = Pairings(0);

	for i in 1:Ms.N, j in 1:Ms.N
		if isAdjacent(Ms.meshes[i], Ms.meshes[j]) push!(adjacent_pairs, (i, j)); end
		if isDiagonal(Ms.meshes[i], Ms.meshes[j]) push!(diagonal_pairs, (i, j)); end
	end

	pairs = vcat(adjacent_pairs, diagonal_pairs);

	return pairs;
end
#=
function addAllMatches!(Ms, imageArray::SharedArray)

pairs = getAllOverlaps(Ms)

@time for k in 0:num_procs:length(pairs)
	toFetch = @sync @parallel for l in 1:num_procs
	ind = k + l;
	if ind > length(pairs) return; end
	(i, j) = pairs[ind];
	return Meshes2Matches(imageArray[:, :, i], Ms.meshes[i], imageArray[:, :, j], Ms.meshes[j], block_size, search_r, min_r);
	end
	for i = 1:length(toFetch)
		M = fetch(toFetch[i])
		if typeof(M) == Void continue; end
		addMatches2MeshSet!(M, Ms);
	end
	end
	return Ms;
end
=#

function addAllMatches!(Ms, images::SharedArray)

pairs = getAllOverlaps(Ms);
n = length(pairs);
i = 1;
nextidx() = (idx=i; i+=1; idx);
matchesArray = cell(n);

@sync begin
	for p in 1:num_procs
		if p != myid() || num_procs == 1
			@async begin
				while true
					idx = nextidx();
						if idx > n
							break
						end
					(a, b) = pairs[idx];
					A_rr = RemoteRef();
					B_rr = RemoteRef();
					put!(A_rr, Ms.meshes[a]);
					put!(B_rr, Ms.meshes[b]);
					if is_pre_aligned(Ms.meshes[a].index)
matchesArray[idx] = remotecall_fetch(p, Meshes2Matches, images[:, :, a], A_rr, images[:, :, b], B_rr, block_size_alignment, search_r_alignment, min_r_alignment);
					else
					matchesArray[idx] = remotecall_fetch(p, Meshes2Matches, images[:, :, a], A_rr, images[:, :, b], B_rr, block_size, search_r, min_r);
					end
				end
			end
		end
	end
end


for k in 1:n
		M = matchesArray[k]
		if typeof(M) == Void || M == Void continue; end
		addMatches2MeshSet!(M, Ms);
end
	return Ms;
end



#=
function addAllMatches!(Ms, imageArray)#::Array{Array{Float64, 2}, 1})

pairs = getAllOverlaps(Ms)

@time for ind in 1:length(pairs)
	(i, j) = pairs[ind];
	M = Meshes2Matches(imageArray[i], Ms.meshes[i], imageArray[j], Ms.meshes[j], block_size_alignment, search_r_alignment, min_r_alignment);
	if typeof(M) == Void continue; end
	addMatches2MeshSet!(M, Ms);
	end
	return Ms;
end
=#

function makeSectionMeshSet(session, section_num)
	Ms = makeNewMeshSet();
	indices = find(i -> session[i,2][2] == section_num, 1:size(session, 1))
	for i in 1:length(indices)
	name = session[i, 1];
	index = session[i, 2];
	dx = session[i, 3];
	dy = session[i, 4];
	addMesh2MeshSet!(Tile2Mesh(name, index, dy, dx, false, mesh_length, mesh_coeff), Ms);
	end
	return Ms;
end

function make_stack_meshset(wafer_num, section_range)
	Ms = makeNewMeshSet();
	for i in section_range
	name = getName(wafer_num, i);
	index = (wafer_num, i, 0, 0);
	dx = 0;
	dy = 0;
	addMesh2MeshSet!(Tile2Mesh(name, index, dy, dx, false, mesh_length_alignment, mesh_coeff_alignment), Ms);
	end
end

function get_matched_points(Ms::MeshSet, k)

src_p = Points(0);
dst_p = Points(0);

for i in 1:Ms.matches[k].n
		w = Ms.matches[k].dst_weights[i];
		t = Ms.matches[k].dst_triangles[i];
		p = Ms.matches[k].src_pointIndices[i];
		src = Ms.meshes[findIndex(Ms, Ms.matches[k].src_index)]
		dst = Ms.meshes[findIndex(Ms, Ms.matches[k].dst_index)]
		p1 = src.nodes[p];
		p2 = dst.nodes[t[1]] * w[1] + dst.nodes[t[2]] * w[2] + dst.nodes[t[3]] * w[3];
		push!(src_p, p1);
		push!(dst_p, p2);
end
	return src_p, dst_p;
end

function get_matched_points_t(Ms::MeshSet, k)

src_p = Points(0);
dst_p = Points(0);

for i in 1:Ms.matches[k].n
		w = Ms.matches[k].dst_weights[i];
		t = Ms.matches[k].dst_triangles[i];
		p = Ms.matches[k].src_pointIndices[i];
		src = Ms.meshes[findIndex(Ms, Ms.matches[k].src_index)]
		dst = Ms.meshes[findIndex(Ms, Ms.matches[k].dst_index)]
		p1 = src.nodes_t[p];
		p2 = dst.nodes_t[t[1]] * w[1] + dst.nodes_t[t[2]] * w[2] + dst.nodes_t[t[3]] * w[3];
		push!(src_p, p1);
		push!(dst_p, p2);
end
	return src_p, dst_p;
end

function loadSection(session, section_num)
	indices = find(i -> session[i,2][2] == section_num, 1:size(session, 1));
	Ms = makeNewMeshSet();
	num_tiles = length(indices);
	paths = Array{String, 1}(num_tiles);

	imageArray = SharedArray(UInt8, tile_size, tile_size, num_tiles);

	ind = 1;

	for i in indices
		name = session[i, 1];
		index = session[i, 2];
		dx = session[i, 3];
		dy = session[i, 4];
		image = getImage(getPath(name));
		addMesh2MeshSet!(Tile2Mesh(name, image, index, dy, dx, false, mesh_length, mesh_coeff), Ms);
		imageArray[:, :, ind] = image;
		ind+=1;
	end

	return Ms, imageArray;
end

function load_section(session, section_num)
	indices = find(i -> session[i,2][2] == section_num, 1:size(session, 1));
	Ms = makeNewMeshSet();
	num_tiles = length(indices);
	paths = Array{String, 1}(num_tiles);

	for i in indices
		name = session[i, 1];
		index = session[i, 2];
		dx = session[i, 3];
		dy = session[i, 4];
		addMesh2MeshSet!(Tile2Mesh(name, tile_size, index, dy, dx, false, mesh_length, mesh_coeff), Ms);
	end

	return Ms;
end

function load_section_images(session, section_num)
	indices = find(i -> session[i,2][2] == section_num, 1:size(session, 1));
	num_tiles = length(indices);
	paths = Array{String, 1}(num_tiles);

	imageArray = SharedArray(UInt8, tile_size, tile_size, num_tiles);

	ind = 1;

	for i in indices
		name = session[i, 1];
		imageArray[:, :, ind] = getImage(getPath(name));
		ind+=1;
	end

	return imageArray;
end


function load_stack(offsets, wafer_num, section_range)
	indices = find(i -> offsets[i, 2][1] == wafer_num && offsets[i,2][2] in section_range, 1:size(offsets, 1));
	Ms = makeNewMeshSet();
	imageArray = Array{Array{UInt8, 2}, 1}(0);# SharedArray(Int64, tile_size, tile_size, num_tiles);

	i_max = 0;
	j_max = 0;

	for i in indices
	name = offsets[i, 1];
	index = offsets[i, 2];
	dx = offsets[i, 4];
	dy = offsets[i, 3];
	image = getImage(getPath(name));
	addMesh2MeshSet!(Tile2Mesh(name, image, index, dy, dx, false, mesh_length_alignment, mesh_coeff), Ms);
	
	push!(imageArray, image);
	i_l = size(imageArray[i], 1);
	j_l = size(imageArray[i], 2);
	if i_max < i_l i_max = i_l end;
	if j_max < j_l j_max = j_l end;
	end
	

	images = SharedArray(UInt8, i_max, j_max, length(imageArray));

	for i in 1:length(imageArray)
	i_l = size(imageArray[i], 1);
	j_l = size(imageArray[i], 2);
	images[1:i_l, 1:j_l, i] = imageArray[i];
	end

	return Ms, images;
end

function printResidualStats(Ms)
	residuals_t = Points(0);
	for k in 1:Ms.M
		for i in 1:Ms.matches[k].n
			w = Ms.matches[k].dst_weights[i];
			t = Ms.matches[k].dst_triangles[i];
			p = Ms.matches[k].src_pointIndices[i];
			src = Ms.meshes[findIndex(Ms, Ms.matches[k].src_index)]
			dst = Ms.meshes[findIndex(Ms, Ms.matches[k].dst_index)]
			p1 = src.nodes_t[p];
			p2 = dst.nodes_t[t[1]] * w[1] + dst.nodes_t[t[2]] * w[2] + dst.nodes_t[t[3]] * w[3]
			push!(residuals_t, p2-p1);
		end
	end
	res_norm = map(norm, residuals_t);
	avg = mean(res_norm);
	sig = std(res_norm);
	max = maximum(res_norm);
	println("Residuals after solving elastically: mean: $avg, sigma = $sig, max = $max");
end



