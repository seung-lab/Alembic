#using Julimaps

#include("Mesh.jl")
#include("convolve.jl")

type Matches
	src_index::Index					# source mesh index
	dst_index::Index					# destination mesh index

	n::Int64						# number of matches

	src_pointIndices::Array{Int64, 1}			# index of points in src_mesh.points
	dst_points::Points					# location of points in the destination
	dst_triangles::Triangles				# index of the triangles
	dst_weights::Weights					# weights at each

	dispVectors::Points					# displacement vector src->dest
end

# i, j are in world coordinates (as per the mesh coordinate specification)
function getBlockMatchAtPoint(A, A_rr::RemoteRef, i, j, B, B_rr::RemoteRef, block_size, search_r)
Am = take!(A_rr);
Bm = take!(B_rr);
return getBlockMatchAtPoint(A, Am, i, j, B, Bm, block_size, search_r);
end

function getBlockMatchAtPoint(A, Am::Mesh, i, j, B, Bm::Mesh, block_size, search_r)
	b_rad = block_size + search_r;

	# convert to matrix coordinates
	Ai = round(Int64, i - Am.disp[1]);
	Aj = round(Int64, j - Am.disp[2]);

	Bi = round(Int64, i - Bm.disp[1]);
	Bj = round(Int64, j - Bm.disp[2]);

	if (!isInternal(A, Ai, Aj, block_size) || !isInternal(B, Bi, Bj, b_rad))
		return noMatch;
	end

	Ai_range = ceil(Ai)-block_size:ceil(Ai)+block_size;
	Aj_range = ceil(Aj)-block_size:ceil(Aj)+block_size;
	Bi_range = ceil(Bi)-b_rad:ceil(Bi)+b_rad;
	Bj_range = ceil(Bj)-b_rad:ceil(Bj)+b_rad;

	xc = normxcorr2(sub(A, Ai_range, Aj_range), sub(B, Bi_range, Bj_range));
	r_max = maximum(xc); 

	ind = findfirst(r_max .== xc);

	if ind == 0 return noMatch; end
	(i_max, j_max) = (rem(ind, size(xc, 1)), cld(ind, size(xc, 1)));
	if i_max == 0 i_max = size(xc, 1); end

	return [i_max - 1 - search_r; j_max - 1 - search_r; r_max];
	
end

function Meshes2Matches(A, A_rr::RemoteRef, B, B_rr::RemoteRef, block_size, search_r, min_r)
Am = take!(A_rr);
Bm = take!(B_rr);
return Meshes2Matches(A, Am, B, Bm, block_size, search_r, min_r)
end

function Meshes2Matches(A, Am::Mesh, B, Bm::Mesh, block_size, search_r, min_r)

	src_index = Am.index;
	dst_index = Bm.index;
	p1 = Am.name;
	p2 = Bm.name;

	n = 0;
	n_total = 0;
	n_lowr = 0;
	n_outlier = 0;
	n_noTriangle = 0;
	n_upperbound = Am.n;

	src_pointIndices = Array(Int64, 0);
	dst_points = Points(0);
	dst_triangles = Triangles(0);
	dst_weights = Weights(0);
	dispVectors = Points(0);
	dispVectors_raw = Array{Array{Float64, 1}, 1}(0);
	dispVectors_mags = Array{Float64, 1}(0);

	if (Am==Bm)
		return Void;
	end

	#=if (is_pre_aligned(Am.index) && is_pre_aligned(Bm.index))
	j = 1;
	nextidx() = (idx=j; j+=1; idx);
	images = SharedArray(UInt8, size(A, 1), size(A, 2), 2);
	images[:, :, 1] = A;
	images[:, :, 2] = B;
	dispVectors_raw_par = Array{Array{Float64, 1}, 1}(n_upperbound);
	@sync begin
	for p in 1:num_procs
		if p != myid() && isodd(p-myid()) || num_procs == 1
			@async begin
				while true
					idx = nextidx();
						if idx > n_upperbound
							break
						end
					(Ai, Aj) = Am.nodes[j];
					A_rr = RemoteRef();
					B_rr = RemoteRef();
					put!(A_rr, Am);
					put!(B_rr, Bm);
					dispVectors_raw_par[j] = remotecall_fetch(p, getBlockMatchAtPoint, images[:, :, 1], A_rr, Ai, Aj, images[:, :, 2], B_rr, block_size, search_r);
					#not in correct order
					
					
				end
			end
		end
	end
	end
	for j in 1:n_upperbound
	v = dispVectors_raw_par[j];
	push!(dispVectors_raw, v);
	if v != noMatch && v[3] >= min_r
	push!(dispVectors_mags, norm(v));
	end
	end
	else=#
	for j in 1:n_upperbound
		(Ai, Aj) = Am.nodes[j];
		v = getBlockMatchAtPoint(A, Am, Ai, Aj, B, Bm, block_size, search_r);
		push!(dispVectors_raw, v);
		if v != noMatch && v[3] >= min_r
			push!(dispVectors_mags, norm(v));
		end
	end

	#end

	if length(dispVectors_mags) == 0
	return Void;
	end


	mu = mean(dispVectors_mags);
	sig = std(dispVectors_mags);
	max = maximum(dispVectors_mags);

	for j in 1:n_upperbound
		v = dispVectors_raw[j];	
		if v == noMatch continue; end	
		n_total +=1	
		if v[3] < min_r; n_lowr +=1; continue; end
		dispVector = v[1:2];
		if norm(dispVector) > mu + 2.5 * sig; n_outlier +=1; continue; end
		dst_point = Am.nodes[j] + dispVector;
		dst_triangle = findMeshTriangle(Bm, dst_point[1], dst_point[2]); 
		if dst_triangle == noTriangle n_noTriangle +=1; continue; end
		n += 1;
		push!(src_pointIndices, j);
		push!(dispVectors, dispVector);
		push!(dst_points, Am.nodes[j] + dispVectors[n]);
		push!(dst_triangles, dst_triangle);
		push!(dst_weights, getTriangleWeights(Bm, dst_triangle, dst_point[1], dst_point[2]));
		#=if !isnan(sum(xc))
		#	imwrite(xcorr2Image(xc), joinpath(".","output_images", "normxcorr", string(join(Am.index, "_"), "_", join(Bm.index, "_"), "_", n, ".jpg")))
		#end
		=#
	end

	if n == 0
	return Void;
	end

	println("###\n$p1 -> $p2: $n_upperbound in mesh, $n_total in overlap, $n accepted.\nRejections: $n_lowr (low r), $n_outlier (outliers), $n_noTriangle (outside triangles).\nDisplacement mean = $mu, sigma = $sig, max = $max\n###");

	matches = Matches(src_index, dst_index, n, src_pointIndices, dst_points, dst_triangles, dst_weights, dispVectors);
	return matches;
end

