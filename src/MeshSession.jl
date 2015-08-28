
function montageSection(n)
	@time Ms, imageArray = loadSection(PRE_MONTAGED_OFFSETS, n);
	@time addAllMatches!(Ms, imageArray);
	@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
	printResidualStats(Ms);
	save(Ms);
	imageArray = 0;
	gc();
end


function align_stack(wafer_num, k::UnitRange{Int64})
	@time Ms, imageArray = load_stack(PRE_ALIGNED_OFFSETS, wafer_num, k);
	@time addAllMatches!(Ms, imageArray);
	@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
	printResidualStats(Ms);
	save(Ms);
	return Ms, imageArray
end


#=

### show matches
function show_matches(Ms, k)
src_p = Points(0);
dst_p = Points(0);
k = 1
for i in 1:Ms.matches[k].n
		w = Ms.matches[k].dst_weights[i];
		t = Ms.matches[k].dst_triangles[i];
		p = Ms.matches[k].src_pointIndices[i];
		src = Ms.meshes[MeshModule.findIndex(Ms, Ms.matches[k].src_index)]
		dst = Ms.meshes[MeshModule.findIndex(Ms, Ms.matches[k].dst_index)]
		p1 = src.nodes[p];
		p2 = dst.nodes[t[1]] * w[1] + dst.nodes[t[2]] * w[2] + dst.nodes[t[3]] * w[3]
		push!(src_p, p1);
		push!(dst_p, p2);
end

#src_im = imageArray[:, :, Ms.matches_pairs[k][1]];
#dst_im = imageArray[:, :, Ms.matches_pairs[k][2]];

src_im = imageArray[Ms.matches_pairs[k][1]];
dst_im = imageArray[Ms.matches_pairs[k][2]];


square_size = 2*block_size + 1;
cols = round(Int64, sqrt(Ms.matches[k].n) * 2);


sbs_rows = Array{Array{UInt8, 2}, 1}(0);
for r in 1:cols:Ms.matches[k].n
sbs_row = zeros(UInt8, 2*square_size + 2, 0);
	for c in 1:cols
		i = r+c;
		if i > Ms.matches[k].n continue; end;
		iv = round(Int64, src_p[i][1]-Ms.meshes[Ms.matches_pairs[k][1]].disp[1])
		jv = round(Int64, src_p[i][2]-Ms.meshes[Ms.matches_pairs[k][1]].disp[2])
		src_block = src_im[iv-block_size:iv+block_size, jv-block_size:jv+block_size];
		iv = round(Int64, dst_p[i][1]-Ms.meshes[Ms.matches_pairs[k][2]].disp[1])
		jv = round(Int64, dst_p[i][2]-Ms.meshes[Ms.matches_pairs[k][2]].disp[2])
		dst_block = dst_im[iv-block_size:iv+block_size, jv-block_size:jv+block_size];
		sbs = vcat(src_block, zeros(UInt8, 2, square_size), dst_block);
		sbs_row = hcat(sbs_row, zeros(UInt8, 2*square_size+2, 5), sbs);
	end
push!(sbs_rows, sbs_row)
end

row_size = size(sbs_rows[1], 2);

sbs_total = zeros(UInt8, 0, row_size);

for i in 1:size(sbs_rows,1)-1

sbs_total = vcat(sbs_total, sbs_rows[i], zeros(UInt8, 5, row_size));

end

lastrow = sbs_rows[size(sbs_rows, 1)];
lastrow = hcat(lastrow, zeros(UInt8, 2 * square_size+2, row_size-size(lastrow, 2)));

sbs_total = vcat(sbs_total, lastrow);

view(sbs_total);

end

=#


