#using Julimaps
#using Images
#importall IO

type Mesh
	name::String;						# name of the image file

	index::Index						# wafer, section, tile index of the mesh				tileIndex = 0 if the tile is a whole section
	disp::Point						# displacement of the tile within the section, y, x			(0, 0) if the tile starts at the top left corner

	dims::Pairing				 		# mesh dimensions in terms of nodes in the i direction, j direction. 	(0, 0) if the mesh is not a regular mesh
	offsets::Point				                # mesh offset in terms of the top left, from the image. 		[0; 0] if the mesh is not a regular mesh
	dists::Point				                # mesh distance in each direction. 					[0; 0] if the mesh is not a regular mesh

	n::Int64 						# number of nodes in the mesh
	m::Int64	 					# number of edges in the mesh

	nodes::Points						# 2-by-n dense matrix of nodes, each column stores [i-coordinate, j-coordinate] in global coordinates
	nodes_t::Points			 			# 2-by-n dense matrix of nodes after transformation, in the same order as nodes - by default same as nodes
	nodes_fixed::BinaryProperty				# 1-by-n dense vector of nodes that denote whether points are fixed

	edges::Edges				                # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
	edge_lengths::FloatProperty	  			# 1-by-m dense vector of floats that stores the rest lengths.
	edge_coeffs::FloatProperty	   			# 1-by-m dense vector of floats that stores the spring coefficients.
end

### IO EXTENSIONS
function getPath(mesh::Mesh)
	return getPath(mesh.name);
end

function getFloatImage(mesh::Mesh)
	return getFloatImage(getPath(mesh.index));
end

function getImage(mesh::Mesh)
	return getImage(mesh.name);
end

function Tile2Mesh(name, size::Int64, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
	return Tile2Mesh(name, size, size, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
end

function Tile2Mesh(name, size_i::Int64, size_j::Int64, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
	(Ai, Aj) = (size_i,size_j);

	dists = [mesh_length * sin(pi / 3); mesh_length];
	dims = (convert(Int64, div(Ai, dists[1]) + 1), convert(Int64, div(Aj, dists[2]) + 1));
 	offsets = [rem(Ai, dists[1])/2; rem(Aj, dists[2])/2];
	disp = [dy; dx];

	n = maximum([getMeshIndex(dims, dims[1], dims[2]); getMeshIndex(dims, dims[1], dims[2]-1)]);
	m = 0;
	m_upperbound = 3 * n;

	nodes = Points(n);
	nodes_fixed = BinaryProperty(n); nodes_fixed[:] = tile_fixed;
	edges = spzeros(Float64, n, m_upperbound);
	edge_lengths = FloatProperty(m_upperbound); edge_lengths[:] = convert(Float64, mesh_length);
	edge_coeffs = FloatProperty(m_upperbound); edge_coeffs[:] = convert(Float64, mesh_coeff);

	for i in 1:dims[1], j in 1:dims[2]
		k = getMeshIndex(dims, i, j); if k == 0 continue; end
		nodes[k] = getMeshCoord(dims, disp+offsets, dists, i, j);
		if (j != 1)
			m += 1;	edges[k, m] = -1; edges[getMeshIndex(dims, i, j-1), m] = 1;
		end

		if (i != 1)
			if iseven(i) || j != dims[2]
				m += 1;	edges[k, m] = -1; edges[getMeshIndex(dims, i-1, j), m] = 1;
			end
			if iseven(i) && (j != dims[2]) 			
				m += 1; edges[k, m] = -1; edges[getMeshIndex(dims, i-1, j+1), m] = 1;
			end
			if isodd(i) && (j != 1)
				m += 1; edges[k, m] = -1; edges[getMeshIndex(dims, i-1, j-1), m] = 1;
			end
			if isodd(i) && ((j == 1) || (j == dims[2]))
				m += 1; edges[k, m] = -1; edges[getMeshIndex(dims, i-2, j), m] = 1;
				edge_lengths[m] = 2 * dists[1];
			end
		end
	end

	edges = edges[:, 1:m];
	edge_lengths = edge_lengths[1:m];
	edge_coeffs = edge_coeffs[1:m];

	return Mesh(name, index, disp, dims, offsets, dists, n, m, nodes, nodes, nodes_fixed, edges, edge_lengths, edge_coeffs);
end


# Tile2Mesh
function Tile2Mesh(name, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
	image = getImage(getPath(name));
	return Tile2Mesh(name, image, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
end

function Tile2Mesh(name, image, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
	(size_i, size_j) = size(image);
	Tile2Mesh(name, size_i::Int64, size_j::Int64, index, dy, dx, tile_fixed, mesh_length, mesh_coeff)
end

function getMeshIndex(dims, i, j)
	ind = 0;
	
	if iseven(i) && (j == dims[2]) return ind; end
	if ((i < 1) || (j < 1) || (i > dims[1]) || (j > dims[2])) return ind; end
	
	ind += div(i-1, 2) * (dims[2] - 1); #even rows
	ind += div(i, 2) * dims[2]; #odd rows
	ind += j;
	ind = convert(Int64, ind);
	return ind;
end

function getMeshCoord(dims, total_offset, dists, i, j)
	if iseven(i) && (j == dims[2]) return (0, 0); end
	
	pi = (i-1) * dists[1] + total_offset[1];

	if iseven(i)	pj = (j-0.5) * dists[2] + total_offset[2];
	else		pj = (j-1) * dists[2] + total_offset[2];
	end

	return [pi; pj];
end

function isInternal(A, i, j, d)
	if i > d && i <= (size(A, 1) - d) && j > d && j <= (size(A, 2) - d)
		return true
	end
	return false
end

# find the triangular mesh indices for a given point in A
function findMeshTriangle(Am, i, j)

	# convert to A's local coordinates, and displace by the mesh offset
	Ai = i - Am.disp[1] - Am.offsets[1];
	Aj = j - Am.disp[2] - Am.offsets[2];

	# find which rows the point belongs to
	i0 = round(Int64, Ai / Am.dists[1] + 1);

	if isodd(i0)
		j0 = round(Int64, Aj / Am.dists[2] + 1);
	else
		j0 = round(Int64, Aj / Am.dists[2] + 0.5);
	end

	ind0 = getMeshIndex(Am.dims, i0, j0);
	
	if ind0 == 0
		return noTriangle;
	end

	node0 = Am.nodes[ind0];

	#println("$i0, $j0 at ($node0)");

	di = i - node0[1];
	dj = j - node0[2];

	#println("$di, $dj");

	theta = abs(atan(di / dj));
	if (theta < pi / 3)
		if (dj >= 0)
			ind1 = getMeshIndex(Am.dims, i0, j0 + 1);
			if (di >= 0)
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0+1, j0);
				else
					ind2 = getMeshIndex(Am.dims, i0+1, j0+1);
				end
			else
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0-1, j0);
				else
					ind2 = getMeshIndex(Am.dims, i0-1, j0+1);
				end
			end
		else
			ind1 = getMeshIndex(Am.dims, i0, j0 - 1);
			if (di >= 0)
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0+1, j0-1);
				else
					ind2 = getMeshIndex(Am.dims, i0+1, j0);
				end
			else
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0-1, j0-1);
				else
					ind2 = getMeshIndex(Am.dims, i0-1, j0);
				end
			end
		end
	else
		if (di >= 0)
			if isodd(i0)
				ind1 = getMeshIndex(Am.dims, i0+1, j0-1);
				ind2 = getMeshIndex(Am.dims, i0+1, j0);
			else
				ind1 = getMeshIndex(Am.dims, i0+1, j0);
				ind2 = getMeshIndex(Am.dims, i0+1, j0+1);
			end
		else
			if isodd(i0)
				ind1 = getMeshIndex(Am.dims, i0-1, j0-1);
				ind2 = getMeshIndex(Am.dims, i0-1, j0);
			else
				ind1 = getMeshIndex(Am.dims, i0-1, j0);
				ind2 = getMeshIndex(Am.dims, i0-1, j0+1);
			end
		end
	end

	if (ind1 == 0 || ind2 == 0)
		return noTriangle;
	end
	return (ind0, ind1, ind2);
end

# Convert Cartesian coordinate to triple of barycentric coefficients
function getTriangleWeights(Am, triangle, pi, pj)
	R = vcat(Am.nodes[triangle[1]]', Am.nodes[triangle[2]]', Am.nodes[triangle[3]]')
	R = hcat(R, ones(Float64, 3, 1));
	r = hcat(pi, pj, 1.0);

	V = r * R^-1;

	return (V[1], V[2], V[3]);
end
