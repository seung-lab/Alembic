type Mesh
	index							# any sort of index associated with the mesh - by default a 4-tuple

	# all coordinates are taken with the image associated with the mesh having its left top corner at (0, 0) 
	src_nodes::Points					# nodes as array of [i, j] coordinates, sorted in i, j order
	dst_nodes::Points		 			# in the same order as src_nodes, after some transformation

	edges::Edges				                # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
end
	     
### IO EXTENSIONS
function get_path(mesh::Mesh)			return get_path(mesh.index);		end

function get_image(mesh::Mesh, dtype = UInt8)	return get_image(mesh.index, dtype);	end

function get_num_nodes(mesh::Mesh)		return size(mesh.src_nodes[1]);		end

function get_num_edges(mesh::Mesh)		return size(mesh.edges[2]);		end

function get_topleft_offset(mesh::Mesh)		return mesh.src_nodes[1];		end

function get_dims_and_dists(mesh::Mesh)
	n = get_num_nodes(mesh);

	i_dist = 0;
	j_dist = mesh.src_nodes[2][2] - mesh.src_nodes[1][2];
	j_dim = 0;

	# calculate j-dim by iterating until i changes
	for ind in 2:get_num_nodes(mesh)
		i_dif = mesh.src_nodes[ind][1] - mesh.src_nodes[ind-1][1]
		if i_dif == 0 j_dim = ind;
		else i_dist = i_dif; break;
		end
	end

	i_dim = convert(Int64, 1 + (src_nodes[get_num_nodes(mesh)][1] - src_nodes[1][1]) / i_dist);
	return (i_dim, j_dim), [i_dist, j_dist];
end

function Mesh(index)
	params = get_params(index);
	meta = get_metadata(index);
	
	# mesh lengths in each dimension
	dists = [params["mesh_length"] * sin(pi / 3); params["mesh_length"]];

	# dimensions of the mesh as a rectangular grid, maximal in each dimension
	# e.g. a mesh with 5-4-5-4-5-4-5-4 nodes in each row will have dims = (8, 5)
	# 1 is added because of 1-indexing (in length)
	# 2 is added to pad the mesh to extend it beyond by one meshpoint in each direction
	dims = div(get_offsets(meta), dists) + 1 + 2;

	# location of the first node (top left)
	# TODO: Julia does not support rem() for two arrays, so divrem() cannot be used
	topleft_offset = (get_offsets(meta) .% dists) / 2 - dists;

	n = maximum([get_mesh_index(dims, dims[1], dims[2]); get_mesh_index(dims, dims[1], dims[2]-1)]);
	m = 0;
	m_upperbound = 3 * n;

	src_nodes = Points(n);
	edges = spzeros(Float64, n, m_upperbound);

	for i in 1:dims[1], j in 1:dims[2]
		k = get_mesh_index(dims, i, j); if k == 0 continue; end
		src_nodes[k] = get_mesh_coord(dims, topleft_offset, dists, i, j);
		if (j != 1)
			m += 1;	edges[k, m] = -1; edges[get_mesh_index(dims, i, j-1), m] = 1;
		end

		if (i != 1)
			if iseven(i) || j != dims[2]
				m += 1;	edges[k, m] = -1; edges[get_mesh_index(dims, i-1, j), m] = 1;
			end
			if iseven(i) && (j != dims[2]) 			
				m += 1; edges[k, m] = -1; edges[get_mesh_index(dims, i-1, j+1), m] = 1;
			end
			if isodd(i) && (j != 1)
				m += 1; edges[k, m] = -1; edges[get_mesh_index(dims, i-1, j-1), m] = 1;
			end
			if isodd(i) && ((j == 1) || (j == dims[2]))
				m += 1; edges[k, m] = -1; edges[get_mesh_index(dims, i-2, j), m] = 1;
			end
		end
	end

	edges = edges[:, 1:m];
	dst_nodes = copy(src_nodes);

	return Mesh(index, src_nodes, dst_nodes, edges);
end


function get_mesh_index(dims, i, j)
	ind = 0;
	
	if iseven(i) && (j == dims[2]) return ind; end
	if ((i < 1) || (j < 1) || (i > dims[1]) || (j > dims[2])) return ind; end
	
	ind += div(i-1, 2) * (dims[2] - 1); #even rows
	ind += div(i, 2) * dims[2]; #odd rows
	ind += j;
	ind = convert(Int64, ind);
	return ind;
end

function get_mesh_coord(dims, total_offset, dists, i, j)
	if iseven(i) && (j == dims[2]) return (0, 0); end
	
	pi = (i-1) * dists[1] + total_offset[1];

	if iseven(i)	pj = (j-0.5) * dists[2] + total_offset[2];
	else		pj = (j-1) * dists[2] + total_offset[2];
	end

	return [pi; pj];
end



# find the triangular mesh indices for a given point in mesh image coordinates
function find_mesh_triangle(mesh::Mesh, point::Point)
	
	dims, dists = get_dims_and_dists(mesh); 
	point_padded = point - get_topleft_offset(mesh);

	indices_raw = point_padded ./ dists

	# find which rows the point belongs to
	i0 = round(Int64, indices_raw[1] + 1);
	
	if isodd(i0)
		j0 = round(Int64, indices_raw[2] + 1);
	else
		j0 = round(Int64, indices_raw[2] + 0.5);
	end

	ind0 = get_mesh_index(dims, i0, j0);
	
	if ind0 == 0
		return NO_TRIANGLE;
	end

	node0 = mesh.nodes[ind0];

	res = point - node0;

	theta = abs(atan(res[1] / res[2]));
	if (theta < pi / 3)
		if (dj >= 0)
			ind1 = get_mesh_index(dims, i0, j0 + 1);
			if (di >= 0)
				if isodd(i0)
					ind2 = get_mesh_index(dims, i0+1, j0);
				else
					ind2 = get_mesh_index(dims, i0+1, j0+1);
				end
			else
				if isodd(i0)
					ind2 = get_mesh_index(dims, i0-1, j0);
				else
					ind2 = get_mesh_index(dims, i0-1, j0+1);
				end
			end
		else
			ind1 = get_mesh_index(dims, i0, j0 - 1);
			if (di >= 0)
				if isodd(i0)
					ind2 = get_mesh_index(dims, i0+1, j0-1);
				else
					ind2 = get_mesh_index(dims, i0+1, j0);
				end
			else
				if isodd(i0)
					ind2 = get_mesh_index(dims, i0-1, j0-1);
				else
					ind2 = get_mesh_index(dims, i0-1, j0);
				end
			end
		end
	else
		if (di >= 0)
			if isodd(i0)
				ind1 = get_mesh_index(dims, i0+1, j0-1);
				ind2 = get_mesh_index(dims, i0+1, j0);
			else
				ind1 = get_mesh_index(dims, i0+1, j0);
				ind2 = get_mesh_index(dims, i0+1, j0+1);
			end
		else
			if isodd(i0)
				ind1 = get_mesh_index(dims, i0-1, j0-1);
				ind2 = get_mesh_index(dims, i0-1, j0);
			else
				ind1 = get_mesh_index(dims, i0-1, j0);
				ind2 = get_mesh_index(dims, i0-1, j0+1);
			end
		end
	end

	if (ind1 == 0 || ind2 == 0)
		return NO_TRIANGLE;
	end
	return (ind0, ind1, ind2);
end

# Convert Cartesian coordinate to triple of barycentric coefficients
function get_triangle_weights(m::Mesh, point::Point, triangle::Triangle)
	R = vcat(mesh.nodes[triangle[1]]', mesh.nodes[triangle[2]]', mesh.nodes[triangle[3]]')
	R = hcat(R, ones(Float64, 3, 1));
	r = vcat(point, 1.0);

	V = r * R^-1;

	return (V[1], V[2], V[3]);
end

