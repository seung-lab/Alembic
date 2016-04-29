type Mesh
	index::Index							# any sort of index associated with the mesh - by default a 4-tuple

	# all coordinates are taken with the image associated with the mesh having its left top corner at (0, 0) 
	src_nodes::Points					# nodes as array of [i, j] coordinates, sorted in i, j order
	dst_nodes::Points		 			# in the same order as src_nodes, after some transformation

	edges::Edges				                # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
	properties::Dict{Any, Any}
end

### INDEX.jl EXTENSIONS
function is_adjacent(Am::Mesh, Bm::Mesh)	return is_adjacent(Am.index, Bm.index);		end
function is_diagonal(Am::Mesh, Bm::Mesh)	return is_diagonal(Am.index, Bm.index);		end
function is_preceding(Am::Mesh, Bm::Mesh, within = 1)	return is_preceding(Am.index, Bm.index, within);	end

### PARAMS.jl EXTENSIONS
function get_params(mesh::Mesh)			return get_params(mesh.index);		end
	     
### META.jl EXTENSIONS
function get_offset(mesh::Mesh)		return get_offset(mesh.index);		end
function get_image_sizes(mesh::Mesh)		return get_image_sizes(get_metadata(mesh.index));	end
function get_metadata(mesh::Mesh)		return get_metadata(mesh.index);	end

### IO.jl EXTENSIONS
function get_path(mesh::Mesh)			return get_path(mesh.index);		end
function get_image(mesh::Mesh, scale=1.0, dtype = IMG_ELTYPE)	return get_image(mesh.index, scale, dtype);	end

### counting
function count_nodes(mesh::Mesh)		return size(mesh.src_nodes, 1);		end
function count_edges(mesh::Mesh)		return size(mesh.edges, 2);		end

function get_index(mesh::Mesh)		return mesh.index;		end

### internal
function get_topleft_offset(mesh::Mesh)		return mesh.src_nodes[1];		end
function get_edge_points(mesh::Mesh, ind)
#	src_ind = findfirst(this -> this < 0, mesh.edges[:, ind]);
#	dst_ind = findfirst(this -> this > 0, mesh.edges[:, ind]);
	@fastmath @inbounds inds = findnz(mesh.edges[:, ind])[1];
	return mesh.src_nodes[inds[1]], mesh.src_nodes[inds[2]]
end

# Not sure why this doesn't work
function isless(meshA::Mesh, meshB::Mesh)
  return get_index(meshA) < get_index(meshB)
end

function isequal(meshA::Mesh, meshB::Mesh)
  return get_index(meshA) == get_index(meshB)
end

function get_edge_length(mesh::Mesh, ind)	start_pt, end_pt = get_edge_points(mesh, ind);
      @fastmath return norm(start_pt - end_pt);
end

function get_edge_lengths(mesh::Mesh)		return map(get_edge_length, repeated(mesh), collect(1:count_edges(mesh)));		end

function get_homogenous_edge_lengths(mesh::Mesh)	return fill(get_edge_length(mesh, 1), count_edges(mesh));		end

function get_globalized_nodes_h(mesh::Mesh)
	@fastmath g_src_nodes, g_dst_nodes = get_globalized_nodes(mesh);
    return hcat(g_src_nodes...), hcat(g_dst_nodes...)
end
function get_globalized_nodes(mesh::Mesh)
    @fastmath g_src_nodes = mesh.src_nodes + fill(get_offset(mesh), count_nodes(mesh));
    @fastmath g_dst_nodes = mesh.dst_nodes + fill(get_offset(mesh), count_nodes(mesh));
    return g_src_nodes, g_dst_nodes
end

function get_dims_and_dists(mesh::Mesh)
	n = count_nodes(mesh);
	@fastmath @inbounds begin
	i_dist = 0;
	j_dist = mesh.src_nodes[2][2] - mesh.src_nodes[1][2];
	j_dim = 0;

	# calculate j-dim by iterating until i changes
	for ind in 2:count_nodes(mesh)
		i_dif = mesh.src_nodes[ind][1] - mesh.src_nodes[ind-1][1]
		if i_dif == 0 j_dim = ind;
		else i_dist = i_dif; break;
		end
	end
	i_dim = round(Int64, 1 + (mesh.src_nodes[count_nodes(mesh)][1] - mesh.src_nodes[1][1]) / i_dist);
      end #fmib
	return (i_dim, j_dim), [i_dist, j_dist];
end

### INIT
function Mesh(index, params = get_params(index), fixed=false)
	
	# mesh lengths in each dimension
	dists = [params["mesh"]["mesh_length"] * sin(pi / 3); params["mesh"]["mesh_length"]];

	# dimensions of the mesh as a rectangular grid, maximal in each dimension
	# e.g. a mesh with 5-4-5-4-5-4-5-4 nodes in each row will have dims = (8, 5)
	# 1 is added because of 1-indexing (in length)
	# 2 is added to pad the mesh to extend it beyond by one meshpoint in each direction
	dims = round(Int64, div(get_image_sizes(index), dists)) + 1 + 2;

	# location of the first node (top left)
	# TODO: Julia does not support rem() for two arrays, so divrem() cannot be used
	topleft_offset = (get_image_sizes(index) .% dists) / 2 - dists;

	n = maximum([get_mesh_index(dims, dims[1], dims[2]); get_mesh_index(dims, dims[1], dims[2]-1)]);
	m = 0;
	m_upperbound = 3 * n;

	src_nodes = Points(n);
	edges = spzeros(Float64, n, m_upperbound);

	@fastmath @inbounds for i in 1:dims[1], j in 1:dims[2]
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

	@inbounds edges = edges[:, 1:m];
	dst_nodes = copy(src_nodes);

	properties = Dict{Any, Any}(
				    "params" => params,
				    "fixed" => fixed);

	return Mesh(index, src_nodes, dst_nodes, edges, properties);
end


function get_mesh_index(dims, i, j)
  @fastmath @inbounds begin
	ind = 0;
	
	if iseven(i) && (j == dims[2]) return ind; end
	if ((i < 1) || (j < 1) || (i > dims[1]) || (j > dims[2])) return ind; end
	
	ind += div(i-1, 2) * (dims[2] - 1); #even rows
	ind += div(i, 2) * dims[2]; #odd rows
	ind += j;
	ind = convert(Int64, ind);
      end#fmib
	return ind;
end

function get_mesh_coord(dims, total_offset, dists, i, j)
  @fastmath @inbounds begin
	if iseven(i) && (j == dims[2]) return (0, 0); end
	
	pi = (i-1) * dists[1] + total_offset[1];

	if iseven(i)	pj = (j-0.5) * dists[2] + total_offset[2];
	else		pj = (j-1) * dists[2] + total_offset[2];
	end

      end#fmib

	return [pi; pj];
end



# find the triangular mesh indices for a given point in mesh image coordinates
function find_mesh_triangle(mesh::Mesh, point::Point)
  @inbounds begin	
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

	node0 = mesh.src_nodes[ind0];

	res = point - node0;


	theta = abs(atan(res[1] / res[2]));
	if (theta < pi / 3)
		if (res[2] >= 0)
			if (res[1] >= 0)
				if isodd(i0)
					ind1 = get_mesh_index(dims, i0, j0 + 1);
					ind2 = get_mesh_index(dims, i0+1, j0);
				else
					ind1 = get_mesh_index(dims, i0, j0 + 1);
					ind2 = get_mesh_index(dims, i0+1, j0+1);
					if ind1 == 0
						ind1 = get_mesh_index(dims, i0-1, j0 + 1);
					end
				end
			else
				if isodd(i0)
					ind1 = get_mesh_index(dims, i0, j0 + 1);
					ind2 = get_mesh_index(dims, i0-1, j0);
				else
					ind1 = get_mesh_index(dims, i0, j0 + 1);
					ind2 = get_mesh_index(dims, i0-1, j0+1);
					if ind1 == 0
						ind1 = get_mesh_index(dims, i0+1, j0 + 1);
					end
				end
			end
		else
			if (res[1] >= 0)
				if isodd(i0)
					ind1 = get_mesh_index(dims, i0, j0 - 1);
					ind2 = get_mesh_index(dims, i0+1, j0-1);
				else
					ind1 = get_mesh_index(dims, i0, j0 - 1);
					ind2 = get_mesh_index(dims, i0+1, j0);
					if ind1 == 0
						ind1 = get_mesh_index(dims, i0-1, j0);
					end
				end
			else
				if isodd(i0)
					ind1 = get_mesh_index(dims, i0, j0 - 1);
					ind2 = get_mesh_index(dims, i0-1, j0-1);
				else
					ind1 = get_mesh_index(dims, i0, j0 - 1);
					ind2 = get_mesh_index(dims, i0-1, j0);
					if ind1 == 0
						ind1 = get_mesh_index(dims, i0+1, j0);
					end
				end
			end
		end
	else
		if (res[1] >= 0)
			if isodd(i0)
				ind1 = get_mesh_index(dims, i0+1, j0-1);
				ind2 = get_mesh_index(dims, i0+1, j0);
				if ind1 == 0
					ind1 = get_mesh_index(dims, i0+2, j0);
				end
				if ind2 == 0
					ind2 = get_mesh_index(dims, i0+2, j0);
				end
			else
				ind1 = get_mesh_index(dims, i0+1, j0);
				ind2 = get_mesh_index(dims, i0+1, j0+1);
			end
		else
			if isodd(i0)
				ind1 = get_mesh_index(dims, i0-1, j0-1);
				ind2 = get_mesh_index(dims, i0-1, j0);
				if ind1 == 0
					ind1 = get_mesh_index(dims, i0-2, j0);
				end
				if ind2 == 0
					ind2 = get_mesh_index(dims, i0-2, j0);
				end
			else
				ind1 = get_mesh_index(dims, i0-1, j0);
				ind2 = get_mesh_index(dims, i0-1, j0+1);
			end
		end
	end

	if (ind1 == 0 || ind2 == 0)
	  	#println("Missing triangle");
		#println("$res at $ind0");
		return NO_TRIANGLE;
	end

      end #ib
	return (ind0, ind1, ind2);
end

function find_mesh_triangle(mesh::Mesh, points::Points)
	return Triangles(map(find_mesh_triangle, repeated(mesh), points));
end

# Convert Cartesian coordinate to triple of barycentric coefficients
function get_triangle_weights(mesh::Mesh, point::Point, triangle::Triangle)
	if triangle == NO_TRIANGLE::Triangle return NO_WEIGHTS::Weight; end
	@inbounds R = vcat(mesh.src_nodes[triangle[1]]', mesh.src_nodes[triangle[2]]', mesh.src_nodes[triangle[3]]')
	R = hcat(R, ones(Float64, 3, 1));
	r = vcat(point, 1.0);

	V = r' * R^-1;

	return (V[1], V[2], V[3]);
end

function get_triangle_weights(mesh::Mesh, points::Points, triangles::Triangles)

	return Weights(map(get_triangle_weights, repeated(mesh), points, triangles))

end

function get_tripoint_dst(mesh::Mesh, triangles::Triangles, weights::Weights)
	return Points(map(get_tripoint_dst, repeated(mesh), triangles, weights))
end

function get_tripoint_dst(mesh::Mesh, triangle::Triangle, weight::Weight)
  @fastmath @inbounds begin
	if triangle == NO_TRIANGLE return NO_POINT; end
	dst_trinodes = mesh.dst_nodes[triangle[1]], mesh.dst_nodes[triangle[2]], mesh.dst_nodes[triangle[3]]
	dst_point = dst_trinodes[1] * weight[1] + dst_trinodes[2] * weight[2] + dst_trinodes[3] * weight[3];
      end #fmib
	return dst_point;
end

function fix!(mesh::Mesh)
	mesh.properties["fixed"] = true;
end

function unfix!(mesh::Mesh)
	mesh.properties["fixed"] = false;
end

function is_fixed(mesh::Mesh)
	get(mesh.properties, "fixed", false)
end
