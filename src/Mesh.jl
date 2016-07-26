type Mesh
	index::Index						# any sort of index associated with the mesh - by default a 4-tuple

	# all coordinates are taken with the image associated with the mesh having its left top corner at (0, 0) 
	src_nodes::Points					# nodes as array of [i, j] coordinates, sorted in i, j order
	dst_nodes::Points		 			# in the same order as src_nodes, after some transformation

	edges::Edges				                # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
	properties::Dict{Any, Any}
end

### INDEX.jl EXTENSIONS
function is_adjacent(Am::Mesh, Bm::Mesh)		return is_adjacent(Am.index, Bm.index);			end
function is_diagonal(Am::Mesh, Bm::Mesh)		return is_diagonal(Am.index, Bm.index);			end
function is_preceding(Am::Mesh, Bm::Mesh, within = 1)	return is_preceding(Am.index, Bm.index, within);	end

### PARAMS.jl EXTENSIONS
function get_params(mesh::Mesh)				return get_params(mesh.index);				end
function globalize!(pts::Points, mesh::Mesh)
  offset = get_offset(mesh)
  @simd for i in 1:length(pts) @fastmath @inbounds pts[i] = pts[i] + offset; end
end
	     
### META.jl EXTENSIONS
function get_offset(mesh::Mesh)				return get_offset(mesh.index);				end
function get_image_size(mesh::Mesh)			return get_image_size(get_metadata(mesh.index));	end
function get_metadata(mesh::Mesh)			return get_metadata(mesh.index);			end

### IO.jl EXTENSIONS
function get_image_path(mesh::Mesh)			return get_image_path(mesh.index);			end
function get_image(mesh::Mesh; kwargs...)		return get_image(mesh.index; kwargs...);		end
#function get_name(mesh::Mesh)				return get_name(mesh.index);				end

### retrieval
function get_index(mesh::Mesh)				return mesh.index;					end
function get_nodes(mesh::Mesh; globalized::Bool = false, use_post::Bool=false)
	nodes = use_post ? copy(mesh.dst_nodes) : copy(mesh.src_nodes);
	globalized ? globalize!(nodes, mesh) : nothing
	return nodes
end

### IO

### counting
function count_nodes(mesh::Mesh)			return size(mesh.src_nodes, 1);				end
function count_edges(mesh::Mesh)			return size(mesh.edges, 2);				end

### internal
function get_topleft_offset(mesh::Mesh)			return mesh.src_nodes[1];				end

function get_dims_and_dists(mesh::Mesh)
	n = count_nodes(mesh);
	@fastmath @inbounds begin
	i_dist = 0;
	j_dist = mesh.src_nodes[2][2] - mesh.src_nodes[1][2];
	j_dim = 0;

	# calculate j-dim by iterating until i changes
	for ind in 2:count_nodes(mesh)
		i_dif = mesh.src_nodes[ind][1] - mesh.src_nodes[ind-1][1]
		if i_dif == 0 j_dim = ind; else i_dist = i_dif; break; end
	end
	i_dim = round(Int64, 1 + (mesh.src_nodes[count_nodes(mesh)][1] - mesh.src_nodes[1][1]) / i_dist);
      end #fmib
	return (i_dim, j_dim), [i_dist, j_dist];
end

function get_mesh_index(dims, i, j)
	ind = 0;
	if iseven(i) && (j == dims[2]) return ind; end
	if ((i < 1) || (j < 1) || (i > dims[1]) || (j > dims[2])) return ind; end
	@fastmath @inbounds ind += div(i-1, 2) * (dims[2] - 1); #even rows
	@fastmath @inbounds ind += div(i, 2) * dims[2]; #odd rows
	ind += j;
	ind = convert(Int64, ind);
	return ind;
end

function get_mesh_coord(dims, total_offset, dists, i, j)
	if iseven(i) && (j == dims[2]) return (0, 0); end
	@fastmath @inbounds pi = (i-1) * dists[1] + total_offset[1];
	if iseven(i)	@fastmath @inbounds pj = (j-0.5) * dists[2] + total_offset[2];
	else		@fastmath @inbounds pj = (j-1) * dists[2] + total_offset[2];
	end
	return [pi; pj];
end

### edge-length related

# returns the endpoints of the edges for all edges - returns two Points arrays, where each edge runs from endpoint_a[i] to endpoint_b[i]
function get_edge_endpoints(mesh::Mesh; globalized::Bool = false, use_post::Bool = false)
	num_edges = count_edges(mesh);
	endpoints_a = Points(num_edges);
	endpoints_b = Points(num_edges);

  	for ind in 1:num_edges
	@inbounds node_inds = findnz(mesh.edges[:, ind])[1];
	@inbounds endpoints_a[ind] = use_post ? mesh.src_nodes[node_inds[1]] : mesh.dst_nodes[node_inds[1]]
	@inbounds endpoints_b[ind] = use_post ? mesh.src_nodes[node_inds[2]] : mesh.dst_nodes[node_inds[2]]
      	end

	globalized ? globalize!(endpoints_a, mesh) : nothing
	globalized ? globalize!(endpoints_b, mesh) : nothing
	return endpoints_a, endpoints_b
end

# returns the midpoints of the edges for all edges - kwargs are parsed to get_edge_endpoints
function get_edge_midpoints(mesh::Mesh; kwargs...)
  	endpoints_a, endpoints_b = get_edge_endpoints(mesh; kwargs...);
	@simd for i in 1:length(endpoints_a)
		@fastmath @inbounds endpoints_a[i] = endpoints_a[i] + endpoints_b[i]
		@fastmath @inbounds endpoints_a[i] = endpoints_a[i] / 2
	end
      return endpoints_a;
end

# returns the lengths of the edges for all edges - kwargs are parsed to get_edge_endpoints, though globalized should not matter at all
function get_edge_lengths(mesh::Mesh; kwargs...)
  	endpoints_a, endpoints_b = get_edge_endpoints(mesh; kwargs...);
        edgelengths = similar(endpoints_a, Float64)
	@simd for i in 1:length(endpoints_a)
		@fastmath @inbounds endpoints_a[i] = endpoints_a[i] - endpoints_b[i]
	        @fastmath @inbounds edgelengths[i] = norm(endpoints_a[i])
	end
      return edgelengths
end

function remove_edge!(mesh::Mesh, ind)
	if !haskey(mesh.properties, "removed_indices")
		mesh.properties["removed_indices"] = Set{Int64}()
	end
	push!(mesh.properties["removed_indices"], ind)
	# edges_to_keep = get_edge_indices(mesh)
	# deleteat!(edges_to_keep, ind)
	# mesh.edges = mesh.edges[:, [edges_to_keep]]
end

function get_removed_edge_indices(mesh::Mesh)
	if haskey(mesh.properties, "removed_indices")
		return collect(mesh.properties["removed_indices"])
	end
	return []
end

function get_fixed_edge_indices(mesh::Mesh)
	if haskey(mesh.properties, "fixed_indices")
		return collect(mesh.properties["fixed_indices"])
	end
	return []
end

function get_edge_indices(mesh::Mesh)
	return collect(1:count_edges(mesh))
end

# Not sure why this doesn't work
function isless(meshA::Mesh, meshB::Mesh)
  return get_index(meshA) < get_index(meshB)
end

function isequal(meshA::Mesh, meshB::Mesh)
  return get_index(meshA) == get_index(meshB)
end


### INIT
function make_mesh(index, params = get_params(index), fixed=false)
	println("Creating mesh for $index")
	# mesh lengths in each dimension
	dists = [params["mesh"]["mesh_length"] * sin(pi / 3); params["mesh"]["mesh_length"]];

	# dimensions of the mesh as a rectangular grid, maximal in each dimension
	# e.g. a mesh with 5-4-5-4-5-4-5-4 nodes in each row will have dims = (8, 5)
	# 1 is added because of 1-indexing (in length)
	# 2 is added to pad the mesh to extend it beyond by one meshpoint in each direction
	dims = round(Int64, div(get_image_size(index), dists)) + 1 + 2;

	# location of the first node (top left)
	# TODO: Julia does not support rem() for two arrays, so divrem() cannot be used
	topleft_offset = (get_image_size(index) .% dists) / 2 - dists;

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

function get_triangle_weights!(mesh::Mesh, point::Point, triangle::Triangle, R, r)

	if triangle == NO_TRIANGLE::Triangle return NO_WEIGHTS::Weight; end

	@inbounds R[1:2, 1] = mesh.src_nodes[triangle[1]];
	@inbounds R[1:2, 2] = mesh.src_nodes[triangle[2]];
	@inbounds R[1:2, 3] = mesh.src_nodes[triangle[3]];
	
	@inbounds r[1:2] = point;

	@fastmath V = R \ r;

	return (V[1], V[2], V[3]);
end

function get_triangle_weights(mesh::Mesh, points::Points, triangles::Triangles)
  	R = ones(Float64, 3, 3)
  	r = ones(Float64, 3)

	return Weights(map(get_triangle_weights!, repeated(mesh), points, triangles, repeated(R), repeated(r)))
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

function deform(points::Points, mesh::Mesh)
  	triangles = find_mesh_triangle(mesh, points);
	weights = get_triangle_weights(mesh, points, triangles);
	return get_tripoint_dst(mesh, triangles, weights);
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

function reset!(mesh::Mesh)
	mesh.dst_nodes = copy(mesh.src_nodes)
end

function get_edges(mesh::Mesh, node_index)
	return findnz(mesh.edges[node_index, :])
end

function get_edges(mesh::Mesh, triangle::Triangle)
	edges = Set()
	for node in triangle
		edges = intersect(edges, get_edges(mesh, node))
	end
	return edges
end
