struct Mesh{T} <: AbstractMesh
	index							# any sort of index associated with the mesh - by default a 4-tuple

	# all coordinates are taken with the image associated with the mesh having its left top corner at (0, 0) 
	src_nodes::Points{T}					# nodes as array of [i, j] coordinates, sorted in i, j order
	dst_nodes::Points{T}		 			# in the same order as src_nodes, after some transformation

	# n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
	# represented by three arrays - use sparse(edges_I, edges_J, edges_V)
	edges_I::Array{Int64, 1}
	edges_J::Array{Int64, 1}
	edges_V::Array{T, 1}

	properties::Dict{Symbol, Any}				# properties field
end

function Base.deepcopy(m::Mesh; index=m.index, src_nodes=m.src_nodes, 
				dst_nodes=m.dst_nodes, edges_I=m.edges_I, edges_J=m.edges_J, edges_V=m.edges_V, properties=m.properties)
	return Mesh(index, deepcopy(src_nodes), deepcopy(dst_nodes), deepcopy(edges_I), deepcopy(edges_J), deepcopy(edges_V), deepcopy(properties))
end


### INDEX.jl EXTENSIONS
function is_adjacent(Am::Mesh, Bm::Mesh)		
	return is_adjacent(Am.index, Bm.index);			
end

function is_diagonal(Am::Mesh, Bm::Mesh)		
	return is_diagonal(Am.index, Bm.index);			
end

function is_preceding(Am::Mesh, Bm::Mesh, within = 1)	
	return is_preceding(Am.index, Bm.index, within);	
end

function globalize!{T}(pts::Points{T}, mesh::Mesh{T})
  offset = get_offset(mesh)
  globalize!(pts, offset)
end
    
### META.jl EXTENSIONS
function get_offset(mesh::Mesh)				
	return get_offset(:match_image, mip=get_mip(:match_image));				
end

function get_image_size(mesh::Mesh)			
	return get_image_size(:match_image, mip=get_mip(:match_image));			
end

function get_metadata(mesh::Mesh)			
	return get_metadata(mesh.index);										
end

### IO.jl EXTENSIONS
function get_image(mesh::Mesh, obj_name=:match_image; mip=get_mip(:match_image))
	return get_image(get_index(mesh), obj_name, mip=mip);	
end

### retrieval
function get_index(mesh::Mesh)				
	return mesh.index;	
end
function is_subsection(mesh::Mesh)				
	return is_subsection(mesh.index);	
end

function get_z(mesh::Mesh)
	return get_z(mesh.index)
end

function get_nodes(mesh::Mesh; globalized::Bool = false, use_post::Bool=false, scale=1.0)
	nodes = use_post ? deepcopy(mesh.dst_nodes) : deepcopy(mesh.src_nodes);
	globalized ? globalize!(nodes, mesh) : nothing
	if scale != 1.0
		nodes *= scale
	end
	return nodes
end

function get_bbox(mesh::Mesh; globalized::Bool = false, use_post::Bool=false, scale=1.0)
	nodes = get_nodes(mesh, globalized=globalized, use_post=use_post, scale=scale)
	return ImageRegistration.find_mesh_bb(nodes)
end

### counting
function count_nodes(mesh::Mesh)			return size(mesh.src_nodes, 2);				end
function count_edges(mesh::Mesh)			return div(length(mesh.edges_V), 2);				end

### internal
function get_topleft_offset(mesh::Mesh)			return mesh.src_nodes[:, 1];				end

function get_dims_and_dists(mesh::Mesh)
	n = count_nodes(mesh);
	@fastmath @inbounds begin
	i_dist = 0;
	j_dist = mesh.src_nodes[2,2] - mesh.src_nodes[2, 1];
	j_dim = 0;

	# calculate j-dim by iterating until i changes
	for ind in 2:count_nodes(mesh)
		i_dif = mesh.src_nodes[1, ind] - mesh.src_nodes[1,ind-1]
		if i_dif == 0 j_dim = ind; else i_dist = i_dif; break; end
	end
	i_dim = round(Int64, 1 + (mesh.src_nodes[1,count_nodes(mesh)] - mesh.src_nodes[1, 1]) / i_dist);
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
function get_edge_endpoints{T}(mesh::Mesh{T}; globalized::Bool = false, use_post::Bool = false)
	num_edges = count_edges(mesh);
	endpoints_a = Points{T}(2, num_edges);
	endpoints_b = Points{T}(2, num_edges);

  	for ind in 1:num_edges
	@fastmath node_loc = ind * 2 - 1
	@inbounds node_inds = view(mesh.edges_I, node_loc : node_loc+1);
	@inbounds endpoints_a[:, ind] = use_post ? mesh.dst_nodes[:, node_inds[1]] : mesh.src_nodes[:, node_inds[1]]
	@inbounds endpoints_b[:, ind] = use_post ? mesh.dst_nodes[:, node_inds[2]] : mesh.src_nodes[:, node_inds[2]]      	
        end

	globalized ? globalize!(endpoints_a, mesh) : nothing
	globalized ? globalize!(endpoints_b, mesh) : nothing
	return endpoints_a, endpoints_b
end

# returns the midpoints of the edges for all edges - kwargs are parsed to get_edge_endpoints
function get_edge_midpoints(mesh::Mesh; kwargs...)
  	endpoints_a, endpoints_b = get_edge_endpoints(mesh; kwargs...);
	@simd for i in 1:size(endpoints_a, 2)
		@fastmath @inbounds endpoints_a[1, i] = endpoints_a[1, i] + endpoints_b[1, i]
		@fastmath @inbounds endpoints_a[1, i] = endpoints_a[1, i] / 2
		@fastmath @inbounds endpoints_a[2, i] = endpoints_a[2, i] + endpoints_b[2, i]
		@fastmath @inbounds endpoints_a[2, i] = endpoints_a[2, i] / 2
	end
      return endpoints_a;
end

# returns the lengths of the edges for all edges - kwargs are parsed to get_edge_endpoints, though globalized should not matter at all
function get_edge_lengths{T}(mesh::Mesh{T}; kwargs...)
  	endpoints_a, endpoints_b = get_edge_endpoints(mesh; kwargs...);
	edgelengths = zeros(T, size(endpoints_a, 2))
	@simd for i in 1:size(endpoints_a, 2)
		@fastmath @inbounds endpoints_a[1, i] = endpoints_a[1, i] - endpoints_b[1, i]
		@fastmath @inbounds endpoints_a[2, i] = endpoints_a[2, i] - endpoints_b[2, i]
	        @fastmath @inbounds edgelengths[i] = norm(endpoints_a[:, i])
	end
      return edgelengths
end
#=
# removes the edge at ind
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
=#

### INIT
function Mesh(index, fixed=false; T=Float64)
	println("Creating mesh for $index")
	# mesh lengths in each dimension
	dists = [PARAMS[:mesh][:mesh_length] * sin(pi / 3); PARAMS[:mesh][:mesh_length]];

	# dimensions of the mesh as a rectangular grid, maximal in each dimension
	# e.g. a mesh with 5-4-5-4-5-4-5-4 nodes in each row will have dims = (8, 5)
	# 1 is added because of 1-indexing (in length)
	# 2 is added to pad the mesh to extend it beyond by one meshpoint in each direction
	dims = round.(Int64, div.(get_image_size(:match_image, mip=get_mip(:match_image)), dists)) + 1 + 2;

	# location of the first node (top left)
	# TODO: Julia does not support rem() for two arrays, so divrem() cannot be used
	topleft_offset = (get_image_size(:match_image, mip=get_mip(:match_image)) .% dists) / 2 - dists;

	n = maximum([get_mesh_index(dims, dims[1], dims[2]); get_mesh_index(dims, dims[1], dims[2]-1)]);
	m = 0;
	m_upperbound = 6 * n;

	src_nodes = Points{T}(2, n);
	edges = spzeros(T, n, m_upperbound);
	
	edges_meshinds = Array{Int64}(m_upperbound)
	edges_edgeinds = Array{Int64}(m_upperbound)
	edges_vals = Array{Float64}(m_upperbound)

	meshindex = fill(0, dims[1], dims[2])
	@fastmath @inbounds for i in 1:dims[1], j in 1:dims[2]
		k = get_mesh_index(dims, i, j); if k == 0 continue; end
		meshindex[i,j] = k
		view(src_nodes,:,k)[:] = get_mesh_coord(dims, topleft_offset, dists, i, j);
		end

	function addedge(m, cur, k, i, j, meshindex, edges_meshinds, edges_edgeinds, edges_vals)
		  edges_meshinds[cur] = k
		  edges_edgeinds[cur] = m
		  edges_vals[cur] = -1.0

		  cur += 1;
		  edges_meshinds[cur] = meshindex[i, j]
		  edges_edgeinds[cur] = m
		  edges_vals[cur] = +1.0
	end

	cur = -1
	@fastmath @inbounds for i in 1:dims[1], j in 1:dims[2]
	  	k = meshindex[i, j]; if k== 0 continue; end
		if (j != 1)
		  	m += 1;
			cur += 2;
			addedge(m, cur, k, i, j-1, meshindex, edges_meshinds, edges_edgeinds, edges_vals);
		end

		if (i != 1)
			if iseven(i) || j != dims[2]
		  	m += 1;
			cur += 2;
			addedge(m, cur, k, i-1, j, meshindex, edges_meshinds, edges_edgeinds, edges_vals);
			end
			if iseven(i) && (j != dims[2]) 			
		  	m += 1;
			cur += 2;
			addedge(m, cur, k, i-1, j+1, meshindex, edges_meshinds, edges_edgeinds, edges_vals);
			end
			if isodd(i) && (j != 1)
		  	m += 1;
			cur += 2;
			addedge(m, cur, k, i-1, j-1, meshindex, edges_meshinds, edges_edgeinds, edges_vals);
			end
			if isodd(i) && ((j == 1) || (j == dims[2]))
		  	m += 1;
			cur += 2;
			addedge(m, cur, k, i-2, j, meshindex, edges_meshinds, edges_edgeinds, edges_vals);
			end
		end
	end
	cur += 1;
	edges_I = edges_meshinds[1:cur];
	edges_J = edges_edgeinds[1:cur];
	edges_V = edges_vals[1:cur];
	#edges = sparse(edges_meshinds[1:cur], edges_edgeinds[1:cur], edges_vals[1:cur], n, m);

	#= legacy code
      @time begin
	@fastmath @inbounds for i in 1:dims[1], j in 1:dims[2]
		k = get_mesh_index(dims, i, j); if k == 0 continue; end
		view(src_nodes,:,k)[:] = get_mesh_coord(dims, topleft_offset, dists, i, j);
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

      end


	@inbounds edges = edges[:, 1:m];
	=#
	dst_nodes = deepcopy(src_nodes);

	properties = Dict{Symbol, Any}(
				    :params => PARAMS,
				    :fixed => fixed);

	return Mesh(index, src_nodes, dst_nodes, edges_I, edges_J, edges_V, properties);
end

function remesh!(mesh::Mesh)
	newmesh = Mesh(mesh.index, mesh.properties[:params], mesh.properties[:fixed]; rotated = true);
	mesh.src_nodes = deepcopy(newmesh.src_nodes);
	mesh.dst_nodes = deepcopy(newmesh.dst_nodes);
	mesh.edges_I = deepcopy(newmesh.edges_I);
	mesh.edges_J = deepcopy(newmesh.edges_J);
	mesh.edges_V = deepcopy(newmesh.edges_V);
	newmesh = 0;
	return mesh
end






# find the triangular mesh indices for a given point in mesh image coordinates
function find_mesh_triangle{T}(mesh::Mesh{T}, point::Union{Point{T}, SubArray{T, 1}})
  @inbounds begin	
	dims, dists = get_dims_and_dists(mesh); 
	@fastmath point_padded = point - get_topleft_offset(mesh);
	@fastmath indices_raw = point_padded ./ dists

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

	node0 = mesh.src_nodes[:, ind0];
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

function find_mesh_triangle{T}(mesh::Mesh{T}, points::Points{T})
	return Triangles(map(find_mesh_triangle, repeated(mesh), columnviews(points)));
end

# Convert Cartesian coordinate to triple of barycentric coefficients
function get_triangle_weights{T}(mesh::Mesh{T}, point::Union{Point{T}, SubArray{T, 1}}, triangle::Triangle)
	if triangle == NO_TRIANGLE::Triangle return NO_WEIGHTS::Weight; end
	@inbounds R = vcat(mesh.src_nodes[:, triangle[1]]', mesh.src_nodes[:, triangle[2]]', mesh.src_nodes[:, triangle[3]]')
	R = hcat(R, ones(T, 3, 1));
	r = vcat(point, 1.0);

	V = r' * R^-1;

	return (V[1], V[2], V[3]);
end

function get_triangle_weights!{T}(mesh::Mesh{T}, point::Union{Point{T}, SubArray{T, 1}}, triangle::Triangle, R, r)

	if triangle == NO_TRIANGLE::Triangle return NO_WEIGHTS::Weight; end

	@inbounds R[1:2, 1] = mesh.src_nodes[:, triangle[1]];
	@inbounds R[1:2, 2] = mesh.src_nodes[:, triangle[2]];
	@inbounds R[1:2, 3] = mesh.src_nodes[:, triangle[3]];
	
	@inbounds r[1:2] = point;

	@fastmath V = R \ r;

	return (V[1], V[2], V[3]);
end

function get_triangle_weights{T}(mesh::Mesh{T}, points::Points{T}, triangles::Triangles)
  	R = ones(T, 3, 3)
  	r = ones(T, 3)

	return Weights(map(get_triangle_weights!, repeated(mesh), columnviews(points), triangles, repeated(R), repeated(r)))
end


function get_tripoint_dst(mesh::Mesh, triangles::Triangles, weights::Weights)
	return hcat(map(get_tripoint_dst, repeated(mesh), triangles, weights)...)
end

function get_tripoint_dst(mesh::Mesh, triangle::Triangle, weight::Weight)
	if triangle == NO_TRIANGLE return NO_POINT; end
	@fastmath @inbounds dst_trinodes = mesh.dst_nodes[:, triangle[1]], mesh.dst_nodes[:, triangle[2]], mesh.dst_nodes[:, triangle[3]]
	@fastmath @inbounds dst_point = dst_trinodes[1] * weight[1] + dst_trinodes[2] * weight[2] + dst_trinodes[3] * weight[3];
	return dst_point;
end

function deform(points::Points, mesh::Mesh)
  triangles = find_mesh_triangle(mesh, points);
	weights = get_triangle_weights(mesh, points, triangles);
	return get_tripoint_dst(mesh, triangles, weights);
end

function fix!(mesh::Mesh)
	mesh.properties[:fixed] = true;
end

function unfix!(mesh::Mesh)
	mesh.properties[:fixed] = false;
end

function is_fixed(mesh::Mesh)
	return mesh.properties[:fixed];
end

function reset!(mesh::Mesh)
	mesh.dst_nodes = deepcopy(mesh.src_nodes)
end

function get_edges(mesh::Mesh)
  return sparse(mesh.edges_I, mesh.edges_J, mesh.edges_V)
end

# broken
#= 
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
=#
