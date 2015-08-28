# using Julimaps
# using Params
# using MeshModule
# using IO
# using Bounding: BoundingBox, snap_bb, find_mesh_bb

type Tile
	name::String
	affine
	mesh
	id::Index # wafer, section, row, col
end

Tile(name::String) = Tile(name, nothing, nothing, parsename(name))

"""
Create dictionary of waferpaths in bucket from file
"""
function waferpaths2dict(wafer_path_file)
	warray = readdlm(wafer_path_file)
	wdict = Dict()
	for (idx, path) in zip(warray[:,1], warray[:,2])
		wdict[idx] = path
	end
	return wdict
end

const BUCKET = "."
const AFFINE_DIR = "~/seungmount/research/150502_piriform/affine_transforms"
const WAFER_DIR = waferpaths2dict("./datasets/piriform/wafer_paths.txt")

"""
Extract wafer no, section no, row, and col for tile from filename
"""
function parsename(path)
	# m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", path)
	m = match(r"(\d*)-(\d*)_montage", path)
	# return parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2]) # wafer, section, row, column
	return parse(Int, m[1]), parse(Int, m[2]), 0, 0
end

"""
Load original image for tile, using the WAFER_DIR and filename
"""
function load_image(tile::Tile)
	# section_folder = string("S2-W00", tile.id[1], "_Sec", tile.id[2], "_Montage")
	# path = joinpath(homedir(), WAFER_DIR[tile.id[1]], section_folder, string(tile.name, ".tif"))
	# section_folder = string("W00", tile.id[1], "_Sec", tile.id[2])
	path = joinpath(".", "output_images", string("(", tile.id[1], ",", tile.id[2], ")_montage.tif"))
	# path = joinpath(".", "input_images", "sections", string("S2-W00", tile.id[1], "_Sec", tile.id[2], ".tif"))
	println(path)
	return getFloatImage(path)
end

"""
Load affine transformation matrix for provided tile
"""
function load_affine!(tile::Tile)
	path = joinpath(AFFINE_DIR, string(tile.name, ".csv"))
	tile.affine = readcsv(path)
end

"""
Set affine transformation matrix for tile
"""
function set_affine!(tile::Tile, affine::Array{Real,(3,3)})
	tile.affine = affine
end

"""
Set affine transformation matrix for tile
"""
function set_mesh!(tile::Tile, mesh)
	tile.mesh = mesh
end

"""
Parse Mesh object to retrieve nodes, edges, and spatial reference

Args:

* mesh: Mesh object

Returns:

* src_nodes: 2xN array of original mesh nodes
* dst_nodes: 2xN array of mesh nodes after elastic solving
* mesh.edges: MxN array of edges as edge-node incidence matrix 
* offset: global offset of the mesh
"""
function parse_mesh(mesh)
    src_nodes = hcat(mesh.nodes...)
    dst_nodes = hcat(mesh.nodes_t...)
    offset = mesh.disp

    # src_nodes = xy2yx(src_nodes .- offset)
    # dst_nodes = xy2yx(dst_nodes .- offset)
    return src_nodes', dst_nodes', mesh.edges, offset
end


function meshwarp(tile::Tile)
    @assert tile.mesh != nothing
    img = load_image(tile)
	src_nodes, dst_nodes, incidence, offset = parse_mesh(tile.mesh)
    node_dict = incidence2dict(incidence)
    triangles = dict2triangles(node_dict)
    return @time meshwarp(img, src_nodes, dst_nodes, triangles, offset)
end

function get_meshwarp_bb(tile)
	nodes = tile.mesh.nodes_t
	bb = Bounding.snap_bb(Bounding.find_mesh_bb(nodes))
	bb.h += 1
	bb.w += 1
	return bb
end


function load_section(dir_path)
# Returns:
# 	Array{Tile, 1}
end

"""
Create tile array with meshes from mesh_set
"""
function load_tiles(mesh_set)
	tiles = []
	for mesh in mesh_set.meshes
		t = Tile(mesh.name)
		set_mesh!(t, mesh)
		push!(tiles, t)
	end
	return tiles
end
