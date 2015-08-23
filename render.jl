# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150807
#
# Functions to load meshes, warp images via piecewise affine transforms, and
# display images with meshes.

using Julimaps
include("incidence2triangles.jl")
include("meshwarp.jl")
include("visualize.jl")
include("Tile.jl")

using JLD
using Images
using ImageView
using Color
using FixedPointNumbers
using MeshModule
using Base.Test

const BUCKET = "."

xy2yx(r) = [r[:,2] r[:,1]]
xy2ij(nodes, height) = [height .- nodes[:, 2] nodes[:, 1]]
rawdata(img) = convert(Array{Float64, 2}, data(separate(img)))

"""
Load matches variables from JLD mesh file generated by Mesh.jl Mesh2JLD

Args:

* mesh_set: MeshSet object loaded from JLD mesh file
* idx: index of mesh matches in mesh_set

Returns:

* src_pts: 2xN array of original mesh nodes
* dst_pts: 2xN array of corresponding block matches
"""
function load_matches(matches)
    src_offset = convert(Array{Int64,1}, matches.src_mesh.disp)
    src_nodes = hcat(matches.src_mesh.nodes...) .- src_offset
    src_idx = matches.src_pointIndices
    src_pts = src_nodes[:,src_idx]
    dst_offset = convert(Array{Int64,1}, matches.dst_mesh.disp)
    dst_pts = hcat(matches.dst_points...) .- dst_offset
    return src_pts, dst_pts
end

"""
Load matches variables from JLD mesh file and calcualte vectors

Args:

* mesh_set: MeshSet object loaded from JLD mesh file
* idx: index of mesh matches in mesh_set

Returns:

* src_pts: 2xN array of original mesh nodes
* dst_pts: 2xN array of corresponding block matches
"""
function load_vectors(matches)
    dst_offset = convert(Array{Int64,1}, matches.dst_mesh.disp)
    src_nodes = hcat(matches.src_mesh.nodes...) .- dst_offset
    src_idx = matches.src_pointIndices
    src_pts = src_nodes[:,src_idx]
    dst_pts = hcat(matches.dst_points...) .- dst_offset
    return vcat(dst_pts, src_pts)
end 

function meshwarp(img::Array, mesh)
    return mesh_warp(img, parse_mesh(mesh)...)
end

function meshwarp(path::String, mesh)
    img = rawdata(imread(path))
    return mesh_warp(img, parse_mesh(mesh)...)
end

function meshwarp(tile::Tile)
    @assert tile.mesh != nothing
    img = load_image(tile)
    return mesh_warp(img, parse_mesh(tile.mesh)...)
end

function bounds2padding(sz, xlow, ylow, xhigh, yhigh)
    xlow *= xlow < 0
    ylow *= ylow < 0
    xlow = abs(xlow)
    ylow = abs(ylow)
    xhigh = xhigh-sz[2]
    yhigh = yhigh-sz[1]
    xhigh *= xhigh > 0
    yhigh *= yhigh > 0
    return [Int(xlow), Int(ylow)], [Int(xhigh), Int(yhigh)]
end

function test_bounds2padding()
    bounds = bounds2padding((10, 10), (3, 1, 9, 11)...)
    @test bounds == ([0, 0], [0, 1])

    bounds = bounds2padding((10, 10), (-1, -1, 1, 1)...)
    @test bounds == ([1, 1], [0, 0])
end

"""
Apply piecewise affine to an image based on a mesh deformation

Args:

* img: 2D array
* mesh: Mesh object with src nodes and dst nodes

Returns:

* img_warped: 2D array with piecewise affine
* bb: BoundingBox object indicating image location in global space

    img_warped, bb = function mesh_warp(img, mesh)
"""
function mesh_warp(img, initial_nodes, final_nodes, incidence, offset)

    println(initial_nodes)
    println(final_nodes)

    bb = find_mesh_bb(final_nodes)
    low, high = bounds2padding(size(img), minsandmax(bb)...)
    img_padded = padimage(img, low..., high...)

    println(low, ", ", high)

    src_nodes = xy2yx((initial_nodes .+ low)')
    dst_nodes = xy2yx((final_nodes .+ low)')

    node_dict = incidence2dict(incidence)
    triangles = dict2triangles(node_dict)

    println(src_nodes)
    println(dst_nodes)

    img_warped, bb = meshwarp(img, src_nodes, dst_nodes, triangles)
    # imgc = draw_mesh(make_isotropic(img_warped), dst_nodes, node_dict)
    # bb = BoundingBox(offset-low..., reverse(size(img_warped))...)
    return img_warped, bb
end

"""
Stitch array of tiles into one image, with max pixel blending
"""
function render_section(tiles)
    tile_imgs = []
    bbs = []
    for tile in tiles
        img, bb = meshwarp(tile)
        push!(tile_imgs, img)
        push!(bbs, bb)
    end
    global_ref = sum(bbs)
    section_img = zeros(global_ref.w, global_ref.h)
    for (idx, (img, bb)) in enumerate(zip(tile_imgs, bbs))
        println(idx)
        i = bb.x - global_ref.x+1
        j = bb.y - global_ref.y+1
        w = bb.w-1
        h = bb.h-1
        section_img[i:i+w, j:j+h] = max(section_img[i:i+w, j:j+h], img')
        tile_imgs[idx] = 0
        gc()
    end
    return section_img
end

# function imfuse_section(tiles)
# # Stitch together Tile objects
#     tile_reds = []
#     tile_greens = []
#     bbs = []
#     for tile in tiles
#         img, bb = mesh_warp_tile(tile)
#         if 
#         push!(tile_imgs, img)
#         push!(bbs, bb)
#     end
#     global_ref = sum(bbs)
#     section_red = zeros(global_ref.w, global_ref.h)
#     for (idx, (img, bb)) in enumerate(zip(tile_imgs, bbs))
#         println(idx)
#         i = bb.x - global_ref.x+1
#         j = bb.y - global_ref.y+1
#         w = bb.w-1
#         h = bb.h-1
#         section_img[i:i+w, j:j+h] = max(section_img[i:i+w, j:j+h], img')
#         tile_imgs[idx] = 0
#         gc()
#     end
#     return section_img
# end    

function resample_to_new_bb(img, BB)
# Shift image to pixel grid of its global space (get an integer spatial ref)
end

"""
Overlay two images on top of each other using their offsets. Colors one
image red, the other green, and the overlap yellow.
Uses rounded interpolation.

Args:

* A: image A (2D array)
* BB_A: spatial reference of image A (2 element offset vector)
* B: image B (2D array)
* BB_B: spatial reference of image B (2 element offset vector)

Returns:

* O: Image object combining both image A & B

    imfuse(A, BB_A, B, BB_B)
"""
function imfuse(A, BB_A, B, BB_B)
    # pad to common origin
    BB_C = BB_B - BB_A
    if BB_C[1] > 0
        B = padimage(B, BB_C[1], 0, 0, 0)
    elseif BB_C[1] < 0
        A = padimage(A, -BB_C[1], 0, 0, 0)
    end 
    if BB_C[2] > 0
        B = padimage(B, 0, BB_C[2], 0, 0)
    elseif BB_C[2] < 0
        A = padimage(A, 0, -BB_C[2], 0, 0)
    end 
    # pad to match sizes
    szA = collect(size(A))
    szB = collect(size(B))
    szC = szB - szA
    if szC[1] > 0
        A = padimage(A, 0, 0, 0, szC[1])
    elseif szC[1] < 0
        B = padimage(B, 0, 0, 0, -szC[1])
    end 
    if szC[2] > 0
        A = padimage(A, 0, 0, szC[2], 0)
    elseif szC[2] < 0
        B = padimage(B, 0, 0, -szC[2], 0)
    end
    O = Overlay((A,B), (RGB(1,0,0), RGB(0,1,0)))
    BB_O = min(BB_A, BB_B)
    return O, BB_O
end

"""
Pad image exterior to meet new_sz dimensions  

     _________________________________  
    |                                 |  
    |             ylow                |  
    |         ______________          |  
    |        |              |         |  
    |  xlow  |     img      |  xhigh  |  
    |        |              |         |  
    |        |______________|         |  
    |                                 |  
    |             yhigh               |  
    |_________________________________|  

Args:

* img: 2D or 3D array
* xlow: amount to pad in x prior to the image
* ylow: amount to pad in y prior to the image
* xhigh: amount to pad in x after to the image
* yhigh: amount to pad in y after to the image

Returns:

* img: original img, extended with rows and columns of zeros

    padimage(img, xlow::Int64, ylow::Int64, xhigh::Int64, yhigh::Int64)
"""
function padimage(img, xlow::Int64, ylow::Int64, xhigh::Int64, yhigh::Int64)
    sz = size(img)
    img = vcat(img, zeros(yhigh, sz[2]))
    sz = size(img)
    img = hcat(img, zeros(sz[1], xhigh))
    sz = size(img)
    img = vcat(zeros(ylow, sz[2]), img)
    sz = size(img)
    img = hcat(zeros(sz[1], xlow), img)
    return img
end

function demo_meshwarp()
# Demo the updated meshwarp function that runs faster than original package
    img = imread(joinpath(BUCKET, "test_images", "turtle.jpg"))
    img = convert(Array{Float64, 3}, data(separate(img)))[:,:,1]
    initial_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    320.0 290.0]
    final_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    400.0 460.0]
    incidence = [1 1 1 0 0 0 0 0;
                -1 0 0 1 1 0 0 0;
                0 0 0 -1 0 1 1 0;
                0 -1 0 0 0 0 -1 1;
                0 0 -1 0 -1 -1 0 -1]
    triangles = [1 2 5;
                1 4 5;
                2 3 5;
                3 4 5];
    src_nodes = xy2yx(initial_nodes')
    dst_nodes = xy2yx(final_nodes')
    node_dict = incidence2dict(incidence)
    draw_mesh(img, src_nodes, node_dict)
    println(size(img))

    warp = meshwarp(img, src_nodes, dst_nodes, triangles)
    draw_mesh(warp, dst_nodes, node_dict)
    println(size(warp))
end

function demo_section()
    mesh_set = load(joinpath(BUCKET, "input_images", "W001_sec21_100_0.001_50000_100.jld"))["MeshSet"]
    matches = mesh_set.matches
    img_set = []
    for match in mesh_set.matches
        src_mesh =  match.src_mesh
        dst_mesh =  match.dst_mesh
        if src_mesh.index[3] > dst_mesh.index[3]
            tile_pathA = joinpath(BUCKET, "input_images", "W001_sec21", src_mesh.path[14:end])
            A, BB_A = meshwarp(tile_pathA, src_mesh) 

            tile_pathB = joinpath(BUCKET, "input_images", "W001_sec21", dst_mesh.path[14:end])
            B, BB_B = meshwarp(tile_pathA, dst_mesh) 

            O, BB_O = imfuse(A, BB_A, B, BB_B)
            imwrite(O, joinpath("/usr/people/tmacrina/Desktop/W001_sec21", string(src_mesh.path[14:end-4],"_",dst_mesh.path[14:end-4],".jpg") )) 
        end
    end
end

function demo_two_tiles_from_mesh_set()
# Load JLD file for MeshSet, warp images, then display
    mesh_set = load(joinpath(BUCKET, "input_images", "Test.jld"))["MeshSet"]
    # r4_c2
    tile_pathA = joinpath(BUCKET, "input_images", "W001_sec20", "Tile_r4-c2_S2-W001_sec20.tif")
    A, BB_A = meshwarp(tile_pathA, mesh_set.meshes[1]) 
    
    # r4_c3
    tile_pathB = joinpath(BUCKET, "input_images", "W001_sec20", "Tile_r4-c3_S2-W001_sec20.tif")
    B, BB_B = meshwarp(tile_pathB, mesh_set.meshes[2])  

    src_pts, dst_pts = load_matches(mesh_set.matches[1])
    # draw_points(make_isotropic(rawdata(imread(tile_pathA))), src_pts)
    # imgc, img2, annotation = draw_points(make_isotropic(rawdata(imread(tile_pathB))), dst_pts)
    O, BB_O = imfuse(A, BB_A, B, BB_B) 
    view(make_isotropic(O))
    return A, BB_A, B, BB_B
end

function demo_render_section()
    fn = "solvedMesh(1, 21,0)_1E-3_1E-3_.5_1E-7"
    mesh_set = load(joinpath(BUCKET, "input_images", string(fn, ".jld")))["MeshSet"]
    tiles = load_tiles(mesh_set)
    section_img = render_section(tiles)
    img = grayim(section_img)
    imwrite(img, joinpath(BUCKET, "output_images", string(fn, "2.jpg")))
end

function test_padimage()
    o = ones(5,2)
    po = padimage(o, 0, 0, 0, 0)
    @test size(po, 1) == 5
    @test size(po, 2) == 2

    o = ones(5,2)
    po = padimage(o, 1, 2, 3, 4)
    @test size(po, 1) == 11
    @test size(po, 2) == 6
    @test po[1,1] == 0

    o = ones(5,2)
    po = padimage(o, 0, 1, 1, 0)
    @test size(po, 1) == 6
    @test size(po, 2) == 3
    @test po[6,3] == 0    

    o = convert(Array{Int64, 2}, ones(5,2))
    po = padimage(o, 1, 2, 3, 4)
    @test size(po, 1) == 11
    @test size(po, 2) == 6
end

function test_imfuse()
    A = rand(5,5)
    B = rand(5,5)
    BB_A = [0, 0]
    BB_B = [0, 0]
    O, BB_O = imfuse(A, BB_A, B, BB_B)
    @test O.channels[1] == A
    @test O.channels[2] == B
    @test BB_O == BB_A
    @test BB_O == BB_B
end
