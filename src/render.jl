# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150807
#
# Functions to load meshes, warp images via piecewise affine transforms, and
# display images with meshes.

using JLD
using Images
using ImageView
using Color
using FixedPointNumbers

rawdata(img) = convert(Array{Float64, 2}, data(separate(img)))

"""
`MERGE_IMAGES` - Place images in global reference image

    merged_img, [bb.i, bb.j] = merge_images(imgs, offsets)

* `imgs`: 1D array, images (2D arrays)
* `offsets`: 1D array, 2-element array positions of corresponding image 
	in 'imgs` in global space

""" 
function merge_images(imgs, offsets)
    bbs = []
    for (img, offset) in zip(imgs, offset)
        push!(bbs, BoundingBox(offset..., size(img)...))
    end
    global_ref = sum(bbs)
    merged_img = zeros(global_ref.h, global_ref.w)
    for (idx, (img, bb)) in enumerate(zip(imgs, bbs))
        println(idx)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        merged_img[i:i+h, j:j+w] = max(merged_img[i:i+h, j:j+w], img)
        imgs[idx] = 0
        gc()
    end
    return merged_img, [global_ref.i, global_ref.j]
end

"""
`MERGE_IMAGES` - Using array of Mesh objects
"""
function merge_images(meshes)
    imgs = []
    offsets = []
    for mesh in meshes
        img, offset = meshwarp(mesh)
        push!(imgs, img)
        push!(offsets, offset)
    end
    return merge_images(imgs, offests)
end

"""
`IMFUSE_SECTION` - Stitch section with seam overlays.

INCOMPLETE
"""
function imfuse_section(meshes)
    img_reds = []
    img_greens = []
    bbs_reds = []
    bbs_greens = []
    for mesh in meshes
        img, bb = meshwarp(mesh)
        img = restrict(img)
        bb = BoundingBox(bb/2..., size(img)[1], size(img)[2])
        if (mesh.index[3] + mesh.index[4]) % 2 == 1
            push!(img_reds, img)
            push!(bbs_reds, bb)
        else
            push!(img_greens, img)
            push!(bbs_greens, bb)
        end
    end
    global_ref = sum(bbs_reds) + sum(bbs_greens)
    red_img = zeros(Int(global_ref.h), Int(global_ref.w))
    for (idx, (img, bb)) in enumerate(zip(img_reds, bbs_reds))
        println(idx)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        red_img[i:i+h, j:j+w] = max(red_img[i:i+h, j:j+w], img)
        img_reds[idx] = 0
        gc()
    end
    green_img = zeros(Int(global_ref.h), Int(global_ref.w))
    for (idx, (img, bb)) in enumerate(zip(img_greens, bbs_greens))
        println(idx)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        green_img[i:i+h, j:j+w] = max(green_img[i:i+h, j:j+w], img)
        img_greens[idx] = 0
        gc()
    end
    O, O_bb = imfuse(red_img, [0,0], green_img, [0,0])
    view(make_isotropic(O))
end    

"""
Overlay two images on top of each other using their offsets. Colors one
image red, the other green, and the overlap yellow.
Uses rounded interpolation.

Args:

* A: image A (2D array)
* offset_A: 2-element array for global position of image A
* B: image B (2D array)
* offset_B: 2-element array for global position of image A

Returns:

* O: Image object combining both image A & B

    imfuse(A, offset_A, B, offset_B)
"""
function imfuse(A, offset_A, B, offset_B)
    # pad to common origin
    BB_C = offset_B - offset_A
    if BB_C[1] > 0
        B = padimage(B, 0, BB_C[1], 0, 0)
    elseif BB_C[1] < 0
        A = padimage(A, 0, -BB_C[1], 0, 0)
    end 
    if BB_C[2] > 0
        B = padimage(B, BB_C[2], 0, 0, 0)
    elseif BB_C[2] < 0
        A = padimage(A, -BB_C[2], 0, 0, 0)
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
    BB_O = min(offset_A, offset_B)
    return O, BB_O
end

"""
Pad image exterior to meet new_sz dimensions 

Images are column-major.

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

function padimage(img, offset, bb)
	ylow = offset[1] - bb.i
	xlow = offset[2] - bb.j
	yhigh = bb.h - size(img, 1)
	xhigh = bb.w - size(img, 2)
	xlow *= xlow > 0
	ylow *= ylow > 0
	xhigh *= xhigh > 0
	yhigh *= yhigh > 0
    return padimage(img, xlow, ylow, xhigh, yhigh)
end

function render_montage_for_directory()
    for fn in readdir(MONTAGED_DIR)[4:end]
        if fn[end-2:end] == "jld"
            println("Rendering ", fn[1:end-4])
            meshset = load(joinpath(MONTAGED_DIR, fn))["MeshSet"]
            img, offset = render_section(meshset.meshes)
            img = grayim(img)
            img["spatialorder"] = ["y", "x"]
            println("Writing ", fn[1:end-4])
            @time imwrite(img, joinpath(MONTAGED_DIR, string(fn[1:end-4], ".tif")))
            img = 0
            gc()
        end
    end
end

function render_prealignment_for_directory()
    img_filenames = filter(x -> x[end-2:end] == "tif", readdir(MONTAGED_DIR))
    for filename_A in img_filenames[5:end]
        img_preceding = filter(x -> parseName(x)[2]-1 == parseName(filename_A)[2], img_filenames)
        if length(img_preceding) > 0
            filename_B = img_preceding[1]
            println("Pre-aligning ", filename_B[1:end-4])
            A = getFloatImage(joinpath(MONTAGED_DIR, filename_A))
            B = getFloatImage(joinpath(MONTAGED_DIR, filename_B))
            @time tform = AffineAlignSections(A, B, 1.0)
            A = 0
            gc()
            println("Rendering ", filename_B[1:end-4])
            @time B, B_offset = imwarp(B, inv(tform), [0.0, 0.0])
            println("Writing ", filename_B[1:end-4])
            imwrite(B, joinpath(PRE_ALIGNED_DIR, string(filename_B[1:end-4], "_ij_", B_offset[1], "_", B_offset[2], ".tif")))
        end
    end
end

function render_alignment_for_directory()
	filename = joinpath(ALIGNED_DIR, "1,11-1,12_aligned.jld")
    println("Rendering meshes in ", filename[1:end-4])
    meshset = load(joinpath(ALIGNED_DIR, filename))["MeshSet"]
    bbs = []
    println("Calculating global bounding box")
    for mesh in meshset.meshes
    	nodes = hcat(mesh.nodes_t...)'
    	push!(bbs, snap_bb(find_mesh_bb(nodes)))
    end
    global_bb = sum(bbs)
    for mesh in meshset.meshes
    	println("Warping ", mesh.name)
	    @time img, offset = meshwarp(mesh)
	    img = padimage(img, offset, global_bb)
	    println("Writing ", mesh.name)
	    @time imwrite(img, joinpath(ALIGNED_DIR, string(mesh.name, ".tif")))
	    img = 0
	    gc()
	end
end

function demo_layer_matches()
    fn = "1,1-1,1_alignment"
    meshset = load(joinpath(ALIGNED_DIR, string(fn, ".jld")))["MeshSet"]
    @time img, dst_offset = meshwarp(meshset.meshes[1])
    img = restrict(img)
    src_nodes, dst_nodes = get_matched_points_t(meshset, 1)
    src_pts = hcat(src_nodes...) .- dst_offset
    dst_pts = hcat(dst_nodes...) .- dst_offset
    vectors = vcat(dst_pts, src_pts) / 2
    imgc, img2 = draw_vectors(make_isotropic(img), vectors)
end
