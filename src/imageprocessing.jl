"""
`MERGE_IMAGES` - Place images in global reference image

    merged_img, [bb.i, bb.j] = merge_images(imgs, offsets)

* `imgs`: 1D array, images (2D arrays)
* `offsets`: 1D array, 2-element array positions of corresponding image 
  in 'imgs` in global space

""" 
function merge_images(imgs, offsets)
    T = typeof(imgs[1][1])
    bbs = []
    for (img, offset) in zip(imgs, offsets)
        push!(bbs, BoundingBox(offset..., size(img)...))
    end
    global_ref = sum(bbs)
    merged_img = zeros(T, global_ref.h, global_ref.w)
    for (idx, (img, bb)) in enumerate(zip(imgs, bbs))
        println("Merging tile ", idx)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        merged_img[i:i+h, j:j+w] = max(merged_img[i:i+h, j:j+w], img)
        imgs[idx] = 0
    end
    return merged_img, [global_ref.i, global_ref.j]
end

"""
`UINT32OVERLAY` - Fuse two UInt8 images together forming UInt32 image for color
"""
function uint32overlay(A::Array{UInt8,2}, B::Array{UInt8,2})
    return convert(Array{UInt32,2},A) .<< 16 + convert(Array{UInt32,2},B) .<< 8
end

"""
`IMFUSE` - Overlay two images on top of each other using their offsets. Colors 
one image red, the other green, and the overlap yellow.
Uses rounded interpolation.

Args:

* A: image A (2D array)
* offset_A: 2-element array for global position of image A
* B: image B (2D array)
* offset_B: 2-element array for global position of image A

Returns:

* O: Image object combining both image A & B

    `imfuse(A, offset_A, B, offset_B)`
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
    O = uint32overlay(A, B)
    # O = Overlay((A,B), (RGB(1,0,0), RGB(0,1,0)))
    BB_O = min(offset_A, offset_B)
    # return O[1:end, 1:end], BB_O
    return O, BB_O
end

"""
`PADIMAGE` - Specify image padding in each of four directions
    
    `new_img = padimage(img, xlow, ylow, xhigh, yhigh)`

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

* new_img: original img, extended with rows and columns of zeros
"""
function padimage{T}(img::Array{T}, xlow::Int64, ylow::Int64, xhigh::Int64, yhigh::Int64)
    h = ylow + size(img, 1) + yhigh
    w = xlow + size(img, 2) + xhigh
    z = zeros(T, h, w)
    z[ylow+1:ylow+size(img,1), xlow+1:xlow+size(img,2)] = img
    return z
end

"""
`RESCOPE` - Crop/pad an image to fill a bounding box
    
    new_img = rescope(img, offset, boundingbox)

Args:

* img: 2D or 3D array
* offset: 2-element array, specifying i & j offset from global origin
* bb: bounding box object in the global reference space

Returns:

* new_img: original img, cropped &/or extended with rows and columns of zeros
"""
function rescopeimage{T}(img::Array{T}, offset, bb)
  z = zeros(T, bb.h+1, bb.w+1)
  imgbb = BoundingBox(offset..., size(img,1)-1, size(img,2)-1)
  xbb = imgbb - bb
  if !isnan(xbb.i) || !isnan(xbb.j) || !isnan(xbb.h) || !isnan(xbb.h)
    crop_img = xbb.i-offset[1]+1 : xbb.i-offset[1]+1+xbb.h, 
                  xbb.j-offset[2]+1 : xbb.j-offset[2]+1+xbb.w
    crop_z = xbb.i-bb.i+1:xbb.i-bb.i+1+xbb.h, xbb.j-bb.j+1:xbb.j-bb.j+1+xbb.w
    z[crop_z...] = img[crop_img...]
  end
  return z
end

"""
`MATCH_PADDING` - For each dimension, pad images to match the larger of the two

  `imgA, imgB = match_padding(imgA, imgB)`
"""
function match_padding(imgA, imgB)
  szA = collect(size(imgA))
  szB = collect(size(imgB))
  szC = max(szA, szB)
  imgA = padimage(imgA, 0, 0, reverse(szC-szA)...)
  imgB = padimage(imgB, 0, 0, reverse(szC-szB)...)
  return imgA, imgB
end