"""
Count number of displacement vector lengths k-sigma from the group mean

Args:

* vectors: 4xN array of points (start and end points of displacement vectors)
* k: number of sigmas to set the threshold for counting

Returns:

* integer for number of displacement vectors above the k-sigma threshold

  num = count_outliers(vectors, k)
"""
function count_outliers(vectors, k)
  d = sum((vectors[1:2,:] - vectors[3:4,:]).^2, 1).^(1/2)
  d_mean = mean(d)
  d_std = mean((d.-d_mean).^2).^(1/2)
  return sum((d.-d_mean)./d_std .> k)
end

# """
# Transform Nx3 pts by 3x3 tform matrix
# """
# function warp_pts(tform, pts)
#     pts = hcat(pts, ones(size(pts,1)))
#     tpts = pts * tform
#     return tpts[:,1:2]
# end

"""
Make N array of 1x2 points into Nx3 homogenous points
"""
function homogenize_points(pts)
  pts_new = hcat(pts...)
  pts_new = [pts_new; ones(eltype(pts_new), 1, size(pts_new,2))]'
  return pts_new;
end

"""
Create array of alphabetized filenames that have file extension in directory
"""
function sort_dir(dir, file_extension="jld")
    files_in_dir = filter(x -> x[end-2:end] == file_extension, readdir(dir))
    return sort(files_in_dir, by=x->parse_name(x))
end

"""
Create scaling transform and apply to image
"""
function imscale(img, scale_factor; kwargs...)
  tform = [scale_factor 0 0; 0 scale_factor 0; 0 0 1];
  return imwarp(img, tform; parallel=true, kwargs...);
end
#=
function imscale!(result, img, scale_factor)
  tform = [scale_factor 0 0; 0 scale_factor 0; 0 0 1];
  bb = BoundingBox{Float64}(offset..., size(img, 1), size(img, 2))
  sbb = scale_bb(stack_bb, scale)
  scaled_size = sbb.h, sbb.w
  warped_img = zeros(T, tbb.h, tbb.w)

  return imwarp!(result, img, tform);
end=#

"""
Create rotation transform and apply to image (degrees)
"""
function imrotate(img, angle; kwargs...)
  tform = make_rotation_matrix(angle)
  return imwarp(img, tform; kwargs...)[1];
end

function make_rotation_matrix_from_index(index::FourTupleIndex)
	rotation = get_rotation(index);
	image_size = get_image_size(index);
	return make_rotation_matrix(rotation, image_size);
end

"""
Read in HDF5 image for contrast adjustment in blocks

Hard-coded sample slice
"""
function adjust_contrast_inplace(fn)
  f = h5open(fn, "r+")
  dset = f["img"]
  maxval = 1.0
  minval = 1.0/255

  # generate sample for hist sample
  smp = dset[10000:10:20000,10000:10:20000]
  hist = nquantile(smp[:] / 255, 20);
  maxin = hist[end-1];
  minin = hist[2];

  stride = 256
  for i in 1:stride:size(dset,1)
    for j in 1:stride:size(dset,2)
      if i + stride - 1 <= size(dset,1)
        i_slice = range(i,stride)
      else
        i_slice = colon(i,size(dset,1))
      end
      if j + stride - 1 <= size(dset,2)
        j_slice = range(j,stride)
      else
        j_slice = colon(j,size(dset,2))
      end
      print((i_slice, j_slice))
      img = dset[i_slice, j_slice]
      img = img / 255
      img = min(maxval, max(minval, (img-minin) / (maxin-minin)))
      img = Array{UInt8, 2}(round(UInt8, img*255))
      dset[i_slice, j_slice] = img
    end
  end
  close(f)
end

"""
Create upsampled image, using whole number factors only.

Since imscale is based on imwarp, the center of the pixel movements 
introduces interpolation and border issues for pure-pixel operations. 
"""
function upsample{T}(src_img::SharedArray{T,2}, factor::Int64; parallel=true)
  dst_img = SharedArray{T}(size(src_img,1)*factor, size(src_img,2)*factor)

  if parallel && nprocs() != 1 
    @sync for p in procs()
      j_range = proc_range(p, 1:size(dst_img, 2))
      @async remotecall_wait(upsample_columns!, p, src_img, dst_img, factor, j_range)
    end
  else
    upsample_columns!(src_img, dst_img, factor, 1:size(dst_img,2))
  end

  return dst_img
end

function upsample_columns!{T}(src_img::SharedArray{T,2}, dst_img::SharedArray{T,2}, factor::Int64, j_range::UnitRange{Int64})
  if length(j_range) > 0
    x = 0
    for i in 1:size(dst_img,1)
      if i % factor == 1
        x += 1
      end
      y = floor(Int64, j_range[1]/factor)
      for j in j_range
        if j % factor == 1
          y += 1
        end
        @inbounds dst_img[i,j] = src_img[x,y]
      end
    end
  end
end

function proc_range(idx, arr)
  worker_procs = setdiff(procs(), myid());
  nchunks = length(worker_procs);
  if nchunks == 0 return 1:length(arr); end
  if idx == myid() return 1:0; end
  splits = [round(Int64, s) for s in linspace(0, length(arr), nchunks + 1)];
  return splits[findfirst(worker_procs, idx)]+1:splits[findfirst(worker_procs, idx) + 1]
end