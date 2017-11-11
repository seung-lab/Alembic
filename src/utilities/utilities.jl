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
  return imwarp(img, tform);
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



    

