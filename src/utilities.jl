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

"""
Transform Nx3 pts by 3x3 tform matrix
"""
function warp_pts(tform, pts)
    pts = hcat(pts, ones(size(pts,1)))
    tpts = pts * tform
    return tpts[:,1:2]
end

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
function imscale(img, scale_factor)
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

function imscale_register_a!(img, scale_factor)
  tform = [scale_factor 0 0; 0 scale_factor 0; 0 0 1];
  bb = BoundingBox{Float64}(0,0, size(img, 1), size(img, 2))
  wbb = tform_bb(bb, tform)
  tbb = snap_bb(wbb)
  if size(REGISTER_A) != (tbb.h, tbb.w)
   global REGISTER_A = zeros(eltype(img), tbb.h, tbb.w)
   end
  return ImageRegistration.imwarp!(REGISTER_A, img, tform);
end

function imscale_register_b!(img, scale_factor)
  tform = [scale_factor 0 0; 0 scale_factor 0; 0 0 1];
  bb = BoundingBox{Float64}(0,0, size(img, 1), size(img, 2))
  wbb = tform_bb(bb, tform)
  tbb = snap_bb(wbb)
  if size(REGISTER_B) != (tbb.h, tbb.w)
   global REGISTER_B = zeros(eltype(img), tbb.h, tbb.w)
   end
  return ImageRegistration.imwarp!(REGISTER_B, img, tform);
end

"""
Create rotation transform and apply to image (degrees)
"""
function imrotate(img, angle)
  tform = make_rotation_matrix(angle)
  return imwarp(img, tform);
end

function user_approves(m="Are you sure?")
  println(m)
  println("(type 'yes')")
  a = readline()
  return chomp(a) == "yes"
end
#=
function flatten_dict(dict::Dict)
	for key in keys(dict)
		if typeof(key) :< Dict

	end
end=#



