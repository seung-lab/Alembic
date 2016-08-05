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
  return imwarp(img, tform; kwargs...);
end

function make_rotation_matrix(angle)
  angle = deg2rad(angle)
  rot = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1]
  rot[abs(rot) .< eps] = 0
  return rot;
end

function make_translation_matrix(offset)
  return [1 0 0; 0 1 0; offset[1] offset[2] 1]
end

function make_scale_matrix(scale)
  return [scale 0 0; 0 scale 0; 0 0 1]
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



