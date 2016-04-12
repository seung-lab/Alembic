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

function reset_offset(index::Index)
  update_offset(index, [0,0])
end

"""
Edit the offset_log text file associated with an index

index: 4-element tuple for section identifier
offset: 2-element collection for the i,j offset
sz: 2-element collection for the i,j height and width
"""
function update_offset(index::Index, offset, sz=[0, 0], needs_render = false)
  
  if is_montaged(index) registry_fp = montaged_registry_path;
  elseif is_prealigned(index) registry_fp = prealigned_registry_path;
  elseif is_aligned(index) registry_fp = aligned_registry_path;
  else registry_fp = premontaged_registry_path; end

  image_fn = string(get_name(index));

  println("Updating registry for ", image_fn, " in:\n", registry_fp, ": offset is now ", offset)

  if !isfile(registry_fp)
    f = open(registry_fp, "w")
    close(f)
    registry = [image_fn, offset..., sz..., needs_render]'
  else  
    registry = readdlm(registry_fp)
    idx = findfirst(registry[:,1], image_fn)
    if idx != 0
      registry[idx, 2:3] = collect(offset)
      if sz != [0, 0]
        registry[idx, 4:5] = collect(sz)
      end
    else
      registry_line = [image_fn, offset..., sz..., needs_render]
      registry = vcat(registry, registry_line')
    end
  end
  registry = registry[sortperm(registry[:, 1], by=parse_name), :];
  writedlm(registry_fp, registry)
  
  if is_montaged(index) global REGISTRY_MONTAGED = parse_registry(registry_fp);
  elseif is_prealigned(index) global REGISTRY_PREALIGNED = parse_registry(registry_fp);
  elseif is_aligned(index) global REGISTRY_ALIGNED = parse_registry(registry_fp);
  else global REGISTRY_PREMONTAGED = parse_registry(registry_fp);
  end
end

function update_offset(name, offset, sz=[0, 0], needs_render = false)
  update_offset(parse_name(name), offset, sz, needs_render);
end

function expunge(index)
#	image = get_path(index);
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

"""
Create rotation transform and apply to image (degrees)
"""
function imrotate(img, angle)
  tform = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
  return imwarp(img, tform);
end

function make_translation_matrix(offset)
  return [1 0 0; 0 1 0; offset[1] offset[2] 1]
end

function make_scale_matrix(scale)
  return [scale 0 0; 0 scale 0; 0 0 1]
end
