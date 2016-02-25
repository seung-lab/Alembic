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
Edit the offset_log text file associated with an index

index: 4-element tuple for section identifier
offset: 2-element collection for the i,j offset
sz: 2-element collection for the i,j height and width
"""
function update_offset(index::Index, offset, sz=[0, 0])
  
  if is_montaged(index) registry_fp = montaged_registry_path;
  elseif is_prealigned(index) registry_fp = prealigned_registry_path;
  elseif is_aligned(index) registry_fp = aligned_registry_path;
  else registry_fp = premontaged_registry_path; end

  image_fn = string(get_name(index));

  println("Updating registry for ", image_fn, " in:\n", registry_fp, ": offset is now ", offset)

  if !isfile(registry_fp)
    f = open(registry_fp, "w")
    close(f)
    registry = [image_fn, offset..., sz...]'
  else  
    registry = readdlm(registry_fp)
    idx = findfirst(registry[:,1], image_fn)
    if idx != 0
      registry[idx, 2:3] = collect(offset)
      if sz != [0, 0]
        registry[idx, 4:5] = collect(sz)
      end
    else
      registry_line = [image_fn, offset..., sz...]
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

function update_offset(name, offset, sz=[0, 0])
  update_offset(parse_name(name), offset, sz);
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


function make_training_data(ms)

       X = Array{Float64, 2}(count_correspondences(ms), 4)
       Y = fill(1, count_correspondences(ms))
       current = 0

       for m in ms.matches
              X[(1:count_correspondences(m)) + current, 1] = get_properties(m, "norm")[:]';
              X[(1:count_correspondences(m)) + current, 2] = get_properties(m, "r_val")[:]';
              X[(1:count_correspondences(m)) + current, 3] = get_properties(m, "src_normalized_dyn_range")[:]';
       	   X[(1:count_correspondences(m)) + current, 4] = get_properties(m, "src_kurtosis")[:]';  
       	   Y[(collect(Int64, get_rejected_indices(m))) + current] = -1
       	   current = current + count_correspondences(m);
       end
       return X, Y
       end

