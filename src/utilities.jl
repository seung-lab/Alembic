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
function update_offsets(index::Index, offset, sz=[0, 0])
  
  if is_montaged(index) offsets_fp = montaged_offsets_path;
  elseif is_prealigned(index) offsets_fp = prealigned_offsets_path;
  elseif is_aligned(index) offsets_fp = aligned_offsets_path;
  else offsets_fp = premontaged_offsets_path; end

  image_fn = string(get_name(index));

  println("Updating offsets for ", image_fn, " in:\n", offsets_fp)

  if !isfile(offsets_fp)
    f = open(offsets_fp, "w")
    close(f)
    offsets = [image_fn, offset..., sz...]'
  else  
    offsets = readdlm(offsets_fp)
    idx = findfirst(offsets[:,1], image_fn)
    if idx != 0
      offsets[idx, 2:3] = collect(offset)
      if sz != [0, 0]
        offsets[idx, 4:5] = collect(sz)
      end
    else
      offsets_line = [image_fn, offset..., sz...]
      offsets = vcat(offsets, offsets_line')
    end
  end
  offsets = offsets[sortperm(offsets[:, 1], by=parse_name), :];
  writedlm(offsets_fp, offsets)
  
  if is_montaged(index) global MONTAGED_OFFSETS = parse_offsets(offsets_fp);
  elseif is_prealigned(index) global PREALIGNED_OFFSETS = parse_offsets(offsets_fp);
  elseif is_aligned(index) global ALIGNED_OFFSETS = parse_offsets(offsets_fp);
  else global PREMONTAGED_OFFSETS = parse_offsets(offsets_fp);
  end
end

function update_offsets(name, offset, sz=[0, 0])
  update_offsets(parse_name(name), offset, sz);
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
