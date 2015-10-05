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
function points_to_Nx3_matrix(points)
  points = hcat(points...)
  if size(points,1)==2
    points = [points; ones(eltype(points), 1, size(points,2))]'
  end
  return points
end

"""
Edit the offset_log text file associated with prealignment and alignment

log_path: file path of the log file, a .txt file
image_name: string not including the file extension
offset: 2-element collection for the i,j offset
sz: 2-element collection for the i,j height and width
"""
function update_offset_log!(log_path, image_name, offset, sz)
  if !isfile(log_path)
    f = open(log_path, "w")
    close(f)
    offset_log = [image_name, offset..., sz...]'
  else  
    offset_log = readdlm(log_path)
    idx = findfirst(offset_log[:,1], image_name)
    if idx != 0
      offset_log[idx, 2:3] = collect(offset)
      offset_log[idx, 4:5] = collect(sz)
    else
      log_line = [image_name, offset..., sz...]
      offset_log = vcat(offset_log, log_line')
    end
  end
  writedlm(log_path, offset_log)
end

function update_offsets(index::Index, offset, sz=[0, 0])
  
  if is_montaged(index) offsets_fp = montaged_offsets_path;
  elseif is_prealigned(index) offsets_fp = prealigned_offsets_path;
  elseif is_aligned(index) offsets_fp = aligned_offsets_path;
  else offsets_fp = premontaged_offsets_path; end
  
  if !isfile(offsets_fp)
    f = open(offsets_fp, "w")
    close(f)
    offsets = [get_name(index), offset..., sz...]'
  else  
    offsets = readdlm(offsets_fp)
    idx = findfirst(offsets[:,1], get_name(index))
    if idx != 0
      offsets[idx, 2:3] = collect(offset)
      if sz != 0
        offsets[idx, 4:5] = collect(sz)
      end
    else
      offsets_line = [get_name(index), offset..., sz...]
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

function update_offsets(name, offset, sz)
  update_offsets(get_index(name), offset, sz);
end
"""
Create array of alphabetized filenames that have file extension in directory
"""
function sort_dir(dir, file_extension="jld")
    files_in_dir = filter(x -> x[end-2:end] == file_extension, readdir(dir))
    return sort(files_in_dir, by=x->parse_name(x))
end




