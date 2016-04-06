function clear_all_filters!(indexA, indexB)
        indices = get_index_range(indexA, indexB)
        for ind in indices
                meshset = load(ind)
                for match in meshset.matches
                        clear_filters!(match)
                end
                save(meshset)
        end
end

function print_filter_eval(indexA, indexB)
        indices = get_index_range(indexA, indexB)
        for ind in indices
                meshset = load(ind)
                for (i, match) in enumerate(meshset.matches);
            fp, fn, tp, total = eval_filters(match, [("sigma_5", >, 5, 0)] ,:)
            p = (100 * tp / (fp + tp))
            r = (100 * tp / (fn + tp))
            if !isnan(r) && r < 100.0
                    println(ind, " #", i, " ", p, " / ", r)
            end
       end
   end
end

function update_meshset_with_edits!(meshset, log)
  for i in 1:size(log,1)
    match_index = log[i,5]
    indices_to_remove = [log[i,6]]
    if typeof(indices_to_remove[1]) != Int64
      indices_to_remove = readdlm(IOBuffer(log[i,6]), ',', Int)
    end
    if indices_to_remove[1] != 0
      remove_matches_from_meshset!(meshset, indices_to_remove, match_index)
    end
  end
end

function crop_meshset!(meshset, a, b)
  meshset.indices = meshset.indices[a:b]
  meshset.meshes = meshset.meshes[a:b]
  meshset.N = length(meshset.meshes)
  meshset.nodes_indices = meshset.nodes_indices[a:b]
  meshset.nodes_indices -= meshset.nodes_indices[1]
  matches_mask = map(x -> x>=(a,a) && x<=(b,b), meshset.matches_pairs)
  matches_pairs = meshset.matches_pairs[matches_mask]
  meshset.matches_pairs = map(x->(x[1]-a+1, x[2]-a+1), matches_pairs)
  meshset.matches = meshset.matches[matches_mask]
  meshset.M = length(meshset.matches)

  meshset.n = 0
  meshset.m = 0
  meshset.m_i = 0
  meshset.m_e = 0
  for mesh in meshset.meshes
    meshset.n += mesh.n
    meshset.m += mesh.m
    meshset.m_i += mesh.m

  end
  for matches in meshset.matches
    meshset.m += matches.n
    meshset.m_e += matches.n
  end
end

"""
Append 'EDITED' with a timestamp to update filename
"""
function update_filename(fn, ts="")
  if ts == ""
    ts = Dates.format(now(), "yyyymmddHHMMSS")
  end
  return string(fn[1:end-4], "_EDITED_", ts, fn[end-3:end])
end

"""
Boolean if bounding boxes intersect
"""
function intersects(bbA::BoundingBox, bbB::BoundingBox)
  bb = bbA - bbB
  return !isequal(bb.i, NaN)
end

"""
Given list of bounding boxes, calculate pairs of indices with overlaps
"""
function find_overlaps(boundingboxes)
  bbs = copy(boundingboxes)
  overlap_tuples = []
  i = length(bbs)
  while length(bbs) != 0
    bbi = pop!(bbs)
    for (j, bbj)  in enumerate(bbs)
      if intersects(bbi, bbj)
        push!(overlap_tuples, (i,j))
      end
    end
    i -= 1
  end
  return overlap_tuples
end

function find_boundingboxes(meshset)
  nodes = [get_globalized_nodes_h(mesh)[2]' for mesh in meshset.meshes]
  return map(ImageRegistration.find_mesh_bb, nodes)
end

"""
Crop image with offset to bounding box
"""
function imcrop(img, offset, bb)
  o = zeros(eltype(img), ceil(bb.h)+1, ceil(bb.w)+1)
  ibb = BoundingBox(offset..., size(img)...)
  d = bb - ibb
  o_start = abs(bb.i-d.i)+1:abs(bb.i-d.i) + d.h
  o_end = abs(bb.j-d.j)+1:abs(bb.j-d.j) + d.w
  im_start = abs(ibb.i-d.i)+1:abs(ibb.i-d.i) + d.h
  im_end = abs(ibb.j-d.j)+1:abs(ibb.j-d.j) + d.w
  o[o_start, o_end] = img[im_start, im_end]
  return o
end

"""
`WRITE_SEAMS` - Write out overlays of montaged seams
""" 
function write_seams(meshset, imgs, offsets, indices, flagged_only=true)
  bbs = []
  for (img, offset) in zip(imgs, offsets)
      push!(bbs, BoundingBox(offset..., size(img)...))
  end
  overlap_tuples = find_overlaps(bbs)
  total_seams = flagged_only ? count_flags(meshset) : length(overlap_tuples)
  for (i,j) in overlap_tuples
    src_index, dst_index = indices[i], indices[j]
    k = find_match_index(meshset::MeshSet, src_index, dst_index)
    if !flagged_only || is_flagged(meshset.matches[k])
      println("Writing match #", k, " of ", total_seams, " seams")
      path = get_review_path(src_index, dst_index)
      try 
        img, fuse_offset = imfuse(imgs[i], offsets[i], imgs[j], offsets[j])
        bb = bbs[i] - bbs[j]
        img_cropped = imcrop(img, fuse_offset, bb)
        f = h5open(path, "w")
        chunksize = min(50, min(size(img_cropped)...))
        @time f["img", "chunk", (chunksize,chunksize)] = img_cropped
        f["offset"] = [bb.i, bb.j]
        f["scale"] = 1.0
        close(f)
      catch e
        idx = (indices[i], indices[j])
        log_render_error(MONTAGED_DIR, idx, e)
      end
    end
  end
end