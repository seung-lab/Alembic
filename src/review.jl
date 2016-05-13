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
  return !isnan(bb.i)
end

"""
Get bounding box offset
"""
function get_offset(bb::BoundingBox)
  return [bb.i, bb.j]
end

"""
Shift bounding box by 2-element array
"""
function translate_bb(bb::BoundingBox, offset)
  return BoundingBox(bb.i + offset[1], bb.j + offset[2], bb.h, bb.w)
end

"""
Convert bounding box to tuple of ranges for easy array slicing
"""
function bb_to_slice(bb::BoundingBox)
  return Int64(bb.i+1) : Int64(bb.i+bb.h), Int64(bb.j+1) : Int64(bb.j+bb.w)
end

"""
Convert tuple of ranges to bounding box
"""
function slice_to_bb(slice)
  return BoundingBox(slice[1][1], slice[2][1], slice[1][end]-slice[1][1], slice[2][end]-slice[2][1])
end

"""
Is point contained within the bounding box (border included)?
"""
function point_is_contained(bb::BoundingBox, pt::Point)
  return (bb.i <= pt[1] <= (bb.i+bb.h)) && (bb.j <= pt[2] <= (bb.j+bb.w))
end

"""
At least one point of line is contained within bounding box (border included)?
"""
function line_is_contained(bb::BoundingBox, line)
  return point_is_contained(bb, line[1:2]) || point_is_contained(bb, line[3:4])
end

"""
Given list of bounding boxes, calculate pairs of indices with overlaps, symmetric
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
        push!(overlap_tuples, (j,i))
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
  overlap_tuples = find_overlaps(bbs) # could include tag for asymmetric list
  total_seams = flagged_only ? count_flags(meshset) : length(overlap_tuples)
  for (i,j) in overlap_tuples
    src_index, dst_index = indices[i], indices[j]
    k = find_match_index(meshset, src_index, dst_index)
    if !flagged_only || is_flagged(meshset.matches[k])
      println("Writing match #", k, " of ", total_seams, " seams")
      path = get_review_path(src_index, dst_index)
      img, fuse_offset = imfuse(imgs[i], offsets[i], imgs[j], offsets[j])
      bb = bbs[i] - bbs[j]
      img_cropped = imcrop(img, fuse_offset, bb)
      f = h5open(path, "w")
      chunksize = min(50, min(size(img_cropped)...))
      @time f["img", "chunk", (chunksize,chunksize)] = img_cropped
      f["offset"] = [bb.i, bb.j]
      f["scale"] = 1.0
      close(f)
    end
  end
end

"""
meshset, area, slice, username, path = load_stack_params("hmcgowan")
review_stack(username, meshset, area, slice, 1, true)
"""
function load_stack_params(username)
  meshset = load((1,2,-3,-3), (1,167,-3,-3))
  area = BoundingBox(5000,5000,28000,28000)
  slice = [400, 400]
  path = get_stack_errors_path(meshset, username)
  return meshset, area, slice, username, path
end

function get_stack_errors_path(meshset, username)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[end].index
  fn = string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","),
                "_aligned_stack_errors.txt")
  fn = update_filename(fn, username)
  return joinpath(INSPECTION_DIR, fn)
end

function review_stack(username, meshset, area, slice, k; auto=false, fps=12)
  mov, slice_range = go_to(meshset, area, slice, k; include_reverse=true)
  println("Reviewing stack @ column ", k)
  errors, escape, fps = mark_stack(mov; fps=fps, include_reverse=true)
  path = get_stack_errors_path(meshset, username)
  store_stack_errors(path, username, slice_range, k, errors)
  println("Last reviewed stack @ column ", k)
  if auto & !escape
    return review_stack(username, meshset, area, slice, k+1; auto=true, fps=fps)
  end
end

"""
Stores all slice reviews in chronological order - no overwriting
"""
function store_stack_errors(path, username, slice_range, k, errors)
  ts = Dates.format(now(), "yymmddHHMMSS")
  i, j = slice_range[1][1], slice_range[2][1]
  n, m = slice_range[1][end]-slice_range[1][1], 
                    slice_range[2][end]-slice_range[2][1]
  error_line = [ts, username, i, j, n, m, k, join(errors, ",")]'
  if !isfile(path)
    f = open(path, "w")
    close(f)
    stack_errors = error_line
  else  
    stack_errors = readdlm(path)
    stack_errors = vcat(stack_errors, error_line)
  end
  stack_errors = stack_errors[sortperm(stack_errors[:, 3]), :]
  println("Saving stack_errors:\n", path)
  writedlm(path, stack_errors)
end

"""
Retrieves all slice ranges, calculating count for last review
"""
function get_stack_errors(path, area)
  s = 0.1
  z = 0
  if isfile(path)
    stack_errors = readdlm(path)
    z = zeros(Int64, round(Int64, area.h*s), round(Int64, area.w*s))
    for k in 1:size(stack_errors, 1)
      i = round(Int64, (stack_errors[k, 3] - area.i)*s)+1
      j = round(Int64, (stack_errors[k, 4] - area.j)*s)+1
      iz = round(Int64, stack_errors[k, 5]*s)
      jz = round(Int64, stack_errors[k, 6]*s)
      errors = stack_errors[k, 8]
      if typeof(errors) != Int64
        errors = readdlm(IOBuffer(stack_errors[k, 8]), ',', Int)
      end
      l = length(errors)
      z[i:i+iz, j:j+jz] = ones(Int64, iz+1, jz+1)*l
    end
  end
  return z
end

function get_stack_errors_groundtruth_path()
  fn = "1,2-1,167_aligned_stack_errors_EDITED_tmacrina_baseline.txt"
  return joinpath(inspection_storage_path, fn)
end

function print_stack_errors_report(meshset, path)
  dC = compare_stack_errors(meshset, path)
  report = ["k" "1_agree" "1_disagree" "2_agree" "2_disagree" "3_agree" "3_disagree"]
  for k in sort(collect(keys(dC)))
    agree1 = join(push!(dC[k][1],0), ",")
    disagree1 = join(push!(dC[k][2], dC[k][3]..., 0), ",")
    agree2 = join(push!(dC[k][4],0), ",")
    disagree2 = join(push!(dC[k][5], dC[k][6]..., 0), ",")
    agree3 = join(push!(dC[k][7],0), ",")
    disagree3 = join(push!(dC[k][8], dC[k][9]..., 0), ",")
    report = vcat(report, [k agree1 disagree1 agree2 disagree2 agree3 disagree3])
  end
  path = string(path[1:end-4], "_report.txt")
  println("Saving report:\n", path)
  writedlm(path, report)
  return report
end

function compare_stack_errors(meshset, pathA, pathB=get_stack_errors_groundtruth_path())
  dC = Dict()
  dA = dict_of_stack_errors(meshset, pathA)
  dB = dict_of_stack_errors(meshset, pathB)
  sections = intersect(Set(keys(dB)), Set(keys(dA)))
  for k in sections
    assert(dA[k][4] == dB[k][4])
    A1, A2, A3 = Set(dA[k][1]), Set(dA[k][2]), Set(dA[k][3])
    B1, B2, B3 = Set(dB[k][1]), Set(dB[k][2]), Set(dB[k][3])
    # [TP in A, TN in A, FP in A, FN in A] # TN: match properly removed
    dC[k] = [intersect(A1, B1),
              setdiff(A1, B1),
              setdiff(B1, A1),
              intersect(A2, B2),
              setdiff(A2, B2),
              setdiff(B2, A2),
              intersect(A3, B3),
              setdiff(A3, B3),
              setdiff(B3, A3)]
  end
  return dC
end

function dict_of_stack_errors(meshset, path)
  d = Dict()
  pts = readdlm(path)
  indices = 1:length(meshset.meshes)
  for i in 1:size(pts,1)
    match_index = pts[i,7]
    frames = readdlm(IOBuffer(pts[i,8]), ',', Int)
    d[match_index] = []
    push!(d[match_index], push!(indices'[frames .== 1], 0))
    push!(d[match_index], push!(indices'[frames .== 2], 0))
    push!(d[match_index], push!(indices'[frames .== 3], 0))
    push!(d[match_index], pts[i,4:7])
  end
  return d
end

function normalize_to_uint8(a)
  assert(minimum(a) == 0)
  mx = maximum(a)
  a /= mx
  return convert(Array{UInt8}, round(a*255))
end

function create_stack_colormap()
  return vcat(linspace(RGB(0.0,0.0,0.0), RGB(0.2,0.2,0.2), 127), 
              linspace(RGB(0.2,0.2,0.2), RGB(1.0,0.0,0.0), 128))
end

function view_errors(path, area)
  z = get_stack_errors(path, area)
  a = normalize_to_uint8(z)
  cm = create_stack_colormap()
  b = apply_colormap(a, cm)
  imgc, img2 = view(b, pixelspacing=[1,1])
  c = canvas(imgc)
  win = Tk.toplevel(c)
  fnotify = ImageView.Frame(win)
  lastrow = 1
  ImageView.grid(fnotify, lastrow+=1, 1, sticky="ew")
  xypos = ImageView.Label(fnotify)
  imgc.handles[:pointerlabel] = xypos
  ImageView.grid(xypos, 1, 1, sticky="ne")
  ImageView.set_visible(win, true)
  c.mouse.motion = (path,x,y)-> updatexylabel(xypos, imgc, z-1, x, y)
end

function go_to(meshset, area, slice, k; include_reverse=false)
  assert(k != 0)
  n, m = round(Int64, area.h/slice[1]), round(Int64, area.w/slice[2])
  i = (k-1)%n + 1
  j = ceil(Int64, k/n)
  section_range = 1:length(meshset.meshes)
  islice = ((i-1)*slice[1]:i*slice[1]) + area.i
  jslice = ((j-1)*slice[1]:j*slice[1]) + area.j
  stack = make_image_stack(meshset, section_range, (islice, jslice); 
                                            include_reverse=include_reverse)
  return stack, (islice, jslice)
end

function get_frame(imgc)
  state = imgc.navigationstate
  t = state.t
  if t > N
    t = 2*N - t
  end
  println(t)
  return t
end