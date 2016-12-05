"""
meshset, area, slice, username, path = load_stack_params("hmcgowan")
review_stack(username, meshset, area, slice, 1, true)
"""
function load_stack_params(username)
  meshset = load((1,2,-3,-3), (1,167,-3,-3))
  area = ImageRegistration.BoundingBox(5000,5000,28000,28000)
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
    push!(d[match_index], push!(indices'[frames .== 1], 0)) #'
    push!(d[match_index], push!(indices'[frames .== 2], 0)) #'
    push!(d[match_index], push!(indices'[frames .== 3], 0)) #'
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
