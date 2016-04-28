"""
Ripped from ImageView > navigation.jl
"""
function stop_playing!(state::ImageView.NavigationState)
    if state.timer != nothing
        close(state.timer)
        state.timer = nothing
    end
end

"""
Ripped from ImageView > navigation.jl
"""
function updatet(ctrls, state)
  Tk.set_value(ctrls.editt, string(state.t))
  Tk.set_value(ctrls.scalet, state.t)
  enableback = state.t > 1
  Tk.set_enabled(ctrls.stepback, enableback)
  Tk.set_enabled(ctrls.playback, enableback)
  enablefwd = state.t < state.tmax
  Tk.set_enabled(ctrls.stepfwd, enablefwd)
  Tk.set_enabled(ctrls.playfwd, enablefwd)
end

"""
Ripped from ImageView > navigation.jl
"""
function incrementt(inc, ctrls, state, showframe)
    state.t += inc
    updatet(ctrls, state)
    showframe(state)
end

"""
Ripped from ImageView > navigation.jl
"""
function playt(inc, ctrls, state, showframe)
    if !(state.fps > 0)
        error("Frame rate is not positive")
    end
    stop_playing!(state)
    dt = 1/state.fps
    state.timer = Timer(timer -> stept(inc, ctrls, state, showframe), dt, dt)
end

"""
Ripped from ImageView > navigation.jl
"""
function stept(inc, ctrls, state, showframe)
    if 1 <= state.t+inc <= state.tmax
        incrementt(inc, ctrls, state, showframe)
    else
        stop_playing!(state)
    end
end

"""
Create loop of the image
"""
function start_loop(imgc, img2, fps=6)
  state = imgc.navigationstate
  ctrls = imgc.navigationctrls
  showframe = state -> ImageView.reslice(imgc, img2, state)
  set_fps!(state, fps)

  if !(state.fps > 0)
      error("Frame rate is not positive")
  end
  stop_playing!(state)
  dt = 1/state.fps
  state.timer = Timer(timer -> loopt(ctrls, state, showframe), dt, dt)
end

"""
Higher level call for ImageView outputs
"""
function stop_loop(imgc)
  stop_playing!(imgc.navigationstate)
end

"""
Endlessly repeat forward loop (continuous stept)
"""
function loopt(ctrls, state, showframe)
  inc = 1
  if state.t+inc < 1 || state.tmax < state.t+inc
      state.t = 0
  end
  incrementt(inc, ctrls, state, showframe)
end

"""
Building on ImageView > navigation.jl
"""
function set_fps!(state, fps)
  state.fps = fps
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

function make_image_stack(meshset, section_range=(1:2), slice=(1:200,1:200); 
                                            include_reverse=false, perm=[1,2,3])
  imgs = []
  for mesh in meshset.meshes[section_range]
    index = (mesh.index[1:2]..., mesh.index[3]-1, mesh.index[4]-1)
    img = reinterpret(UFixed8, get_slice(get_path(index), slice))
    push!(imgs, img)
  end
  if include_reverse
    img_stack = cat(3, imgs..., reverse(imgs[2:end-1])...)
  else
    img_stack = cat(3, imgs...)
  end
  return Image(permutedims(img_stack, perm), timedim=3)
end

"""
"""
function mark_stack(mov; fps=10, include_reverse=true)
  e = Condition()

  N = size(mov, 3)
  if include_reverse
    N = round(Int64, (N + 2) / 2)
  end
  frame_errors = zeros(Int64, N)
  escape = false
  paused = false

  imgc, img2 = view(mov, pixelspacing=[1,1])
  state = imgc.navigationstate
  set_fps!(state, fps)
  start_loop(imgc, img2, fps)

  c = canvas(imgc)
  win = Tk.toplevel(c)
  bind(win, "<KP_1>", path->mark_frame(1))
  bind(win, "<KP_2>", path->mark_frame(2))
  bind(win, "<KP_3>", path->mark_frame(3))
  bind(win, "<KP_Enter>", path->destroy())
  bind(win, "<Up>", path->adjust_fps(1))
  bind(win, "<Down>", path->adjust_fps(-1))
  bind(win, "<Escape>", path->end_auto())
  bind(win, "<space>", path->toggle_loop())
  # bind(win, "<Delete>", path->destroy())
  bind(win, "<Destroy>", path->end_count())

  function mark_frame(n)
    f = get_frame()
    frame_errors[f] = n
  end

  function get_frame()
    t = state.t
    if t > N
      t = 2*N - t
    end
    println(t)
    return t
  end

  function adjust_fps(d)
    fps += d
    fps = min(fps, 50)
    fps = max(fps, 1)
    println("Adjust FPS: ", fps)
    stop_loop(imgc)
    set_fps!(state, fps)
    start_loop(imgc, img2, fps)
  end

  function toggle_loop()
    if paused
      start_loop(imgc, img2, fps)
    else
      stop_loop(imgc)
    end
    paused = !paused
  end

  function destroy()
    stop_loop(imgc)
    end_count()
    Tk.destroy(win)
  end

  function end_auto()
    escape = !escape
    println("End auto mode: ", escape)
    destroy()
  end

  function end_count()
    notify(e)
    bind(win, "<KP_1>", path->path)
    bind(win, "<KP_2>", path->path)
    bind(win, "<KP_3>", path->path)
    bind(win, "<KP_Enter>", path->path)
    bind(win, "<Up>", path->path)
    bind(win, "<Down>", path->path)
    bind(win, "<Escape>", path->path)
    bind(win, "<space>", path->path)
    # bind(win, "<Delete>", path->path)
    bind(win, "<Destroy>", path->path)
  end
  
  wait(e)
  return frame_errors, escape, fps
end

"""
Cycle through sections of the stack stack, with images staged for easier viewing
"""
function view_stack(indexA::Index, indexB::Index, slice=(1:200,1:200);
                                        include_reverse=false, perm=[1,2,3])
  imgs = []
  for index in get_index_range(indexA, indexB)
    if is_finished(indexA)
      index = finished(index)
    end
    if is_montaged(indexA)
      index = prealigned(index)
    end
    print(string(join(index[1:2], ",") ,"|"))
    img = reinterpret(UFixed8, get_slice(get_path(index), slice))
    push!(imgs, img)
  end
  return view_stack(imgs; include_reverse=include_reverse, perm=perm)
end

function view_stack(imgs; include_reverse=false, perm=[1,2,3])
  img_stack = permutedims(cat(3, imgs...), perm)
  if include_reverse
    img_stack = cat(3, img_stack, img_stack[:,:,end:-1:1])
  end
  imgc, img2 = view(Image(img_stack, timedim=3))

  c = canvas(imgc)
  win = Tk.toplevel(c)
  fnotify = ImageView.Frame(win)
  lastrow = 2
  ImageView.grid(fnotify, lastrow+=1, 1, sticky="ew")
  xypos = ImageView.Label(fnotify)
  imgc.handles[:pointerlabel] = xypos
  ImageView.grid(xypos, 1, 1, sticky="ne")
  ImageView.set_visible(win, true)
  c.mouse.motion = (path,x,y)-> updatexylabel(xypos, imgc, img2, x, y, slice[2][1], slice[1][1])

  bind(win, "<Up>", path->adjust_fps(imgc, img2, 1))
  bind(win, "<Down>", path->adjust_fps(imgc, img2, -1))
  bind(win, "<Escape>", path->exit_stack(imgc, img2))
  bind(win, "<Destroy>", path->exit_stack(imgc, img2))

  start_loop(imgc, img2, 30)
  return imgs
end

function updatexylabel(xypos, imgc, img2, x, y, x_off, y_off)
  w = width(imgc.c)
  xu, yu = ImageView.device_to_user(Tk.getgc(imgc.c), x, y)
  # Image-coordinates
  xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
  if ImageView.isinside(imgc.canvasbb, x, y)
    val = img2[xi,yi]
    xo, yo = xi + x_off, yi + y_off
    str = "$xo, $yo ($xi, $yi): $val"
    if length(str)*10>w
      ImageView.set_value(xypos, "$xo, $yo ($xi, $yi)")
    else
      ImageView.set_value(xypos, str)
    end
  else
    ImageView.set_value(xypos, "$xo, $yo ($xi, $yi)")
  end
end

function adjust_fps(imgc, img2, d)
  state = imgc.navigationstate
  fps = state.fps
  fps += d
  if 1 <= fps <= 50
    stop_loop(imgc)
    set_fps!(state, fps)
    start_loop(imgc, img2, fps)
  end
end

function exit_stack(imgc, img2)
  println("exit stack")
  stop_loop(imgc)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  bind(win, "<Up>", path->path)
  bind(win, "<Down>", path->path)
  bind(win, "<Destroy>", path->path)
  bind(win, "<Escape>", path->path)
  destroy(win)
end
