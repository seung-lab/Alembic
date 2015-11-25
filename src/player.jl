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

function go_up_slice(meshset, section_range, slice)
  new_slice = (slice[1][1]-length(slice[1]):slice[1][1], slice[2])
  view_stack(meshset, section_range, new_slice)
  return new_slice
end

function go_down_slice(meshset, section_range, slice)
  new_slice = (slice[1][end]:slice[1][end]+length(slice[1])-1, slice[2])
  view_stack(meshset, section_range, new_slice)
  return new_slice
end

function go_left_slice(meshset, section_range, slice)
  new_slice = (slice[1], slice[2][1]-length(slice[2])-1:slice[2][1])
  view_stack(meshset, section_range, new_slice)
  return new_slice
end

function go_right_slice(meshset, section_range, slice)
  new_slice = (slice[1], slice[2][end]:slice[2][end]+length(slice[2])-1)
  view_stack(meshset, section_range, new_slice)
  return new_slice
end

function load_stack_params(username)
  meshset = load_aligned((1,2,-3,-3), (1,16,-3,-3))
  area = BoundingBox(3500,3500,32000,32000)
  slice = [400, 400]
  return meshset, area, slice, username
end

function get_stack_errors_path(meshset, username)
  firstindex = meshset.meshes[1].index
  lastindex = meshset.meshes[end].index
  fn = string(join(firstindex[1:2], ","), "-", join(lastindex[1:2], ","),
                "_aligned_stack_errors.txt")
  fn = update_filename(fn, username)
  return joinpath(INSPECTION_DIR, fn)
end

function review_stack(username, meshset, area, slice, k, auto=false)
  mov, slice_range = go_to(meshset, area, slice, k)
  errors, escape = mark_stack(mov)
  path = get_stack_errors_path(meshset, username)
  store_stack_errors(path, username, slice_range, errors)
  if auto & !escape
    return review_stack(username, meshset, area, slice, k+1, true)
  end
end

"""
Stores all slice reviews in chronological order - no overwriting
"""
function store_stack_errors(path, username, slice_range, errors)
  ts = Dates.format(now(), "yymmddHHMMSS")
  i, j = slice_range[1][1], slice_range[2][1]
  n, m = slice_range[1][end]-slice_range[1][1], 
                    slice_range[2][end]-slice_range[2][1]
  error_line = [ts, username, i, j, n, m, join(Set(errors), ",")]'
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
      i = round(Int64, (stack_errors[k, 3] - area.i)*s)
      j = round(Int64, (stack_errors[k, 4] - area.j)*s)
      iz = round(Int64, stack_errors[k, 5]*s)
      jz = round(Int64, stack_errors[k, 6]*s)
      errors = stack_errors[k, 7]
      if typeof(errors) != Int64
        errors = readdlm(IOBuffer(stack_errors[k, 7]), ',', Int)
      end
      l = length(errors)
      z[i:i+iz, j:j+jz] = ones(Int64, iz+1, jz+1)*l
    end
  end
  return z
end

function normalize_to_uint8(a)
  assert(minimum(a) == 0)
  mx = maximum(a)
  a /= m
  return convert(Array{UInt8}, round(a*255))
end

function create_stack_colormap()
  return vcat(linspace(RGB(0,0,0), RGB(1,0,0), 100), 
              linspace(RGB(1,0,0), RGB(1,1,0), 100), 
              linspace(RGB(1,1,0), RGB(1,1,1), 55))
end

function view_errors(path, area)
  z = get_stack_errors(path, area)
  a = normalize_to_uint8(z)
  cm = create_stack_colormap()
  b = apply_colormap(a, cm)
  view(b, pixelspacing=[1,1])
end

function go_to(meshset, area, slice, k)
  assert(k != 0)
  n, m = round(Int64, area.h/slice[1]), round(Int64, area.w/slice[2])
  i = (k-1)%n + 1
  j = ceil(Int64, k/n)
  section_range = 1:length(meshset.meshes)
  islice = ((i-1)*slice[1]:i*slice[1]) + area.i
  jslice = ((j-1)*slice[1]:j*slice[1]) + area.j
  return make_stack(meshset, section_range, (islice, jslice)), (islice, jslice)
end

function make_stack(meshset, section_range=(1:2), slice=(1:200,1:200), perm=[1,2,3])
  imgs = []
  for mesh in meshset.meshes[section_range]
    index = (mesh.index[1:2]..., mesh.index[3]-1, mesh.index[4]-1)
    img = reinterpret(UFixed8, get_h5_slice(get_h5_path(index), slice))
    push!(imgs, img)
  end
  img_stack = cat(3, imgs...)
  return Image(permutedims(img_stack, perm), timedim=3)
end

"""
"""
function mark_stack(mov, fps=10)
  e = Condition()

  error_frames = [0]
  escape = false

  imgc, img2 = view(mov, pixelspacing=[1,1])
  start_loop(imgc, img2, fps)
  state = imgc.navigationstate
  set_fps!(state, fps)

  c = canvas(imgc)
  win = Tk.toplevel(c)
  bind(win, "<KP_Enter>", path->count())
  bind(win, "<Escape>", path->end_auto())
  bind(win, "<Destroy>", path->end_count())

  function count()
    t = state.t
    println(t)
    push!(error_frames, t)
  end

  function end_auto()
    escape = !escape
    println("End auto mode: ", escape)
  end

  function end_count()
    notify(e)
    bind(win, "<KP_Enter>", path->path)
    bind(win, "<Escape>", path->path)
    bind(win, "<Destroy>", path->path)
  end
  
  wait(e)
  return error_frames, escape
end

function loop_stack(mov, fps=12)
  imgc, img2 = view(mov, pixelspacing=[1,1])
  start_loop(imgc, img2, fps)

  e = Condition()
  c = canvas(imgc)
  win = Tk.toplevel(c)

  bind(win, "<Destroy>", path->exit_stack())
  bind(win, "<Escape>", path->exit_stack())

  function exit_stack()
    stop_loop(imgc)
    destroy(win)
    bind(win, "<Destroy>", path->path)
    bind(win, "<Escape>", path->path)
    notify(e)
  end

  wait(e)
end

"""
Cycle through sections of the stack stack, with images staged for easier viewing
"""
function view_stack(meshset, section_range=(1:2), slice=(1:200,1:200), perm=[1,2,3])
  imgs = []
  for mesh in meshset.meshes[section_range]
    index = (mesh.index[1:2]..., mesh.index[3]-1, mesh.index[4]-1)
    img = reinterpret(UFixed8, get_h5_slice(get_h5_path(index), slice))
    push!(imgs, img)
  end

  println(slice)
  img_stack = cat(3, imgs..., reverse(imgs[2:end-1])...)  # loop it
  # img_stack = cat(3, imgs...)
  # img_stack = permutedims(cat(3, imgs...), perm)
  # img_stack = Image(img_stack, timedim=3)
  # if perm != [1,2,3]
    # imgc, img2 = view(Image(permutedims(img_stack, [3,2,1]), timedim=3))
    # start_loop(imgc, img2, 10)
    # imgc, img2 = view(Image(permutedims(img_stack, [1,3,2]), timedim=3))
    # start_loop(imgc, img2, 10)
  imgc, img2 = view(Image(permutedims(img_stack, perm), timedim=3), pixelspacing=[1,1])
  start_loop(imgc, img2, 6)
  # end
  # start_loop(imgc, img2, 10)

  e = Condition()
  c = canvas(imgc)
  win = Tk.toplevel(c)

  bind(win, "<Destroy>", path->exit_stack())
  bind(win, "<Escape>", path->exit_stack())

  function exit_stack()
    stop_loop(imgc)
    destroy(win)
    bind(win, "<Destroy>", path->path)
    bind(win, "<Escape>", path->path)
    notify(e)
  end

  wait(e)
  return slice
end