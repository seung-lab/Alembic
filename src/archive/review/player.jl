"""
Ripped from ImageView > navigation.jl
"""
function stop_playing!(state::ImageView.Navigation.NavigationState)
    if state.timer != nothing
        close(state.timer)
        state.timer = nothing
    end
end

# function update_player_annotations(imgc, img2, annotations, t)
#   if length(annotations) > 0
#     display_list = annotations["display"]
#     for k in display_list
#       data = annotations["annotations"][k]["data"][:, :, t]
#       if annotations["annotations"][k]["hash"] == nothing
#         func = annotations["annotations"][k]["function"]
#         handle = eval(func)(imgc, img2, data)
#         annotations["annotations"][k]["hash"] = handle.hash
#         annotations["annotations"][k]["symbol"] = fieldnames(handle.ann.data)[1]
#       else
#         h = annotations["annotations"][k]["hash"]
#         sym = annotations["annotations"][k]["symbol"]
#         imgc.annotations[h].data.(sym) = data
#       end
#     end
#     ImageView.redraw(imgc)
#   end
# end

"""
Ripped from ImageView > navigation.jl
"""
function updatet(ctrls, state, t_display)
  Tk.set_value(ctrls.editt, t_display)
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
function incrementt(inc, ctrls, state, showframe, annotations)
    state.t += inc
    t_display = string(state.t)
    if length(annotations) > 0
      index = annotations["indices"][state.t]
      t_display = string(join(index[1:2], ","))
    end
    updatet(ctrls, state, t_display)
    # update_player_annotations(imgc, img2, annotations, state.t)
    showframe(state)
end

"""
Ripped from ImageView > navigation.jl
"""
function playt(inc, ctrls, state, showframe, annotations)
    if !(state.fps > 0)
        error("Frame rate is not positive")
    end
    stop_playing!(state)
    dt = 1/state.fps
    state.timer = Timer(timer -> stept(inc, ctrls, state, showframe, annotations), dt, dt)
end

"""
Ripped from ImageView > navigation.jl
"""
function stept(inc, ctrls, state, showframe, annotations)
    if 1 <= state.t+inc <= state.tmax
        incrementt(inc, ctrls, state, showframe, annotations)
    else
        stop_playing!(state)
    end
end

"""
Endlessly repeat forward loop (continuous stept)
"""
function loopt(ctrls, state, showframe, annotations)
  inc = 1
  if state.t+inc < 1 || state.tmax < state.t+inc
      state.t = 0
  end
  incrementt(inc, ctrls, state, showframe, annotations)
end

"""
Building on ImageView > navigation.jl
"""
function set_fps!(state, fps)
  state.fps = fps
end

"""
Create loop of the image
"""
function start_loop(imgc, img2, showframe, fps=6, annotations=Dict())
  state = imgc.navigationstate
  ctrls = imgc.navigationctrls
  set_fps!(state, fps)

  if !(state.fps > 0)
      error("Frame rate is not positive")
  end
  stop_playing!(state)
  dt = 1/state.fps
  state.timer = Timer(timer -> loopt(ctrls, state, showframe, annotations), dt, dt)
end

"""
Higher level call for ImageView outputs
"""
function stop_loop(imgc)
  stop_playing!(imgc.navigationstate)
end

"""
Cycle through sections of the stack stack, with images staged for easier viewing
"""
function view_stack(firstindex::Index, lastindex::Index, center, radius;
                                        include_reverse=false, perm=[1,2,3], scale=1.0, thumb=false)
  slice = make_slice(center, radius)
  return view_stack(firstindex, lastindex, slice, include_reverse=include_reverse, perm=perm, scale=scale, thumb=thumb)
end

function view_stack(firstindex::Index, lastindex::Index, slice=(1:200,1:200); include_reverse=false, perm=[1,2,3], scale=1.0, thumb=false)
  stack = make_stack(firstindex, lastindex, slice, scale=scale, thumb=thumb)
  offset = ImageRegistration.get_offset(slice_to_bb(slice))
  indices = get_index_range(firstindex, lastindex)
  annotations = Dict("indices" => [indices..., reverse(indices)...], "slice"=>slice, "scale"=>scale)
  imgc, img2 = view_stack(stack, offset=offset, scale=scale, annotations=annotations, include_reverse=include_reverse, perm=perm)
  return stack, annotations, imgc, img2
end

function view_stack(stack; offset=[0,0], scale=1.0, annotations=Dict(), include_reverse=false, perm=[1,2,3])
  img_stack = permutedims(reinterpret(UFixed8, stack), perm)
  if include_reverse
    img_stack = cat(3, img_stack, img_stack[:,:,end:-1:1])
  end
  imgc, img2 = ImageView.view(Images.Image(img_stack, timedim=3), pixelspacing=DATASET_RESOLUTION[perm[2:-1:1]])
  override_xy_label(imgc, img2, offset, 1/scale)

  c = canvas(imgc)
  win = Tk.toplevel(c)
  state = imgc.navigationstate
  ctrls = imgc.navigationctrls
  showframe = state -> ImageView.reslice(imgc, img2, state)
  incrementt(0, ctrls, state, showframe, annotations)
  bind(win, "<Up>", path->adjust_fps(imgc, img2, showframe, 1, annotations))
  bind(win, "<Down>", path->adjust_fps(imgc, img2, showframe, -1, annotations))
  bind(win, "<Escape>", path->exit_stack(imgc, img2))
  bind(win, "<Destroy>", path->exit_stack(imgc, img2))
  bind(win, "<Alt-Right>", path->stept(1,ctrls,state,showframe, annotations))
  bind(win, "<Alt-Left>", path->stept(-1,ctrls,state,showframe, annotations))
  bind(win, "<Alt-Shift-Right>", path->playt(1,ctrls,state,showframe, annotations))
  bind(win, "<Alt-Shift-Left>", path->playt(-1,ctrls,state,showframe, annotations))
  bind(c, "<Control-Button-3>", (c, x, y)->get_line_index(imgc, img2, 
                                                        parse(Int, x), 
                                                        parse(Int, y), 
                                                        annotations))
  # bind(c, "<Button-3>", (c, x, y)->inspect_match(imgc, img2, 
  #                                                       parse(Int, x), 
  #                                                       parse(Int, y), 
  #                                                       annotations))
  return imgc, img2
end

# function inspect_match(imgc, img2, x, y, matches, prox=0.0125)
#   # prox: 0.0125 = 100/8000
#   t = imgc.navigationstate.t
#   index = annotations["indices"][t]
#   ann = annotations["ann"][index]["match"]
#   vectors = []
#   for (k, v) in ann
#     push!(vectors, v)
#   end
#   vectors = hcat(vectors...)

#   xu, yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
#   xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
#   limit = (img2.zoombb.xmax - img2.zoombb.xmin) * prox 
#   annidx = find_idx_of_nearest_pt(vectors[1:2,:], [xi, yi], limit)
#   if annidx > 0
#     idx = indices[annidx]
#     ptA = vectors[1:2,idx] # - params["src_offset"]
#     ptB = vectors[3:4,idx] # - params["dst_offset"]
#     println(idx, ": ", ptA, ", ", ptB)
#     bm_win = view_blockmatch(matches, idx)
#     detect_blockmatch_removal(imgc, img2, bm_win, matches, idx, params)
#   end


function make_vectors(meshset::MeshSet, match_ind::Int, offset, scale=10)
  src_points, dst_points, filtered_inds = get_correspondences(meshset, match_ind; globalized = true, use_post = true)
  return offset_points(src_points, dst_points, offset), filtered_inds
end

function mask_and_combine_points(src, dst, mask)
  src_masked, dst_masked = src[mask], dst[mask]
  return [[i...] for i in zip(src_masked, dst_masked)]
end

function filter_contained_points(vectors, bb)
  return filter(i->point_is_contained(bb, i[1:2]), vectors)
end

function filter_contained_lines(vectors, bb)
  return filter(i->line_is_contained(bb, i), vectors)
end

function compile_match_annotations(meshset::MeshSet, match::Match, bb::ImageRegistration.BoundingBox)
  match_ind = find_match_index(meshset, match)
  offset = ImageRegistration.get_offset(bb)
  local_bb = translate_bb(bb, -offset)
  (src, dst), accepted_inds = make_vectors(meshset::MeshSet, match_ind::Int, offset)
  if length(src) > 0
    rejected_inds = collect(setdiff(1:length(src), accepted_inds))
    accepted = mask_and_combine_points(src, dst, accepted_inds)
    rejected = mask_and_combine_points(src, dst, rejected_inds)
    accepted_vectors = transpose_vectors(hcat(filter_contained_points(accepted, local_bb)...))
    rejected_vectors = transpose_vectors(hcat(filter_contained_points(accepted, local_bb)...))
    return accepted_vectors, rejected_vectors
  else
    return [], []
  end
end

function compile_mesh_annotations(mesh::Mesh, bb::ImageRegistration.BoundingBox, use_prealigned::Bool)
  offset = ImageRegistration.get_offset(bb)
  local_bb = translate_bb(bb, -offset)
  endpoints_a, endpoints_b = get_edge_endpoints(mesh; use_post = !use_prealigned)
  edges = [[endpoints_a[i] - offset, endpoints_b[i] - offset] for i in 1:length(endpoints_a)]
  edge_indices = get_edge_indices(mesh)
  edges_removed_indices = get_removed_edge_indices(mesh)
  edges_fixed_indices = get_fixed_edge_indices(mesh)
  mask = Bool[map(line_is_contained, repeated(local_bb), edges)...]
  edges_to_display = hcat(edges[mask]...)
  edges_to_display = transpose_vectors(edges_to_display)
  lengths_pre = get_edge_lengths(mesh; use_post = false)
  lengths_post = get_edge_lengths(mesh; use_post = true)
  strain = (lengths_post - lengths_pre) ./ lengths_pre
  strain_to_display = strain[mask]
  edge_indices_included = edge_indices[mask]
  edges_removed_included = intersect(edge_indices_included, edges_removed_indices)
  edges_fixed_included = intersect(edge_indices_included, edges_fixed_indices)
  return edges_to_display, strain_to_display, edge_indices_included, edges_removed_included, edges_fixed_included
end

function compile_annotations(parent_name, firstindex::Index, lastindex::Index, slice)
  ms = compile_partial_meshset(parent_name, firstindex, lastindex)
  return ms, compile_annotations(ms, firstindex, lastindex, slice)
end

function compile_annotations(meshset::MeshSet, center::Tuple{Int,Int}, radius=1024; scale=1.0, use_prealigned=false)
  x, y = center
  slice = (x-radius):(x+radius), (y-radius):(y+radius)
  return compile_annotations(meshset, slice, scale)
end

function compile_annotations(meshset::MeshSet, 
                              slice::Tuple{UnitRange{Int64},UnitRange{Int64}};
                              scale::Float64=1.0, use_prealigned::Bool=false)
  indices = map(get_index, meshset.meshes)
  firstindex, lastindex = minimum(indices), maximum(indices)
  bb = slice_to_bb(slice)
  annotations = Dict("ann" => Dict(), 
                 "meshset" => meshset, 
                 "slice" => slice,
                 "bb" => bb,
                 "scale" => scale,
                 "vector_scale" => 10,
                 "use_prealigned" => use_prealigned)
  for index in indices
    annotations["ann"][index] = Dict()
    annotations["ann"][index]["match"] = Dict()
    annotations["ann"][index]["match_names"] = Dict()
    match_inds = find_match_indices(meshset, index)
    for i in match_inds
      src_index, dst_index = get_src_and_dst_indices(meshset.matches[i])
      if (firstindex <= src_index <= lastindex) && 
                  (firstindex <= dst_index <= lastindex)
        vectors, _ = compile_match_annotations(meshset, meshset.matches[i], bb)
        annotations["ann"][index]["match"][i] = vectors
        annotations["ann"][index]["match_names"][i] = get_name(meshset.matches[i])
      end
    end
    mesh = get_mesh(meshset, index)
    edges, strain, edge_indices, edges_removed, edges_fixed = compile_mesh_annotations(mesh, bb, use_prealigned)
    edges_removed_mask = Bool[i in edges_removed for i in edge_indices]
    edges_fixed_mask = Bool[i in edges_fixed for i in edge_indices]
    annotations["ann"][index]["mesh"] = Dict("edges" => edges, 
                                              "strain" => strain,
                                              "indices" => edge_indices,
                                              "removed" => edges_removed_mask,
                                              "fixed" => edges_fixed_mask)
  end
  if !use_prealigned
    indices = get_index_range(aligned(firstindex), aligned(lastindex))
  end
  annotations["indices"] = [indices, reverse(indices)]
  return annotations
end

function display_roi(imgc, img2, annotations)
  indices = annotations["indices"]
  bb = slice_to_bb(annotations["slice"])
  offset = ImageRegistration.get_offset(bb)
  for (i, index) in enumerate(indices[1:Int(length(indices)/2)])
    println("Display annotations for $index")
    pts = load("import", index)
    start_pts = hcat([reverse(pts[i,:][:] - offset) for i in 1:size(pts,1)]...)
    end_pts = start_pts[:,2:end]
    # show_points(imgc, img2, pts, shape='.', color=RGB(0,1,0), t=i)
    lines = vcat(start_pts[:,1:end-1], end_pts)
    show_lines(imgc, img2, lines, linewidth=2, color=RGB(0,1,0), t=i)
  end
end

function display_annotations(imgc, img2, annotations; include_reverse=false)
  colors = [RGB(0,1,0), RGB(1,0,1), RGB(0,0.8,0), RGB(0.8,0,0.8)]
  bwr = create_bwr_colormap()
  indices = annotations["indices"]
  vector_scale = annotations["vector_scale"]
  scale = annotations["scale"]
  if !include_reverse
    indices = indices[1:end/2]
  end
  all_strains = [annotations["ann"][i]["mesh"]["strain"] for i in indices]
  all_strains = vcat(all_strains...)
  min_s = minimum(all_strains)
  max_s = maximum(all_strains)
  print("\n")
  min_s = min_s == 0 ? 1 : min_s
  for (i, index) in enumerate(indices)
    println("Display annotations for $index")
    data = annotations["ann"][index]["mesh"]["edges"]*scale
    # strain = map(i->RGB(i,i,i), annotations["ann"][index]["mesh"]["strain"])
    s = copy(annotations["ann"][index]["mesh"]["strain"])
    st = annotations["ann"][index]["mesh"]["strain"]
    min_s = minimum(st)
    max_s = maximum(st)
    println("Min/max mesh strain: $min_s / $max_s")
    if min_s != 0
      s[st.<=0] = round(UInt8, -st[st.<=0]./min_s*127+128)
    end
    if max_s != 0
      s[st.>0] = round(UInt8, st[st.>0]./max_s*127+128)
    end
    strain = apply_colormap(s, bwr)
    removed_mask = annotations["ann"][index]["mesh"]["removed"]
    fixed_mask = annotations["ann"][index]["mesh"]["fixed"]
    strain[removed_mask] = [RGB(0,0,0) for i in 1:sum(removed_mask)]
    strain[fixed_mask] = [RGB(1,0,1) for i in 1:sum(fixed_mask)]
    if size(data, 2) > 1
      show_colored_lines(imgc, img2, data, strain, t=i)
    end
    for (k, match_ind) in enumerate(keys(annotations["ann"][index]["match"]))
      if length(annotations["ann"][index]["match"][match_ind]) > 0
        data = change_vector_lengths(annotations["ann"][index]["match"][match_ind], vector_scale)*scale
        name = annotations["ann"][index]["match_names"][match_ind]
        if size(data, 2) > 1
          show_points(imgc, img2, data[1:2, :], shape='o', color=colors[k], t=i)
          show_lines(imgc, img2, data, color=colors[k], t=i)
          fontsize = 40
          show_text(imgc, img2, name, 2*fontsize, 2*fontsize*k, fontsize=fontsize, color=colors[k], t=i)
        end
      end
    end
  end
end

# function view_annotated_stack(parentfirst::Index, parentlast::Index, firstindex::Index, lastindex::Index, slice=(1:200,1:200);
#                                         include_reverse=false, perm=[1,2,3])
#   annotations = compile_annotations(get_name(parentfirst, parentlast), firstindex, lastindex, slice, include_reverse)
#   stack = make_stack(firstindex, lastindex, slice);
#   imgc, img2 = view_annotated_stack(stack, annotations)
#   return stack, imgc, img2, annotations
# end

function make_stack(annotations)
  return make_stack(minimum(annotations["indices"]), maximum(annotations["indices"]), annotations["slice"], scale=annotations["scale"])
end

function view_stack(annotations::Dict; include_reverse=false)
  offset = ImageRegistration.get_offset(annotations["bb"])
  scale = annotations["scale"]
  stack = make_stack(annotations)
  imgc, img2 = view_stack(stack, offset=offset, scale=scale, annotations=annotations, include_reverse=include_reverse)
  return stack, imgc, img2
end

function view_annotated_stack(annotations; include_reverse=false)
  stack = make_stack(annotations)
  imgc, img2 = view_annotated_stack(stack, annotations; include_reverse=include_reverse)
  return stack, imgc, img2
end

function view_annotated_stack(stack, annotations; include_reverse=false)
  offset = ImageRegistration.get_offset(annotations["bb"])
  scale = annotations["scale"]
  imgc, img2 = view_stack(stack, offset=offset, scale=scale, annotations=annotations, include_reverse=include_reverse)
  display_annotations(imgc, img2, annotations, include_reverse=include_reverse)
  return imgc, img2
end

function adjust_fps(imgc, img2, showframe, d, annotations)
  state = imgc.navigationstate
  fps = state.fps
  fps += d
  if 1 <= fps <= 50
    stop_loop(imgc)
    set_fps!(state, fps)
    start_loop(imgc, img2, showframe, fps, annotations)
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
  bind(win, "<Alt-Right>", path->path)
  bind(win, "<Alt-Left>", path->path)
  bind(win, "<Alt-Shift-Right>", path->path)
  bind(win, "<Alt-Shift-Left>", path->path)
  bind(c, "<Control-Button-3>", path->path)
  bind(c, "<Button-3>", path->path)
  destroy(win)
end

function get_line_index(imgc, img2, x, y, annotations, prox=0.0125)
  # prox: 0.0125 = 100/8000
  t = imgc.navigationstate.t
  index = annotations["indices"][t]
  edges = annotations["ann"][index]["mesh"]["edges"]

  xu, yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
  xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
  limit = (img2.zoombb.xmax - img2.zoombb.xmin) * prox 
  idx = find_idx_of_nearest_line(edges, [xi, yi], limit)
  if idx != 0
    ms = annotations["meshset"]
    mesh = get_mesh(ms, prealigned(index));
    remove_edge!(mesh, annotations["ann"][index]["mesh"]["indices"][idx])
    println(annotations["ann"][index]["mesh"]["indices"][idx])
  end
end

"""
Find index of location in array with point closest to the given point (if there
exists such a unique point within a certain pixel limit).

Args:

* pts: 2xN array of coordinates
* pt: 2-element coordinate array of point in interest
* limit: number of pixels that nearest must be within

Returns:

* index of the nearest point in the pts array

  idx = find_idx_of_nearest_pt(pts, pt, limit)
"""
function find_idx_of_nearest_line(lines, pt, limit)
  pts = lines[1:2, :]
  idx = 0
  min_d = Inf
  for i in 1:size(pts, 2)
    closest_pt_vector = [Inf, Inf]
    if lines[3:4, i] == lines[1:2, i]
      closest_pt_vector = lines[1:2, i]
    else
      pt_vector = pt - lines[1:2, i]
      line_vector = lines[3:4, i] - lines[1:2, i]
      line_magnitude = norm(line_vector)
      closest_magnitude = pt_vector'*line_vector / (line_magnitude)^2
      t = max(0, min(1, closest_magnitude[1]))
      closest_pt_vector = lines[1:2, i] + t*line_vector
    end
    displacement_vector = pt - closest_pt_vector
    d = norm(displacement_vector)
    if d < min_d
      min_d = d
      idx = i
    end
  end
  if idx != 0
    line = lines[:, idx]
    m = @sprintf("%0.2f", min_d)
    l = round(line, 1)
    println("$pt is $m px from line $l @ $idx")
 end
 return idx
end

	function compare_alignments(ir, jr)
	so = make_stack(prealigned(1,31), prealigned(1,42), (ir,jr))
	sn = make_stack(aligned(1,31), aligned(1,42), (ir,jr))
	blackbar = fill(UInt8(0), 50, length(jr), 12)
	view_stack(cat(1, so, blackbar, sn))
	end
