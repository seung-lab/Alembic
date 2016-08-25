function cpselect(moving, fixed)
    moving = convert_uint8(moving)
    fixed = convert_uint8(fixed)
    cgrid = canvasgrid(1,2; pad=10)
    opts = Dict(:pixelspacing => [1,1])
    moving_imgc, moving_img2 = ImageView.view(cgrid[1,1], moving; opts...)
    fixed_imgc, fixed_img2 = ImageView.view(cgrid[1,2], fixed; opts...)
    mpts = bind_select_and_remove(moving_imgc, moving_img2, RGB(0,1,0))
    fpts = bind_select_and_remove(fixed_imgc, fixed_img2, RGB(1,0,0))
    return mpts, fpts
end

function select_points(index::Index; thumb=false)
    moving_index = index
    assert(moving_index != NO_INDEX)
    fixed_index = get_preceding(index)
    assert(fixed_index != NO_INDEX)
    return select_points(moving_index, fixed_index, thumb=thumb)
end

"""
Select corresponding points between two images

* Use right-click to add a point in one image
* Use control+right-click to remove the last point added in one image
"""
function select_points(moving_index::Index, fixed_index::Index; thumb=false)
    scale = 1.0
    if thumb
        scale = h5read(get_path("thumbnail", moving_index), "scale")
    end
    moving_img = thumb ? load("thumbnail", moving_index) : load(moving_index)
    fixed_img = thumb ? load("thumbnail", fixed_index) : load(fixed_index)
    moving_img = convert_uint8(moving_img)
    fixed_img = convert_uint8(fixed_img)
    cgrid = canvasgrid(1,2; pad=10)
    opts = Dict(:pixelspacing => [1,1])
    moving_imgc, moving_img2 = ImageView.view(cgrid[1,1], moving_img; opts...)
    fixed_imgc, fixed_img2 = ImageView.view(cgrid[1,2], fixed_img; opts...)
    # ms = MeshSet()
    mpts = bind_select_and_remove(moving_imgc, moving_img2, RGB(0,1,0))
    fpts = bind_select_and_remove(fixed_imgc, fixed_img2, RGB(1,0,0))
    bind_save(moving_imgc, moving_img2, moving_index, mpts, fpts, scale, thumb)
    bind_save(fixed_imgc, fixed_img2, moving_index, mpts, fpts, scale, thumb)
    bind_go(moving_imgc, moving_img2, moving_index, thumb)
    bind_go(fixed_imgc, fixed_img2, moving_index, thumb)
    bind_exit(moving_imgc, moving_img2)
    bind_exit(fixed_imgc, fixed_img2)
end

function bind_select_and_remove(imgc, img2, color=RGB(1,0,0), pts=[], ann=[])
  c = canvas(imgc)
  win = Tk.toplevel(c)
  bind(c, "<Button-3>", (c, x, y)->select_point(imgc, 
                                                img2,
                                                parse(Int64, x), 
                                                parse(Int64, y), 
                                                pts, 
                                                ann,
                                                color))
  bind(c, "<Control-Button-3>", path->remove_last_point(imgc, 
                                                img2,
                                                pts,
                                                ann))
  return pts
end

function bind_save(imgc, img2, index, mpts, fpts, scale, thumb)
    c = canvas(imgc)
    win = Tk.toplevel(c)
    bind(win, "s", path->save_correspondences(win, index, mpts, fpts, scale, thumb))
end

function bind_go(imgc, img2, index, thumb)
    c = canvas(imgc)
    win = Tk.toplevel(c)
    next = get_succeeding(index)
    previous = get_preceding(index)
    bind(win, "<Control-Right>", path->go_to(win, next, thumb))
    bind(win, "<Control-Left>", path->go_to(win, previous, thumb))
end

function bind_exit(imgc, img2)
    c = canvas(imgc)
    win = Tk.toplevel(c)    
    bind(win, "<Escape>", path->exit_cpselect(win))
end

function select_point(imgc, img2, x, y, pts=[], ann=[], color=RGB(1,0,0))
    r = getgc(imgc.c)
    xu, yu = ImageView.device_to_user(r, x, y)
    println((xu, yu))
    push!(pts, [yu, xu])
    n = length(pts)
    push!(ann, annotate!(imgc, img2, AnnotationText(xu, yu, "$n", 
                                                color=color, fontsize=64)))
end

function remove_last_point(imgc, img2, pts, ann)
    if length(pts) > 0 && length(ann) > 0
        pop!(pts)
        a = pop!(ann)
        ImageView.delete!(imgc, a)
    end
end

# function recompute_tform(firstindex::Index, lastindex::Index)
#     for index in get_index_range(firstindex, lastindex)
#         pts = load("import", index)
#         n = size(pts, 1)
#         src_pts = pts[:,1:2]
#         dst_pts = pts[:,3:4]
#         println("Resaving transform for $index")
#         tform_path = get_path("stats", index)
#         tform = compute_tform(src_pts, dst_pts)
#         writedlm(tform_path, tform)
#     end
# end

function save_correspondences(win, index, mpts, fpts, scale, thumb)
    println("Saving correspondences")
    correspondences_path = get_path("correspondence", index)
    n = min(length(mpts), length(fpts))
    mpts = mpts[1:n]
    fpts = fpts[1:n]
    src_pts = hcat(mpts...)'/scale
    dst_pts = hcat(fpts...)'/scale
    pts = hcat(src_pts, dst_pts)
    writedlm(correspondences_path, pts)
    tform = compute_tform(src_pts, dst_pts)
    println("Saving transform")
    # tform_path = get_path("relative_transform", index)
    # writedlm(tform_path, tform)
    update_registry(index, tform)
    next = get_succeeding(index)
    go_to(win, next, thumb)
end

function compute_tform(src_pts, dst_pts, lambda=0)
    n = size(src_pts, 1)
    src_pts = hcat(src_pts, ones(n))
    dst_pts = hcat(dst_pts, ones(n))
    return lambda*calculate_affine(src_pts, dst_pts) + (1-lambda)*calculate_rigid(src_pts, dst_pts)
end

function go_to(win, index, thumb)
    println("Go to $index")
    exit_cpselect(win)
    select_points(index, thumb=thumb)
end

function exit_cpselect(win)
    destroy(win)
end

function convert_uint8(img)
    if typeof(img) == Array{UInt8, 2}
        img = img / 255
    end
    return img
end