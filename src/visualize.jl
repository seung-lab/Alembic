"""
Turn 2D array into Image object with isotropic pixel spacing for display

Args:

* img: 2D array

Returns:

* img: Image object for isotropic ImageView display
"""
function make_isotropic(img)
    img = Image(img)
    img["pixelspacing"] = [1, 1]
    return img
end

function view_isotropic(img)
    return view(img, pixelspacing=[1,1])
end

"""
Display thumbnail image with vectors and match indices overlayed
"""
function view_matches(img, vectors, factor)
  imgc, img2 = view_isotropic(img)
  a, b = draw_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1), factor)
  c = draw_indices(imgc, img2, vectors[1:2,:])
  return imgc, img2
end

"""
Display mesh on image
Args:
  imgc: ImageCanvas object
  img2: ImageSliced2d object  
  nodes: Nx2 array of mesh node positions
  node_dict: dictionary of edges, indexed by node, containing connected nodes
Returns:
  imgc: ImageCanvas object
  img2: ImageSliced2d object 
"""
function show_mesh(imgc, img2, nodes, node_dict, color=RGB(1,1,1))   
    lines = Array(Float64, 4, 0)
    for k in sort(collect(keys(node_dict)))
        for v in node_dict[k]
            a = reverse(vec(nodes[k,:]))
            b = reverse(vec(nodes[v,:]))
            lines = hcat(lines, vcat(a, b))
        end
    end
    an_lines = annotate!(imgc, img2, AnnotationLines(lines, color=color, 
                                                            coord_order="yyxx"))
    return an_lines
end

function show_mesh(img, nodes, node_dict)
    imgc, img2 = view(img)
    return draw_mesh(imgc, img2, nodes, node_dict)
end

"""
Show vectors on image plot
Args:
  imgc: ImageCanvas object
  img2: ImageSliced2d object 
  vectors: 4xN array of vector start and end points
  pt_color: optional color argument for points
  vec_color: optional color argument for vectors
Returns:
  imgc: ImageCanvas object
  img2: ImageSliced2d object
  an_points: annotation object for the points
  an_vectors: annotation object for the vectors
"""
function show_vectors(imgc, img2, vectors, pt_color=RGB(0,0,1), 
                                                vec_color=RGB(1,0,1))
  if length(vectors) > 0
    an_vectors = annotate!(imgc, img2, AnnotationLines(vectors, color=vec_color, 
                                            coord_order="xxyy", linewidth=3))
    an_points = annotate!(imgc, img2, AnnotationPoints(vectors[1:2,:], 
                                                    color=pt_color, shape='*'))
  end
end

function show_vectors(img, vectors)
    imgc, img2 = view(img, pixelspacing=[1,1])
    return show_vectors(imgc, img2, vectors)
end

"""
Display match displacement vectors on images
Args:
  imgc: ImageCanvas object
  img2: ImageSliced2d object 
  nodes: 2xN array of points
Returns:
  an_points: annotation object for the points
"""
function show_points(imgc, img2, points; color=RGB(0,0,1), linewidth=2.0, size=10.0, shape='o', t=NaN)
  an_points = nothing
  if length(points) > 0
    an_points = annotate!(imgc, img2, AnnotationPoints(points, 
                                                        color=color, 
                                                        linewidth=linewidth,
                                                        size=size,
                                                        shape=shape, t=t, scale=true))
  end
  return an_points
end 

function show_lines(imgc, img2, lines; color=RGB(0,0,1), linewidth=3.0, t=NaN)
  an_vectors = nothing
  if length(lines) > 0
    an_vectors = annotate!(imgc, img2, AnnotationLines(lines, color=color, t=t,
                                            coord_order="xxyy", linewidth=linewidth))
  end
  return an_vectors
end 

function show_colored_lines(imgc, img2, lines, colors; linewidth=3.0, t=NaN)
  an_vectors = nothing
  if length(lines) > 0
    an_vectors = annotate!(imgc, img2, AnnotationColoredLines(lines, colors, t=t,
                                            coord_order="xxyy", linewidth=linewidth))
  end
  return an_vectors
end 

function show_text(imgc, img2, str, x, y; color=RGB(0,0,1), fontsize=30, t=NaN)
  return annotate!(imgc, img2, AnnotationText(x, y, str, t=t, color=color, fontsize=fontsize, halign="left"))
end

function imwrite_box(img, point, radius, path, color=RGB(0,1,0), linewidth=1.0)
  upper_left = point - [radius, radius]
  lower_right = point + [radius, radius]
  imgc, img2 = view(img, pixelspacing=[1,1])
  annotate!(imgc, img2, AnnotationBox(tuple(upper_left...), 
                                        tuple(lower_right...),
                                        color=color, 
                                        linewidth=linewidth,
                                        coord_order="yxyx"))
  write_canvas(imgc, path)
  close_image(imgc)
end    

function write_canvas(imgc, path)
  Cairo.write_to_png(imgc.c.back, path)
end

function close_image(imgc)
  destroy(toplevel(imgc))
end

function set_canvas_size(imgc, w, h)
  ImageView.set_size(toplevel(imgc), floor(Int64, w), floor(Int64, h))
end

function create_image(bb)
    z = ones(bb.h+1, bb.w+1)
    return view_isotropic(z), [bb.i, bb.j]
end

function show_box(imgc, img2, upper_left, lower_right, color=RGB(0,1,0), linewidth=1.0)
    an_box = annotate!(imgc, img2, AnnotationBox(tuple(upper_left...), 
                                                    tuple(lower_right...),
                                                    color=color, 
                                                    linewidth=linewidth,
                                                    coord_order="yxyx"))
  return an_box
end

function show_points(img, points)
    imgc, img2 = view(img)
    return draw_points(imgc, img2, points)
end  

function show_reference_vector(imgc, img2, ratio, bb, use_index=true)
  reference = 1
  margin = 0.02
  vectors = [bb.h*(1-margin);
              bb.w*(1-margin) - reference*ratio;
              bb.h*(1-margin);
              bb.w*(1-margin)]
  draw_vectors(imgc, img2, vectors, RGB(1,1,0), RGB(0,1,1), 1)
  if use_index
    draw_indices(imgc, img2, vectors[1:2,:], 18.0, [-reference*ratio/2, 20], (0.5,0.5,0.5))
  end
  return imgc, img2
end

"""
Write thumbnail image, with whatever drawings included
"""
function write_imageview(path, imgc, img2)
  println("Writing ", path)
  write_canvas(imgc, path)
  close_image(imgc)
end

function xcview(xc)
  xc_image = xcorr2Image(xc)
  hot = create_hot_colormap()
  xc_color = apply_colormap(xc_image, hot)
  ImageView.view(xc_color') #, pixelspacing=[1,1])
end

function uview(img::Array{UInt8,2})
  return ImageView.view(reinterpret(UFixed8, img), pixelspacing=[1,1])
end

function override_xy_label(imgc, img2, offset, scale=1.0)
  c = canvas(imgc)
  win = Tk.toplevel(c)
  fnotify = ImageView.Frame(win)
  lastrow = 2
  ImageView.grid(fnotify, lastrow+=1, 1, sticky="ew")
  xypos = ImageView.Label(fnotify)
  imgc.handles[:pointerlabel] = xypos
  ImageView.grid(xypos, 1, 1, sticky="ne")
  ImageView.set_visible(win, true)
  c.mouse.motion = (path,x,y)-> updatexylabel(xypos, imgc, img2, x, y, offset..., scale)
end

function updatexylabel(xypos, imgc, img2, x, y, x_off, y_off, scale=1.0)
  w = width(imgc.c)
  xu, yu = ImageView.device_to_user(Tk.getgc(imgc.c), x, y)
  # Image-coordinates
  xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
  if ImageView.isinside(imgc.canvasbb, x, y)
    val = img2[xi,yi]
    xo, yo = xi + x_off, yi + y_off
    xo, yo = round(Int64, xo / scale), round(Int64, yo / scale)
    str = "$yo, $xo ($yi, $xi): $val"
    if length(str)*10>w
      ImageView.set_value(xypos, "$yo, $xo ($yi, $xi)")
    else
      ImageView.set_value(xypos, str)
    end
  else
    ImageView.set_value(xypos, "($yi, $xi)")
  end
end
