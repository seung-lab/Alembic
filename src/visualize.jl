# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150810

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
function draw_mesh(imgc, img2, nodes, node_dict, color=RGB(1,1,1))   
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

function draw_mesh(img, nodes, node_dict)
    imgc, img2 = view(img)
    return draw_mesh(imgc, img2, nodes, node_dict)
end

"""
Draw vectors on image plot
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
function draw_vectors(imgc, img2, vectors, pt_color=RGB(0,0,1), 
                                                vec_color=RGB(1,0,1), k=100)
    vectors = [vectors[2,:]; 
                vectors[1,:]; 
                (vectors[4,:]-vectors[2,:])*k + vectors[2,:]; 
                (vectors[3,:]-vectors[1,:])*k + vectors[1,:]]
    an_vectors = annotate!(imgc, img2, AnnotationLines(vectors, color=vec_color, 
                                            coord_order="xxyy", linewidth=3))
    an_points = annotate!(imgc, img2, AnnotationPoints(vectors[1:2,:], 
                                                    color=pt_color, shape='*'))
    return an_points, an_vectors
end

function draw_text(imgc, img2, text, point, fontsize=24.0, color=(0,0,0))
  c = canvas(imgc)
  ctx = Cairo.getgc(c)
  Cairo.save(ctx)
  Cairo.move_to(ctx, point)
  Cairo.set_font_size(ctx, fontsize)
  Cairo.set_source_rgb(ctx, color...);
  Cairo.show_text(ctx, text)
  Cairo.stroke(ctx)
  Cairo.restore(ctx)
end

function draw_indices(imgc, img2, points, fontsize=24.0, offset=[-20,-20], color=(0,0,0))
    points = [points[2,:]; points[1,:]] .- offset
    c = canvas(imgc)
    ctx = Cairo.getgc(c)
    Cairo.save(ctx)
    Cairo.set_font_size(ctx, fontsize)
    Cairo.set_source_rgb(ctx, color...);
    for i in 1:size(points,2)
        Cairo.move_to(ctx, points[:,i]...)
        Cairo.show_text(ctx,"$i")
        # an_txt = annotate!(imgc, img2, AnnotationText(points[:,i]..., "$i", 
                                            # fontsize=12, color=RGB(0,0,0)))
    end
    Cairo.stroke(ctx)
    Cairo.restore(ctx)
    return c
end

function draw_indices(imgc, img2, pointsA, pointsB, fontsize, offsetA, offsetB, colorA, colorB)
    pointsA = [pointsA[2,:]; pointsA[1,:]] .- offsetA
    pointsB = [pointsB[2,:]; pointsB[1,:]] .- offsetB
    c = canvas(imgc)
    ctx = Cairo.getgc(c)
    Cairo.save(ctx)
    Cairo.set_font_size(ctx, fontsize)
    Cairo.set_source_rgb(ctx, colorA...)
    for i in 1:size(pointsA,2)
        Cairo.move_to(ctx, pointsA[:,i]...)
        Cairo.show_text(ctx,"$i")
        # an_txt = annotate!(imgc, img2, AnnotationText(points[:,i]..., "$i", 
                                            # fontsize=12, color=RGB(0,0,0)))
    end
    Cairo.set_source_rgb(ctx, colorB...)
    for i in 1:size(pointsB,2)
      Cairo.move_to(ctx, pointsB[:,i]...)
      Cairo.show_text(ctx,"$(i+size(pointsA,2))")
    end
    Cairo.stroke(ctx)
    Cairo.restore(ctx)
    return c
end

function draw_vectors(img, vectors)
    imgc, img2 = view(img)
    return draw_vectors(imgc, img2, vectors)
end

function demo_draw_vectors()
    a = zeros(10,10)
    for i in 1:10
        for j in 1:10
            if i*j > 20
                a[i,j] = 1.0
            end
        end
    end
    vectors = [1.0 1.0 1.0 3.0;
                3.0 1.0 5.0 1.0]'
    a = make_isotropic(a)
    draw_vectors(a, vectors)
end

"""
Display match displacement vectors on images
Args:
  imgc: ImageCanvas object
  img2: ImageSliced2d object 
  nodes: 2xN array of points
Returns:
  imgc: ImageCanvas object
  img2: ImageSliced2d object
  an_points: annotation object for the points
"""
function draw_points(imgc, img2, points, color=RGB(0,0,1), linewidth=1.0, shape='x')
    points = [points[2,:]; points[1,:]]
    an_points = annotate!(imgc, img2, AnnotationPoints(points, 
                                                        color=color, 
                                                        linewidth=linewidth,
                                                        shape=shape))
    return an_points
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

function draw_box(imgc, img2, bb, color=RGB(0,1,0), linewidth=1.0)
    upper_left, lower_right = bb2corners(bb) 
    an_box = annotate!(imgc, img2, AnnotationBox(tuple(upper_left...), 
                                                    tuple(lower_right...),
                                                    color=color, 
                                                    linewidth=linewidth,
                                                    coord_order="yxyx"))
  return an_box
end

function indices2string(indexA, indexB)
  return string(join(indexA[1:2], ","), "-", join(indexB[1:2], ","))
end

function draw_points(img, points)
    imgc, img2 = view(img)
    return draw_points(imgc, img2, points)
end  

function draw_imfuse_meshes(Oc, O2, dst_nodes_A, SR_A, dst_nodes_B, SR_B)
# Incomplete
    SR_A = [0, 0]
    # SR_B = [7184.9, -178.7780] # NEED INTERPOLATION!
    SR_B = [7185, -179]
    SR_C = SR_B - SR_A
    if SR_C[1] > 0
        dst_nodes_B[:, 2] += SR_C[1]
    elseif SR_C[1] < 0
        dst_nodes_A[:, 2] -= SR_C[1]
    end 
    if SR_C[2] > 0
        dst_nodes_B[:, 1] += SR_C[2]
    elseif SR_C[2] < 0
        dst_nodes_A[:, 1] -= SR_C[2]
    end
    draw_mesh(imgc, img2, dst_nodes_B, node_dict_B, RGB(0,0,1))
    draw_mesh(imgc, img2, dst_nodes_A, node_dict_A, RGB(1,0,1))
end

function demo_color_plot()
# Write demo vector color plot to file
    mesh_path = joinpath(BUCKET, "EM_Images", "r4c3_solved.jld")
    write_name = "Tile_r4-c3_S2-W001_sec20_hexplot2.png"
    offset = [29090, 36251]
    write_vector_color_plot(mesh_path, write_name, offset)
end

function demo_mesh()
# Load a mesh and display it on a black background
    mesh_path = joinpath(BUCKET, "test_images", "solvedMesh.jld")
    offset = [21906, 36429]
    v, vt, e = load_mesh(mesh_path, offset)
    incidence = e
    initial_nodes = v
    nodes = xy2yx(initial_nodes')
    node_dict = incidence2dict(incidence)
    sz = round(Int,maximum(nodes[:,1]))+10, round(Int,maximum(nodes[:,2]))+10
    img = zeros(Bool, sz...)

    println("Original images and shapes")
    draw_mesh(make_isotropic(img), nodes, node_dict)
end

function draw_reference_vector(imgc, img2, ratio, bb, use_index=true)
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

# function invisible_view(img)
#   props = Dict{Symbol,Any}()
#   img2 = ImageView.ImageSlice2d(img, props)
#   imgc = ImageView.ImageCanvas(ImageView.cairo_format(img), props)
#   w = ImageView.width(img2)
#   h = ImageView.height(img2)
#   ww, wh = ImageView.rendersize(w, h, imgc.aspect_x_per_y)
#   win = ImageView.Toplevel("ImageView_invisible", ww, wh, false)
#   framec = OS_NAME == :Darwin ? ImageView.Frame(win, padding = 30) : ImageView.Frame(win)
#   ImageView.grid(framec, 1, 1, sticky="nsew")
#   ImageView.grid_rowconfigure(win, 1, weight=1) # scale this cell when the window resizes
#   ImageView.grid_columnconfigure(win, 1, weight=1)
#   c = ImageView.Canvas(framec, ww, wh)
#   imgc.c = c
#   ImageView.grid(c, 1, 1, sticky="nsew")        # fill the edges of its cell on all 4 sides
#   ImageView.grid_rowconfigure(framec, 1, weight=1) # scale this cell when the window resizes
#   ImageView.grid_columnconfigure(framec, 1, weight=1)
#   lastrow = 1
#   fnotify = ImageView.Frame(win)
#   ImageView.grid(fnotify, lastrow+=1, 1, sticky="ew")
#   xypos = ImageView.Label(fnotify)
#   imgc.handles[:pointerlabel] = xypos
#   ImageView.grid(xypos, 1, 1, sticky="ne")
#   ImageView.set_visible(win, false)
#   ImageView.allocate_surface!(imgc, w, h)
#   ImageView.rerender(imgc, img2)
#   ImageView.resize(imgc, img2)
#   Tk.configure(imgc.c)
# end