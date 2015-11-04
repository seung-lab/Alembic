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
  imgc: ImageCanvas object
  img2: ImageSliced2d object
  an_points: annotation object for the points
"""
function show_points(imgc, img2, points, color=RGB(0,0,1), linewidth=1.0, shape='x')
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

function show_box(imgc, img2, bb, color=RGB(0,1,0), linewidth=1.0)
    upper_left, lower_right = bb2corners(bb) 
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

function uview(img::Array{UInt8,2})
  return ImageView.view(convert(Array{Ufixed8}, img))
end