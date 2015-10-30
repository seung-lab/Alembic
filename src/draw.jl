function create_drawing(img::Array{UInt32,2})
  return Cairo.CairoRGBSurface(img)
end

function create_contex(surface::Cairo.CairoSurface)
  return Cairo.CairoContext(surface)
end

function draw_vectors(ctx::Cairo.CairoContext, vectors, color, factor=1)
  Cairo.save(ctx)
  Cairo.set_source_rgb(ctx, color...)
  for line in vectors
    diff = line[2] - line[1]
    Cairo.move_to(ctx, line[1]...)
    Cairo.line_to(ctx, line[1]+diff*factor...)
  end
  Cairo.set_line_width(ctx, 3.0)
  Cairo.stroke(ctx)
  return ctx
end

function draw_line(ctx::Cairo.CairoContext, line, color=(0,0,0), thickness=2.0)
  Cairo.save(ctx)
  Cairo.set_source_rgb(ctx, color...)
  Cairo.move_to(ctx, line[1]...)
  Cairo.line_to(ctx, line[2]...)
  Cairo.set_line_width(ctx, thickness)
  Cairo.stroke(ctx)
  return ctx
end

function draw_points(ctx::Cairo.CairoContext, points, color=(0,0,0), factor=1)
  d = 5
  Cairo.save(ctx)
  Cairo.set_source_rgb(ctx, color...)
  for point in points
    Cairo.move_to(ctx, point+[-d,-d]...)
    Cairo.line_to(ctx, point+[d,d]...)
    Cairo.move_to(ctx, point+[-d,d]...)
    Cairo.line_to(ctx, point+[d,-d]...)
  end
  Cairo.set_line_width(ctx, 3.0)
  Cairo.stroke(ctx)
  return ctx
end

function draw_point(ctx::Cairo.CairoContext, point, color, factor=1)
  d = 5
  Cairo.save(ctx)
  Cairo.set_source_rgb(ctx, color...)
  Cairo.move_to(ctx, point+[-d,-d]...)
  Cairo.line_to(ctx, point+[d,d]...)
  Cairo.move_to(ctx, point+[-d,d]...)
  Cairo.line_to(ctx, point+[d,-d]...)
  Cairo.set_line_width(ctx, 2.0)
  Cairo.stroke(ctx)
  return ctx
end

function draw_text(ctx::Cairo.CairoContext, text, point, offset, fontsize, color)
  fill_color = color
  line_color = [1,1,1] - color
  Cairo.save(ctx)
  Cairo.select_font_face(ctx, "Sans", Cairo.FONT_SLANT_NORMAL, 
                                          Cairo.FONT_WEIGHT_BOLD)
  Cairo.set_font_size(ctx, fontsize)
  Cairo.move_to(ctx, point+offset...)
  Cairo.rotate(ctx, -pi/2)
  Cairo.scale(ctx, -1, 1)
  Cairo.text_path(ctx, text)
  Cairo.set_source_rgb(ctx, fill_color...)
  Cairo.fill_preserve(ctx)
  Cairo.set_source_rgb(ctx, line_color...)
  Cairo.set_line_width(ctx, 1.0)
  Cairo.scale(ctx, -1, 1)
  Cairo.rotate(ctx, pi/2)
  Cairo.stroke(ctx)
  return ctx
end

function draw_indices(ctx::Cairo.CairoContext, points, offset, fontsize, color)
  fill_color = color
  line_color = [1,1,1] - color
  Cairo.save(ctx)
  Cairo.select_font_face(ctx, "Sans", Cairo.FONT_SLANT_NORMAL, 
                                          Cairo.FONT_WEIGHT_BOLD)
  Cairo.set_font_size(ctx, fontsize)
  for (k, point) in enumerate(points)
    Cairo.move_to(ctx, point+offset...)
    Cairo.rotate(ctx, -pi/2)
    Cairo.scale(ctx, -1, 1)
    Cairo.text_path(ctx, "$k")
    Cairo.set_source_rgb(ctx, fill_color...)
    Cairo.fill_preserve(ctx)
    Cairo.set_source_rgb(ctx, line_color...)
    Cairo.set_line_width(ctx, 1.0)
    Cairo.scale(ctx, -1, 1)
    Cairo.rotate(ctx, pi/2)
  end
  Cairo.stroke(ctx)
  return ctx
end

function get_drawing(surface::Cairo.CairoSurface)
  return convert(Array{RGB24}, surface.data)
end

function draw_reference(ctx, h, w, factor)
  reference = 1
  margin = 10
  a = [h-margin, w-margin-reference*factor]
  b = [h-margin, w-margin]
  color = [0.9,0.9,0.9]
  draw_line(ctx, (a,b), color, 3.0)
  draw_text(ctx, "1", a+[-20, reference*factor/2], [0,0], 24.0, color)
  return ctx
end