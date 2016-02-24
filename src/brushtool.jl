# Pulled from ImageView:rubberband.jl
# https://github.com/timholy/ImageView.jl/blob/master/src/rubberband.jl
#
# For rubberband, we draw the selection region on the front canvas, and repair
# by copying from the back. Note that the front canvas has
#     user coordinates = device coordinates,
# so device_to_user doesn't do anything. The back canvas has
#     user coordinates = image pixel coordinates,
# so is the correct reference for anything dealing with image pixels.
type BrushTool
    pos1
    pos2
    moved::Bool
end

function btdraw(r::Graphics.GraphicsContext, bt::BrushTool)
    rectangle(r, bt.pos1[1], bt.pos1[2], bt.pos2[1]-bt.pos1[1], bt.pos2[2]-bt.pos1[2])
    set_line_width(r, 1)
    set_dash(r, [3.0,3.0], 3.0)
    set_source_rgb(r, 1, 1, 1)
    stroke_preserve(r)
    set_dash(r, [3.0,3.0], 0.0)
    set_source_rgb(r, 0, 0, 0)
    stroke_preserve(r)
end

# callback_done is executed when the user finishes drawing the rubberband.
# Its syntax is callback_done(canvas, boundingbox), where the boundingbox is
# in user coordinates.
function brushtool_start(c::Tk.Canvas, x, y, callback_done::Function)
    # Copy the surface to another buffer, so we can repaint the areas obscured by the rubberband
    println("bt start")
    r = Graphics.getgc(c)
    Graphics.save(r)
    reset_transform(r)
    ctxcopy = copy(r)
    bt = BrushTool((x,y), (x,y), false)
    callbacks_old = (c.mouse.motion, c.mouse.button2release)
    c.mouse.motion = (c, x, y) -> brushtool_move(c, bt, x, y, ctxcopy)
    c.mouse.button2release = (c, x, y) -> brushtool_stop(c, bt, x, y, ctxcopy, callbacks_old, callback_done)
end

function brushtool_move(c::Tk.Canvas, bt::BrushTool, x, y, ctxcopy)
    r = Graphics.getgc(c)
    if bt.moved
        # Erase the previous rubberband by copying from back surface to front
        Cairo.set_source(r, ctxcopy)
        # Since the path was already created and preserved, we just modify its properties
        Cairo.set_line_width(r, 2)
        Cairo.set_dash(r, Float64[])
        Cairo.stroke(r)
    end
    bt.moved = true
    # Draw the new rubberband
    bt.pos2 = (x, y)
    btdraw(r, bt)
    Tk.reveal(c)
    Tk.update()
end

function brushtool_stop(c::Tk.Canvas, bt::BrushTool, x, y, ctxcopy, callbacks_old, callback_done)
    println("bt stop")
    c.mouse.motion = callbacks_old[1]
    c.mouse.button2release = callbacks_old[2]
    if !bt.moved
        return
    end
    r = Graphics.getgc(c)
    Cairo.set_source(r, ctxcopy)
    Cairo.set_line_width(r, 2)
    Cairo.stroke(r)
    Tk.reveal(c)
    Cairo.restore(r)
    Tk.update()
    x1, y1 = bt.pos1[1], bt.pos1[2]
    if abs(x1-x) > 2 || abs(y1-y) > 2
        # It moved sufficiently, let's execute the callback
        xu, yu = Graphics.device_to_user(r, x, y)
        x1u, y1u = Graphics.device_to_user(r, x1, y1)
        bb = Graphics.BoundingBox(min(x1u,xu), max(x1u,xu), min(y1u,yu), max(y1u,yu))
        println("moved sufficiently: ", join((x1u, y1u, xu, yu), ", "))
        callback_done(c, bb)
    end
end