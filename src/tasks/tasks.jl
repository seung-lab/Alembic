function match_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    z_indices = unique(vcat(pairs...))
    println("Match Task $(pairs)")
    ms = make_stack(z_indices)
    match!(ms, pairs)
    split_meshset(ms)
    reset_cache()
end

function render_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        ms = load(:mesh, string(z))
        println("Render Task $(z)")
        render(ms)
        make_mips(z; mips=get_downsample_mips(:dst_image))
        reset_cache()
    end
end

function mask_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        println("Mask Task $(z)")
        src_image = get_image(z, :src_image, mip=get_mip(:dst_image), input_mip=get_mip(:src_image))
        src_offset = get_offset(:src_image, mip=get_mip(:dst_image))
        if use_roi_mask()
          src_roi = get_image(z, :roi_mask, mip=get_mip(:dst_image), input_mip=get_mip(:roi_mask))
          roi_value = get_mask_value(:roi_mask)
          unsafe_mask_image!(src_image, src_roi, roi_value, src_image, keep_id=true)
        end
        @time save_image(:dst_image, Array(src_image), src_offset, z, mip=get_mip(:dst_image))
        reset_cache()
    end
end
