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
        dst_z = get(Z_MAP, z, z)
        make_mips(dst_z; mips=get_downsample_mips(:dst_image))
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

function mip_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        println("Mip Task $(z)")
        make_mips(z; mips=get_downsample_mips(:dst_image))
        reset_cache()
    end
end

function mask_union_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        println("Mask Union Task $(z)")
        offset = get_offset(:roi_mask, mip=get_mip(:dst_image))
        roi_mask = get_image(z, :roi_mask, mip=get_mip(:dst_image), input_mip=get_mip(:roi_mask))
        defect_mask = get_image(z, :defect_mask, mip=get_mip(:dst_image), input_mip=get_mip(:defect_mask))
        mask_union = convert(Array{UInt8,2}, (roi_mask + defect_mask) .> 0)
        @time save_image(:dst_image, mask_union, offset, z, mip=get_mip(:dst_image), cdn_cache=false)
        make_mips(z; mips=get_downsample_mips(:dst_image))
        reset_cache()
    end
end  

function dilation_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    kwargs = symbol_dict(params["task"]["kwargs"])
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        println("Dilation Task $(z)")
        dilate(z; kwargs...);
        make_mips(z; mips=get_downsample_mips(:dst_image))
        reset_cache()
    end
end    

function threshold_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    kwargs = symbol_dict(params["task"]["kwargs"])
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        println("Threshold Task $(z)")
        threshold(z; kwargs...);
        make_mips(z; mips=get_downsample_mips(:dst_image))
        reset_cache()
    end
end    

function mask_rectangle_task(params::Dict)
    load_params(params)
    pairs = params["task"]["pairs"]
    kwargs = symbol_dict(params["task"]["kwargs"])
    z_indices = unique(vcat(pairs...))
    for z in z_indices
        println("Mask Rectangle Task $(z)")
        mip = get_mip(:dst_image)
        input_mip = params["task"]["kwargs"]["mip"]
        ranges = parse_ranges(params["task"]["kwargs"]["ranges"], mip=mip, input_mip=input_mip)
        val = UInt8(params["task"]["kwargs"]["val"])
        path = get_path(:dst_image)
        cv = CloudVolume.CloudVolumeWrapper(path, mip=mip, 
                        bounded=true, autocrop=true,
                        cdn_cache=false, progress=false, parallel=PARALLEL,
                        output_to_shared_memory=false, non_aligned_writes=true)
        for r in ranges
            img = ones(UInt8, length(r[1]), length(r[2]), 1, 1)*val
            print((r..., z:z, size(img), img[1]))
            cv[r..., z:z] = img
        end
        make_mips(z; mips=get_downsample_mips(:dst_image))
        reset_cache()
    end
end

"""
Parse task kwargs for ranges: [[x_start, x_stop], [y_start, y_stop]]
"""
function parse_ranges(arr::Array; mip=get_mip(:dst_image), input_mip=get_mip(:dst_image))
    s = get_scale(mip) / get_scale(input_mip)
    o = []
    for a in arr
        x_start = round(Int, a[1][1]*s)
        x_stop = round(Int, a[1][2]*s)
        y_start = round(Int, a[2][1]*s)
        y_stop = round(Int, a[2][2]*s)
        push!(o, (x_start:x_stop, y_start:y_stop))
    end
    return o
end

function symbol_dict(d::Dict)
    return Dict(Symbol(k) => v for (k,v) in d)
end
