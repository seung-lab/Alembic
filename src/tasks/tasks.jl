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

