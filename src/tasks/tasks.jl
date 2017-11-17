function match_task(params::Dict)
    load_params(params)
    ms = make_stack()
    println("Match Task $(get_name(ms))")
    match!(ms)
    save(ms, "meshset")
end

"""
Fix the center section, solve in halves, then resolve in range around center
"""
# function solve_task(params::Dict)
# end

function render_task(params::Dict)
    load_params(params)
    for z in get_z_range()
        ms = load("mesh", string(z))
        println("Render Task $(get_name(ms))")
        render(ms)
    end
end

