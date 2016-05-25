#="""
MeshConjGrad- given spring mesh, solve for equilibrium positions of vertices with nonlinear conjugate gradient

    V = # mesh vertices in R^d
    E = # of springs

    'Vertices' - dxV matrix, columns contain vertex positions

    'Incidence' - VxE generalized oriented incidence matrix
    springs <-> columns
    intratile spring <-> column containing 1 -1
    intertile spring <-> column containing w1 w2 w3 -1 where (w1, w2, w3) represents a weighted combination of vertices

    most functions compute Springs=Vertices*Incidences, dxE matrix
    spring vectors <-> columns

    'Stiffnesses', 'RestLengths' - 1xE vectors specifying spring properties

    'Moving' - integer vector containing indices of moving vertices
 could be changed to 1xE binary vector
"""
=#
#currently defined in Julimaps
#global eps = 1E-16

function Energy(Springs, Stiffnesses, RestLengths)
    # potential energy per spring (normalized)
    Lengths=get_lengths(Springs)
    return Energy_given_lengths(Lengths, Stiffnesses, RestLengths)
end

function get_lengths(Springs)
    @fastmath halflen = div(length(Springs), 2);
    r1 = 1:halflen
    r2 = halflen + 1:halflen * 2
    @fastmath @inbounds return sqrt(Springs[r1] .* Springs[r1] + Springs[r2] .* Springs[r2]) + eps
end

function Energy_given_lengths(Lengths, Stiffnesses, RestLengths)
    @fastmath dLengths = Lengths - RestLengths
    @fastmath return sum(Stiffnesses.*(dLengths.*dLengths))/2/length(Lengths)
end

#="""
gradient of unnormalized potential energy with respect to vertex positions
returns dxV array, same size as Vertices
physically, -gradient is spring forces acting on vertices
"""=#
function Gradient(Springs, Incidence_d, Stiffnesses_d, RestLengths_d)
    print(".");
    Lengths = get_length(Springs)
    return Gradient_given_lengths(Springs, Lengths, Incidence_d, Stiffnesses_d, RestLengths_d)
end

function Gradient_given_lengths(Springs, Lengths, Incidence_t, Stiffnesses_d, RestLengths_d)
    @fastmath Directions = Springs ./ vcat(Lengths, Lengths);
    #Directions[isnan(Directions)] *= 0
#    @fastmath Forces = (Springs.-(Directions .* RestLengths_d)) .* Stiffnesses_d
    @fastmath Forces = (Springs-(Directions .* RestLengths_d)) .* Stiffnesses_d
    @fastmath return (Incidence_t' * Forces)
end

function SolveMeshConjugateGradient!(Vertices, Fixed, Incidence, Stiffnesses, RestLengths, max_iter, ftol)
    # double everything
    Vertices_t = Vertices';
    Vertices_t = vcat(Vertices_t[:, 1], Vertices_t[:, 2])

    Incidence_t = Incidence'
    Incidence_t = vcat(hcat(Incidence_t, spzeros(size(Incidence_t)...)), hcat(spzeros(size(Incidence_t)...), Incidence_t))
    Incidence_d = Incidence_t'

    Moving = vcat(~Fixed, ~Fixed)
    Stiffnesses_d = vcat(Stiffnesses, Stiffnesses)
    RestLengths_d = vcat(RestLengths, RestLengths)


    function cost(x)
        Vertices_t[Moving] = x[:];
        Springs = Incidence_t * Vertices_t;
        return Energy(Springs,Stiffnesses,RestLengths)
    end

    function cost_gradient!(x,storage)
        Vertices_t[Moving] = x[:];
        Springs = Incidence_t * Vertices_t;
        g = Gradient(Springs, Incidence_d, Stiffnesses_d, RestLengths_d)
        storage[:] = g
    end

    function cost_and_gradient!(x,storage)
        @inbounds Vertices_t[Moving] = x[:];
        @fastmath Springs = Incidence_d' * Vertices_t;
    	@fastmath Lengths = get_lengths(Springs);
        @inbounds storage[:] = Gradient_given_lengths(Springs, Lengths, Incidence_t, Stiffnesses_d, RestLengths_d)[Moving]
        return Energy_given_lengths(Lengths,Stiffnesses,RestLengths)
    end

    df = DifferentiableFunction(cost, cost_gradient!, cost_and_gradient!)

    #    function cost_and_gradient!(x,storage)
    #        Vert = copy(Vertices)
    #        Vert[:,Moving]=reshape(x,2,div(length(x),2))
    #        Springs=Vert*Incidence
    #        g = Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    #        storage[:] = g[:,Moving][:]
    #        return Energy(Springs,Stiffnesses,RestLengths)
    #    end
    #    df = DifferentiableFunction(cost, cost_gradient!, cost_and_gradient!)

    @fastmath res = optimize(df,Vertices_t[Moving],method=ConjugateGradient(),show_trace=true,iterations=max_iter,ftol=ftol)
    # return res
    Vertices_t[Moving] = res.minimum[:];
    Vertices[:] = vcat(Vertices_t[1:(length(Vertices_t)/2)]', Vertices_t[1+(length(Vertices_t)/2):end]');
    
end
