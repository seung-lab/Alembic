"""
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

global eps = 1E-8

function Energy( Springs, Stiffnesses, RestLengths)
    # potential energy per spring (normalized)
    Lengths=sqrt(sum(Springs.^2,1))   # spring lengths (row vector)
    sum(Stiffnesses[:].*(Lengths[:]-RestLengths[:]).^2)/2/size(Springs,2)
end

function Gradient( Springs, Incidence, Stiffnesses, RestLengths)
    # gradient of unnormalized potential energy with respect to vertex positions
    # returns dxV array, same size as Vertices
    # physically, -gradient is spring forces acting on vertices
    d=size(Springs,1)
    Lengths=sqrt(sum(Springs.^2,1)) + eps
    Directions=broadcast(/,Springs,Lengths)
    Directions[isnan(Directions)] *= 0
    Forces=broadcast(*,Stiffnesses[:]',Springs-broadcast(*,RestLengths[:]',Directions))
    Forces*Incidence'
end

function SolveMesh2!(Vertices, Fixed, Incidence, Stiffnesses, RestLengths)

    function cost(x)
        Vert = copy(Vertices)
        Vert[:,Moving]=reshape(x,2,div(length(x),2))
        Springs=Vert*Incidence
        return Energy(Springs,Stiffnesses,RestLengths)
    end

    function cost_gradient!(x,storage)
        Vert = copy(Vertices)
        Vert[:,Moving]=reshape(x,2,div(length(x),2))
        Springs=Vert*Incidence
        g = Gradient(Springs, Incidence, Stiffnesses, RestLengths)
        storage[:] = g[:,Moving][:]
    end
    df = DifferentiableFunction(cost, cost_gradient!)

    #    function cost_and_gradient!(x,storage)
    #        Vert = copy(Vertices)
    #        Vert[:,Moving]=reshape(x,2,div(length(x),2))
    #        Springs=Vert*Incidence
    #        g = Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    #        storage[:] = g[:,Moving][:]
    #        return Energy(Springs,Stiffnesses,RestLengths)
    #    end
    #    df = DifferentiableFunction(cost, cost_gradient!, cost_and_gradient!)

    Moving = ~Fixed

    res = optimize(df,Vertices[:,Moving][:],method=:cg,show_trace=true,ftol=1e-8)
    # return res
    Vertices[:,Moving] = reshape(res.minimum, 2, div(length(res.minimum),2))
end

