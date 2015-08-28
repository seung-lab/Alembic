# MeshSolve - given spring mesh, solve for equilibrium positions of vertices
# 
# V = # mesh vertices in R^d
# E = # of springs
#
# 'Vertice' - dxV matrix, columns contain vertex positions
#
# 'Incidence' - VxE generalized oriented incidence matrix
#    springs <-> columns
#    intratile spring <-> column containing 1 -1
#    intertile spring <-> column containing w1 w2 w3 -1 where (w1, w2, w3) represents a weighted combination of vertices
#
#  most functions compute Springs=Vertices*Incidences, dxE matrix
#    spring vectors <-> columns

# 'Stiffnesses', 'RestLengths' - 1xE vectors specifying spring properties

# 'Moving' - integer vector containing indices of moving vertices
# could be changed to 1xE binary vector

import PyPlot

global eps = 1E-8

function Energy( Springs, Stiffnesses, RestLengths)
    # potential energy in springs
    Lengths=sqrt(sum(Springs.^2,1))   # spring lengths (row vector)
    sum(Stiffnesses[:].*(Lengths[:]-RestLengths[:]).^2)/2
end

function Gradient( Springs, Incidence, Stiffnesses, RestLengths)
    # gradient of energy with respect to vertex positions
    # returns dxV array, same size as Vertices
    # physically, -gradient is spring forces acting on vertices
    d=size(Springs,1)
    Lengths=sqrt(sum(Springs.^2,1)) + eps
    Directions=broadcast(/,Springs,Lengths)
    Directions[isnan(Directions)] *= 0
    Forces=broadcast(*,Stiffnesses[:]',Springs-broadcast(*,RestLengths[:]',Directions))
    Forces*Incidence'
end

function Hessian( Springs, Incidence, Stiffnesses, RestLengths)
    # Hessian of the potential energy as an Vd x Vd matrix
    # i.e. VxV block matrix of dxd blocks
    # Note: symmetric positive definite
    V = size(Incidence,1)
    d = size(Springs,1)
    H = zeros(V*d, V*d)

    Lengths=sqrt(sum(Springs.^2,1))

    for a=1:size(Springs,2)
        if Lengths[a]==0
            dH = eye(d)   # technically NaN if RestLengths[a]!=0
        else
            Direction=Springs[:,a]/Lengths[a]
            dH = eye(d)-RestLengths[a]/Lengths[a]*(eye(d)-Direction*Direction')
        end
        dH = Stiffnesses[a]*dH;
        VertexList=find(Incidence[:,a])    # vertices incident on spring a
        for i=VertexList
            for j=VertexList
                # indices of (i,j) block of Hessian
                irange=(i-1)*d+(1:d)
                jrange=(j-1)*d+(1:d)
                H[ irange, jrange ] += Incidence[i,a]*Incidence[j,a]*dH
            end
        end
    end  
    H
end

function Hessian2( Springs, Incidence, Stiffnesses, RestLengths)
    # Hessian of the potential energy as an Vd x Vd matrix
    # i.e. VxV block matrix of dxd blocks
    # Note: symmetric positive definite
    E = size(Incidence,2)   # number of springs
    d = size(Springs,1)     # linear size of blocks in block matrix
    d2 = d^d                # number of elements in block
    # allocate space for sparse matrix
    maxel=16*E*d2   # 16 = square of maximum number of vertices per spring
    II=zeros(Int64,maxel)
    JJ=zeros(Int64,maxel)
    SS=zeros(Float64,maxel)

    Lengths=sqrt(sum(Springs.^2,1))
    numel=0
    for a=1:size(Springs,2)
        if Lengths[a]==0
            dH = eye(d)   # technically NaN if RestLengths[a]!=0
        else
            Direction=Springs[:,a]/Lengths[a]
            dH = eye(d)-RestLengths[a]/Lengths[a]*(eye(d)-Direction*Direction')
        end
        dH = Stiffnesses[a]*dH;
        VertexList=find(Incidence[:,a])    # vertices incident on spring a
        for i=VertexList
            for j=VertexList
                # indices and values of (i,j) block of Hessian
                II[numel+(1:d2)]=[id for id=(i-1)*d+(1:d), jd=(j-1)*d+(1:d)][:]
                JJ[numel+(1:d2)]=[jd for id=(i-1)*d+(1:d), jd=(j-1)*d+(1:d)][:]
                SS[numel+(1:d2)]=(Incidence[i,a]*Incidence[j,a]*dH)[:]
                numel += d2
            end
        end
    end
    sparse(II[1:numel],JJ[1:numel],SS[1:numel])
end

function SolveMesh!(Vertices, Fixed, Incidence, Stiffnesses, RestLengths, eta_gradient, eta_newton, ftol_grad, ftol_newton)
    d=size(Vertices,1)
    V=size(Vertices,2)
    E=size(Incidence,2)
    Lengths=zeros(1,V)
    Moving = ~Fixed
    Moving2=[Moving[:]'; Moving[:]'][:]   # double the dimensionality
    U=Array{Float64, 1}(0);     # energy vs. time
    g=similar(Vertices)  # gradient of potential energy

    iter = 1;

    while true
        Springs=Vertices*Incidence
        g=Gradient(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta_gradient*g[:,Moving]
        push!(U, Energy(Springs,Stiffnesses,RestLengths))
        println(iter," ", U[iter])
        if iter != 1
    	if (U[iter-1] - U[iter]) / U[iter-1] < ftol_grad
    		println("Switching to Newton's Method:");    iter += 1; break;
    	end
    	end
        iter += 1;
    end

    while true
    	Springs=Vertices*Incidence
    	g=Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    	H=Hessian2(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta_newton*reshape(H[Moving2,Moving2]\g[:,Moving][:],2,length(find(Moving)))
    	push!(U, Energy(Springs,Stiffnesses,RestLengths))
        println(iter," ", U[iter])
        if (U[iter-1] - U[iter]) / U[iter-1] < ftol_newton
            println("Converged below ", ftol_newton); break;
        end
        iter+=1;
    end
end

#=
function SolveMesh!(Vertices, Fixed, Incidence, Stiffnesses, RestLengths, eta_grad, grad_threshold, eta_newton, newton_threshold, show_plot)

d=size(Vertices,1)
V=size(Vertices,2)
E=size(Incidence,2)
Lengths=zeros(1,V)
Moving = ~Fixed
Moving2=[Moving[:]'; Moving[:]'][:]   # double the dimensionality
U=zeros(1,niter)     # energy vs. time
g=similar(Vertices)  # gradient of potential energy

for iter=1:niter
    Springs=Vertices*Incidence
    g=Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    if iter<ngrad
        # gradient descent
        Vertices[:,Moving]=Vertices[:,Moving]-eta*g[:,Moving]
    else
        #  Newton's method
	H=Hessian2(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta*reshape(H[Moving2,Moving2]\g[:,Moving][:],2,length(find(Moving)))

    end
    U[iter]=Energy(Springs,Stiffnesses,RestLengths)
    println(iter," ", U[iter])
    #    visualize the dynamics
	if(show_plot)
    PyPlot.subplot(221)
    PyPlot.cla()
    PyPlot.scatter(Vertices[1,:],Vertices[2,:])
    PyPlot.subplot(222)
    PyPlot.plot(1:iter,U[1:iter])
    PyPlot.subplot(223)
    Lengths=sqrt(sum(Springs.^2,1))
    PyPlot.cla()
    PyPlot.plot(1:E,Lengths')
    PyPlot.draw()
	end
end

end
=#

