import PyPlot
#using HDF5
#using JLD

#data=load("r4c2-r4c3.jld")
#Vertices=data["nodes"]
#Incidence=data["edges"]
#RestLengths=data["edge_lengths"]
#Stiffnesses=data["edge_coeffs"]
#Fixed=data["nodes_fixed"]

#fid = jldopen("./solvedMesh.jld", "w")
#fid["n"] = 3317;
#fid["m"] = 9781;
#fid["nodes"] = Vertices[:, 1:3317];

#Vertices=convert(Array{Float64, 2},Vertices)
#Incidence=convert(SparseMatrixCSC{Float64, Int64}, Incidence)
#RestLengths=convert(Array{Float64, 1}, RestLengths)
#Stiffnesses=convert(Array{Float64, 1}, Stiffnesses)

Vertices = hcat(FMs.nodes...);
Incidence = FMs.edges;
RestLengths = FMs.edge_lengths;
Stiffnesses = FMs.edge_coeffs;
Fixed = FMs.nodes_fixed;

d=size(Vertices,1)
V=size(Vertices,2)
E=size(Incidence,2)
Lengths=zeros(1,V)

println("$d, $V, $E")

Moving = ~Fixed

#function MeshSolve(Vertices, Incidence, Stiffnesses, RestLengths, Moving)
Moving2=[Moving[:]'; Moving[:]'][:]   # double the dimensionality

#Vertices[:,Moving]=Vertices[:,Moving]+20*randn(size(Vertices[:,Moving]))
eta=0.25
niter=250
U=zeros(1,niter)     # energy vs. time
g=similar(Vertices)  # gradient of potential energy

for iter=1:niter
    Springs=Vertices*Incidence
    g=Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    if iter<10
        # gradient descent
        Vertices[:,Moving]=Vertices[:,Moving]-eta*g[:,Moving]
    else
        #  Newton's method
        #=
        H=Hessian(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta*reshape(sparse(H[Moving2,Moving2])\g[:,Moving][:],2,length(find(Moving)))=#

H=Hessian2(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta*reshape(H[Moving2,Moving2]\g[:,Moving][:],2,length(find(Moving)))

#=        H=Hessian2(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[Moving]=Vertices[Moving]-eta*reshape(H[Moving2,Moving2]\g[:,Moving][:],2,length(find(Moving)))=#
    end
    U[iter]=Energy(Springs,Stiffnesses,RestLengths)
    println(iter," ", U[iter])
    #    visualize the dynamics
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

#=
fid["nodes_t"] = Vertices[:, 1:3317];
fid["edges"] = Incidence[1:3317, 1:9781]
=#
close(fid);


