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

    'Moving' - 1xE binary vector containing indices of moving vertices
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
    Lengths = Array{Float64}(halflen);
    @fastmath @inbounds @simd for ind in 1:halflen
        @fastmath sec_ind = ind + halflen;
	@fastmath @inbounds Lengths[ind] = sqrt(Springs[ind]^2 + Springs[sec_ind]^2) + eps
    end
   return Lengths; 
end

function get_lengths!(Springs, Lengths)
    @fastmath halflen = div(length(Springs), 2);
    sec_ind = 0;
    @fastmath @inbounds @simd for ind in 1:halflen
        @fastmath sec_ind = ind + halflen;
	@fastmath @inbounds spring_first_sq = Springs[ind]^2;
	@fastmath @inbounds spring_second_sq = Springs[sec_ind]^2;
	@fastmath len = sqrt(spring_first_sq + spring_second_sq);
	@fastmath @inbounds Lengths[ind] = len + eps
    end
   return Lengths; 
end

function Energy_given_lengths(Lengths, Stiffnesses, RestLengths)
    #@fastmath dLengths = Lengths - RestLengths
    energy = 0;
    @time @fastmath @inbounds @simd for ind in 1:length(Lengths)
      dlen = Lengths[ind] - RestLengths[ind];
      dlen = dlen * dlen;
      energy += Stiffnesses[ind] * dlen
    end
    return energy / 2
    #@fastmath return sum(Stiffnesses.*(dLengths.*dLengths))/2/length(Lengths)
end

function Energy_given_lengths!(Lengths, Stiffnesses, RestLengths, Energies)
    @simd for ind in 1:length(Lengths)
      @fastmath @inbounds Energies[ind] = Stiffnesses[ind] * (Lengths[ind] - RestLengths[ind])^2
    end
    @fastmath return sum(Energies) / 2
end

#="""
gradient of unnormalized potential energy with respect to vertex positions
returns dxV array, same size as Vertices
physically, -gradient is spring forces acting on vertices
"""=#
function Gradient(Springs, Incidence_d, Stiffnesses_d, RestLengths_d)
    Lengths = get_length(Springs)
    return Gradient_given_lengths(Springs, Lengths, Incidence_d, Stiffnesses_d, RestLengths_d)
end

# vectorised, but slower because Julia is weird
function Gradient_given_lengths(Springs, Lengths, Incidence_t, Stiffnesses_d, RestLengths_d)
    @fastmath Directions = Springs ./ vcat(Lengths, Lengths);
    #Directions[isnan(Directions)] *= 0
    @fastmath Forces = (Springs-(Directions .* RestLengths_d)) .* Stiffnesses_d
    @fastmath return (Incidence_t' * Forces)
end

function Gradient_given_lengths!(Springs, Lengths, Incidence_d, Stiffnesses_d, RestLengths_d, Forces, Gradients)
  # @time begin
     #Forces = Array{Float64}(length(Springs))
     len = length(Lengths);
     @fastmath @inbounds @simd for ind in 1:len
	sec_ind = ind + len;
	@fastmath @inbounds Forces[ind] = (Springs[ind] - (Springs[ind] / Lengths[ind] * RestLengths_d[ind])) * Stiffnesses_d[ind]
	@fastmath @inbounds Forces[sec_ind] = (Springs[sec_ind] - (Springs[sec_ind] / Lengths[ind] * RestLengths_d[sec_ind])) * Stiffnesses_d[sec_ind]
     end
 #end
  #print("Storage:")
    @fastmath @inbounds A_mul_B!(1.0, Incidence_d, Forces, 0.0, Gradients)
end

function SolveMeshConjugateGradient!(Vertices, Fixed, Incidence, Stiffnesses, RestLengths, max_iter, ftol)

    N = size(Vertices, 2);
    perm = sortperm(Fixed);
    invperm = sortperm(perm);

    # permute so that the first M things are moving = this allows subarray access to be dense, hence vectorizable
    Vertices[:] = Vertices[:,perm]
    Fixed = Fixed[perm]
    Incidence = Incidence[perm,:]
    M = sum(!Fixed)	# number of moving things

    # double everything to make everything 1-dimensional
    Vertices_t = Vertices';
    Vertices_t = vcat(Vertices_t[:, 1], Vertices_t[:, 2])

    Incidence_t = Incidence'
    Incidence_t = vcat(hcat(Incidence_t, spzeros(size(Incidence_t)...)), hcat(spzeros(size(Incidence_t)...), Incidence_t))
    Incidence_d = Incidence_t'

    Lengths = Array{Float64}(div(size(Incidence_t, 1), 2));
    Energies = Array{Float64}(div(size(Incidence_t, 1), 2));
    Springs = Array{Float64}(size(Incidence_t, 1));
    Forces = Array{Float64}(size(Incidence_t, 1));
    Gradients = Array{Float64}(size(Vertices_t, 1));

    Moving = vcat(~Fixed, ~Fixed)
    Stiffnesses_d = vcat(Stiffnesses, Stiffnesses)
    RestLengths_d = vcat(RestLengths, RestLengths)
#=
    cost_iter = 0;
    grad_iter = 0;
    cost_grad_iter = 0;

    cost_time = 0.0;
    grad_time = 0.0;
    cost_grad_time = 0.0;
    =#


    function cost(x)
      @time begin
#	tic()
        @inbounds Vertices_t[Moving] = x;
        #Springs = Incidence_t * Vertices_t;
        @fastmath @inbounds A_mul_B!(1.0, Incidence_t, Vertices_t, 0.0, Springs)
     #   Springs = Incidence_d' * Vertices_t;
	@fastmath get_lengths!(Springs, Lengths);	
        @fastmath energy = Energy_given_lengths!(Lengths,Stiffnesses,RestLengths, Energies)
	print("cost: ")
#	cost_iter += 1;
#	cost_time += toc();
      end
      return energy
    end

    # NOT USED EVER - MIGHT BE INCORRECT
    function cost_gradient!(x,storage)
      @time begin
#	tic()
        @inbounds Vertices_t[Moving] = x;
        @fastmath @inbounds A_mul_B!(1.0, Incidence_t, Vertices_t, 0.0, Springs)
        #Springs = Incidence_t * Vertices_t;
      # Springs = Incidence_d' * Vertices_t;
        g = Gradient(Springs, Incidence_d, Stiffnesses_d, RestLengths_d)
        storage[:] = g
	print("gradient: ")
#	grad_iter += 1;
#	grad_time += toc();
    end
    end

    function cost_and_gradient!(x,storage)
      @time begin
#	tic()
        @inbounds Vertices_t[Moving] = x;
	#print("Springs: ")
        @fastmath @inbounds A_mul_B!(1.0, Incidence_t, Vertices_t, 0.0, Springs)
        #fastmath Springs = Incidence_t * Vertices_t;
      # @fastmath Springs = Incidence_d' * Vertices_t;
    	@fastmath get_lengths!(Springs, Lengths);
        @fastmath @inbounds Gradient_given_lengths!(Springs, Lengths, Incidence_d, Stiffnesses_d, RestLengths_d, Forces, Gradients)
	#@fastmath @inbounds storage[:] = Gradients[Moving]
	@fastmath @inbounds copy!(view(storage, 1:M), view(Gradients,1:M));
	@fastmath @inbounds copy!(view(storage, M+1:2*M), view(Gradients, N+1:N+M));
	@fastmath energy = Energy_given_lengths!(Lengths,Stiffnesses,RestLengths, Energies);
	print("cost_and_gradient: ")
#	cost_grad_iter += 1;
#	cost_grad_time += toc();
      end
        return energy;
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

    #@fastmath res = optimize(df,Vertices_t[Moving],method=ConjugateGradient(),show_trace=true,iterations=max_iter,ftol=ftol)
    @fastmath res = optimize(df,Vertices_t[Moving],method=ConjugateGradient(),show_trace=true,iterations=max_iter,ftol=ftol)
    # return res
   # Vertices_t[Moving] = res.minimum;
    copy!(view(Vertices_t,1:M), view(res.minimum, 1:M));
    copy!(view(Vertices_t,N+1:N+M), view(res.minimum, M+1:2*M));
    Vertices[:] = vcat(Vertices_t[1:div(length(Vertices_t),2)]', Vertices_t[1+div(length(Vertices_t),2):end]');
    Vertices[:] = Vertices[:,invperm];
    #println("cost_iter: $cost_iter, cost_time: $cost_time")
    #println("grad_iter: $grad_iter, grad_time: $grad_time")
    #println("cost_grad_iter: $cost_grad_iter, cost_grad_time: $cost_grad_time")
    
end
