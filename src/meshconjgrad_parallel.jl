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

function Gradient_given_lengths!(Springs_s, Lengths, Incidence_t, Stiffnesses_d, RestLengths_d, Forces_s, Gradients_s)
    @fastmath Directions = Springs_s ./ vcat(Lengths, Lengths);
    #Directions[isnan(Directions)] *= 0
#    @fastmath Forces = (Springs.-(Directions .* RestLengths_d)) .* Stiffnesses_d
    @fastmath Forces_s[:] = (Springs_s-(Directions .* RestLengths_d)) .* Stiffnesses_d
    @fastmath mul_by_inctt!(Forces_s, Gradients_s)
end

#from julia parallel example
function myrange(q::SharedArray)
        idx = indexpids(q)
	    if idx == 0
	              # This worker is not assigned a piece
	return 1:0
		end
	nchunks = length(procs(q))
	splits = [round(Int, s) for s in linspace(0,size(q,1),nchunks+1)]
	return splits[idx]+1:splits[idx+1]
end

function setup_local_incidences(Incidence_d, Springs_s, Incidence_t, Gradients_s)
	global LOCAL_INC_D = sdata(Incidence_d)[:, myrange(Springs_s)]
	global LOCAL_INC_T = sdata(Incidence_t)[:, myrange(Gradients_s)]
end



function mul_by_incdt_chunked!(B, R, ran = myrange(R))
  @fastmath R[ran] = LOCAL_INC_D' * B;
end
function mul_by_incdt!(B::SharedArray, R::SharedArray)
  	@sync for p in procs()
		@async remotecall_wait(p, mul_by_incdt_chunked!, sdata(B), R)
	end
end

function mul_by_inctt_chunked!(B, R, ran = myrange(R))
  @fastmath R[ran] = LOCAL_INC_T' * B;
end

function mul_by_inctt!(B::SharedArray, R::SharedArray)
  	@sync for p in procs()
		@async remotecall_wait(p, mul_by_inctt_chunked!, sdata(B), R)
	end
end

function At_mul_B_to_R_chunked!(A::SharedSparseMatrixCSC, B::SharedArray, R::SharedArray, ran = myrange(R))
	R[ran] = sdata(A)[:, ran]' * B
end


function At_mul_B_to_R!(A::SharedSparseMatrixCSC, B::SharedArray, R::SharedArray)
  	@sync for p in procs()
		@async remotecall_wait(p, mul_by_incdt_chunked!, B, R)
	end
end

function SolveMeshConjugateGradient!(Vertices, Fixed, Incidence, Stiffnesses, RestLengths, max_iter, ftol)
    # double everything
    Vertices_t = Vertices';
    Vertices_t = vcat(Vertices_t[:, 1], Vertices_t[:, 2])

    #Incidence_t = Incidence'
    #Incidence_t = vcat(hcat(Incidence_t, spzeros(size(Incidence_t)...)), hcat(spzeros(size(Incidence_t)...), Incidence_t))
    Incidence_d = hcat(vcat(Incidence, spzeros(size(Incidence)...)), vcat(spzeros(size(Incidence)...), Incidence))
    Incidence_t = Incidence_d';
    Incidence_d = share(Incidence_d);
    Incidence_t = share(Incidence_t);


    Moving = vcat(~Fixed, ~Fixed)
    Stiffnesses_d = vcat(Stiffnesses, Stiffnesses)
    RestLengths_d = vcat(RestLengths, RestLengths)

    Springs_s = SharedArray(Float64, size(Incidence_d, 2))
    Forces_s = SharedArray(Float64, size(Incidence_d, 2))
    Gradients_s = SharedArray(Float64, size(Incidence_d, 1))
	@sync for p in procs()
	  @async remotecall_wait(p, setup_local_incidences, Incidence_d, Springs_s, Incidence_t, Gradients_s)
	end
    println("incidences shared")


    function cost(x)
        Vertices_t[Moving] = x[Moving];
        Springs = Incidence_t * Vertices_t;
        return Energy(Springs,Stiffnesses,RestLengths)
    end

    function cost_gradient!(x,storage)
        Vertices_t[Moving] = x[Moving];
        Springs = Incidence_t * Vertices_t;
        g = Gradient(Springs, Incidence_d, Stiffnesses_d, RestLengths_d)
        storage[:] = g
    end

    function cost_and_gradient!(x,storage)
        #@inbounds Vertices_t[Moving] = x[Moving];
	if typeof(Vertices_t) == Array{Float64, 1}
	  Vertices_t = share(Vertices_t);
	end
        @inbounds Vertices_t[Moving] = x[Moving];
        @fastmath mul_by_incdt!(Vertices_t, Springs_s);
    	@fastmath Lengths = get_lengths(Springs_s);
        Gradient_given_lengths!(Springs_s, Lengths, Incidence_t, Stiffnesses_d, RestLengths_d, Forces_s, Gradients_s)
	@inbounds storage[:] = Gradients_s[:]
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

    @fastmath res = optimize(df,Vertices_t,method=:cg,show_trace=true,iterations=max_iter,ftol=ftol)
    # return res
    Vertices_t[Moving] = res.minimum[Moving];
    Vertices[:] = vcat(Vertices_t[1:(length(Vertices_t)/2)]', Vertices_t[1+(length(Vertices_t)/2):end]');
    
end
