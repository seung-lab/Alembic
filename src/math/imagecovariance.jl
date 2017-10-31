immutable SigmaEnv
  	img::Array{Float64, 2}
	II::Array{Int64, 2}
	JJ::Array{Int64, 2}
	IIx::Array{Float64, 2}
	JJx::Array{Float64, 2}
end

global SIGMA_ENVS = Dict{NTuple{2, Int64}, SigmaEnv}()

function get_sigmaenv(img)
  if haskey(SIGMA_ENVS, size(img)) return SIGMA_ENVS[size(img)] end

  SIGMA_ENVS[size(img)] = SigmaEnv(
  		similar(img), 
		[i for i=1:size(img,1), j=1:size(img,2)], 
		[j for i=1:size(img,1), j=1:size(img,2)],
		similar(img),
		similar(img))
  return SIGMA_ENVS[size(img)]
end

function subtract_beta(img, beta = 0.5)
	return max.(img - beta * maximum(img), 0)
end

function subtract_beta!(x, img, beta = 0.5)
        d = beta * maximum(img)
	@simd for i in 1:length(img)
	   @fastmath @inbounds x[i] = img[i] - d
	 end
	 for i in 1:length(img)
	   @inbounds if x[i] < 0 x[i] = 0 end
	 end
end

function sigma(img, beta = 0.5)
        se::SigmaEnv = get_sigmaenv(img)
	subtract_beta!(se.img, img, beta)
#	x = subtract_beta(img, beta)
    if sum(se.img) == 0
        return Inf
    else
      sum_img = sum(se.img)
      @simd for i in 1:length(img)
	@fastmath @inbounds se.img[i] = se.img[i] / sum_img
        end
        return sqrt(sum(diag(ImageCovariance(se))))
	end
end

function ImageCovariance(se::SigmaEnv)
 	s1 = 0.0
  	s2 = 0.0
	@simd for i in 1:length(se.img)
	  @inbounds t1 = se.II[i]
	  @inbounds t2 = se.JJ[i]
	  @inbounds v = se.img[i]
	  @fastmath s1 += v * t1
	  @fastmath s2 += v * t2
	end
	@simd for i in 1:length(se.img)
	  @inbounds t1 = se.II[i]
	  @inbounds t2 = se.JJ[i]
	  @fastmath @inbounds se.IIx[i] = t1 - s1
	  @fastmath @inbounds se.JJx[i] = t2 - s2
	end
    	C=zeros(2,2)
	c11 = 0.0
	c22 = 0.0
	c12 = 0.0
	@simd for i in 1:length(se.img)
	  @inbounds t1 = se.IIx[i]
	  @inbounds t2 = se.JJx[i]
	  @inbounds v = se.img[i]
	  c11 += v * t1 * t1
	  c22 += v * t2 * t2
	  c12 += v * t1 * t2
	end
    	C[1,1]=c11
    	C[2,2]=c22
    	C[1,2]=c12
    	C[2,1]=c12
    	C
end

function ImageCovariance(img)
    # covariance matrix for a 2D image 'img' regarded as
    # a probability distribution over pixel locations
    # regard location in 2D as a random variable
    # 'img' should be nonnegative and normalized like a prob distribution
    II=[i for i=1:size(img,1), j=1:size(img,2)]
    JJ=[j for i=1:size(img,1), j=1:size(img,2)]
    # subtract the mean values of i and j
    II -= sum(img.*II)
    JJ -= sum(img.*JJ)
    C=zeros(2,2)
    C[1,1]=sum(img.*II.^2)
    C[2,2]=sum(img.*JJ.^2)
    C[1,2]=sum(img.*II.*JJ)
    C[2,1]=C[1,2]
    C
end

function softmax(img,beta)
    # https://en.wikipedia.org/wiki/Softmax_function
    # no protection against overflow or underflow
    # larger beta means that 
    prob = exp(beta*img)
    prob /= sum(prob)
end
#=
# example 1
img=gaussian2d(2,[21,21])
ImageCovariance(img)   # this should be roughly 4*eye(2)

# example 2
normxcorroutput = log(gaussian2d(2,[21 21])) # substitute a real normalized cross correlation here
# treat 'normxcorroutput' as proportional to the log probability
C=ImageCovariance(softmax(normxcorroutput,5))

# is the peak narrow?  trace of covariance matrix
println(sqrt(C[1,1]+C[2,2]))  # rms deviation from mean. small value=>narrow

# is the peak anisotropic?  ratio of eigenvalues. large value=>anisotropic
lambda=eigvals(C)
println(maximum(lambda)/minimum(lambda))
=#
