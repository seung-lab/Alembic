# Immutable type that contains the preallocated space needed for a convolution between two objects of a given size, as well as the precomuted plans for the rfft / irfft.
immutable ConvolveEnv
  size_A::NTuple{2, Int64}
  size_B::NTuple{2, Int64}

  crop::Symbol
  padding::Symbol

  intermediate_A::Array{Float64, 2}
  intermediate_B::Array{Float64, 2}

  complex_intermediate_A::Array{Complex{Float64}, 2}
  complex_intermediate_B::Array{Complex{Float64}, 2}

  rfft_plan
  irfft_plan

  result_ranges::Array{UnitRange{Int64},1}
  result::Array{Float64, 2}
end

immutable Normxcorr2Env
  size_A::NTuple{2, Int64}
  size_B::NTuple{2, Int64}

  conv_dt::Array{Float64, 2}
  conv_imgpad::Array{Float64, 2}

  conv_sum::Array{Float64, 2}
  conv_sum2::Array{Float64, 2}

  local_sum::Array{Float64, 2}
  local_sum2::Array{Float64, 2}
end

global CONVOLVE_ENVS = Dict{Symbol, ConvolveEnv}()
global NORMXCORR2_ENVS = Dict{Symbol, Normxcorr2Env}()

function clear_convolveenvs()
	global CONVOLVE_ENVS = Dict{Symbol, ConvolveEnv}()
end
function clear_normxcorr2envs()
	global CONVOLVE_ENVS = Dict{Symbol, Normxcorr2Env}()
end

function register_convolveenv(ce::ConvolveEnv)
	CONVOLVE_ENVS[Symbol(ce.size_A, ce.size_B, ce.crop, ce.padding)] = ce
end
function register_normxcorr2env(ne::Normxcorr2Env)
	NORMXCORR2_ENVS[Symbol(ne.size_A, ne.size_B)] = ne
end

function get_convolveenv(A, B; crop = :valid, padding = :none)
  if haskey(CONVOLVE_ENVS, Symbol(size(A), size(B), crop, padding)) return CONVOLVE_ENVS[Symbol(size(A), size(B), crop, padding)] end

  ce = ConvolveEnv(A, B, crop = crop, padding = padding)
  register_convolveenv(ce)

  return ce
end
function get_normxcorr2env(A, B)
  if haskey(NORMXCORR2_ENVS, Symbol(size(A), size(B))) return NORMXCORR2_ENVS[Symbol(size(A), size(B))] end

  ne = Normxcorr2Env(A, B)
  register_normxcorr2env(ne)

  return ne
end

function ConvolveEnv(A, B; factorable = Primes.primes(10), crop = :valid, padding = :none)
	if crop == :valid
		result_ranges = [min(a,b):max(a,b) for (a,b) in zip(size(A), size(B))]
	elseif crop == :same
	  	result_ranges = [(min(a,b)-div(min(a,b)-1, 2)):(max(a,b)+div(min(a,b)-1, 2)) for (a,b) in zip(size(A), size(B))]
	elseif crop == :full
		result_ranges = [1:(a+b-1) for (a,b) in zip(size(A), size(B))]
	end
    #common_size=tuple(map(+,size(A),size(B))...)
    common_size = tuple(map(maximum,result_ranges)...)
    common_size = (nextprod(factorable, common_size[1]), nextprod(factorable, common_size[2]))

    size_A = size(A)
    size_B = size(B)
    
    intermediate_A = zeros(Float64, common_size)
    intermediate_B = zeros(Float64, common_size)

    complex_intermediate_A = zeros(Complex{Float64}, div(common_size[1], 2)+1, common_size[2])
    complex_intermediate_B = zeros(Complex{Float64}, div(common_size[1], 2)+1, common_size[2])

    rfft_plan = plan_rfft(intermediate_A)
    irfft_plan = plan_irfft(complex_intermediate_A, common_size[1])


    result_size = tuple(map(length, result_ranges)...)
    result = zeros(Float64, result_size)

    return ConvolveEnv(
  	size_A::NTuple{2, Int64},
  	size_B::NTuple{2, Int64},

  	crop::Symbol,
	padding::Symbol,

  	intermediate_A::Array{Float64, 2},
  	intermediate_B::Array{Float64, 2},

  	complex_intermediate_A::Array{Complex{Float64}, 2},
  	complex_intermediate_B::Array{Complex{Float64}, 2},

  	rfft_plan,
  	irfft_plan,

  	result_ranges::Array{UnitRange{Int64},1},
  	result::Array{Float64, 2}
    )
end

function Normxcorr2Env(A, B)
    (n1,n2)=size(A);
    (m1,m2)=size(B);
    conv_dt = zeros(Float64, n1, n2)
    conv_imgpad = zeros(Float64, m1 + 1, m2 + 1)
    conv_sum = zeros(Float64, m1 + 1, m2 + 1)
    conv_sum2 = zeros(Float64, m1 + 1, m2 + 1)
    local_sum = zeros(Float64, m1 - n1 + 1, m2 - n2 + 1)
    local_sum2 = zeros(Float64, m1 - n1 + 1, m2 - n2 + 1)

    return Normxcorr2Env(
    	(n1,n2),
    	(m1,m2),
	conv_dt,
	conv_imgpad,
	conv_sum,
	conv_sum2,
	local_sum,
	local_sum2
    )
end

function clean!(ce::ConvolveEnv)
  @inbounds begin
	ce.intermediate_A[:] = Float64(0)
	ce.intermediate_B[:] = Float64(0)
	ce.complex_intermediate_A[:] = Complex{Float64}(0)
	ce.complex_intermediate_B[:] = Complex{Float64}(0)
  end
end

#=
immutable ConvolveEnv

CONV_RESULT = Array{Float64, 2}(10,10);
COMPLEX_CONV_INTERMEDIATE_A = Array{Complex{Float64}, 2}(6,10);
COMPLEX_CONV_INTERMEDIATE_B = Array{Complex{Float64}, 2}(6,10);

CONV_INTERMEDIATE_A = Array{Float64, 2}(10,10);
CONV_INTERMEDIATE_B = Array{Float64, 2}(10,10);

IMG_CROPPED_FACTORABLE = Array{Float64, 2}(10,10);
TEMPLATE_CROPPED_FACTORABLE = Array{Float64, 2}(10,10);

FFT_INPLACE_PLAN_A = plan_fft!(COMPLEX_CONV_INTERMEDIATE_A)
FFT_INPLACE_PLAN_B = plan_fft!(COMPLEX_CONV_INTERMEDIATE_B)
IFFT_INPLACE_PLAN_A = plan_ifft!(COMPLEX_CONV_INTERMEDIATE_A)

RFFT_PLAN = plan_rfft(CONV_INTERMEDIATE_A);
IRFFT_PLAN = plan_irfft(COMPLEX_CONV_INTERMEDIATE_A, 10);

CONV_DT = Array{Float64, 2}(0,0);
CONV_IMGPAD = Array{Float64, 2}(0,0);
CONV_SUM = Array{Float64, 2}(0,0);
CONV_SUM2 = Array{Float64, 2}(0,0);
LOCAL_SUM = Array{Float64, 2}(0,0);
LOCAL_SUM2 = Array{Float64, 2}(0,0);

end
=#
#=
function init_Convolve()

global CONV_RESULT = zeros(Float64, 10,10);
global COMPLEX_CONV_INTERMEDIATE_A = zeros(Complex{Float64}, 6,10);
global COMPLEX_CONV_INTERMEDIATE_B = zeros(Complex{Float64}, 6,10);

global CONV_INTERMEDIATE_A = zeros(Float64, 10,10);
global CONV_INTERMEDIATE_B = zeros(Float64, 10,10);

global IMG_CROPPED_FACTORABLE = zeros(Float64, 10,10);
global TEMPLATE_CROPPED_FACTORABLE = zeros(Float64, 10,10);

global FFT_INPLACE_PLAN_A = plan_fft!(COMPLEX_CONV_INTERMEDIATE_A)
global FFT_INPLACE_PLAN_B = plan_fft!(COMPLEX_CONV_INTERMEDIATE_B)
global IFFT_INPLACE_PLAN_A = plan_ifft!(COMPLEX_CONV_INTERMEDIATE_A)

global RFFT_PLAN = plan_rfft(CONV_INTERMEDIATE_A);
global IRFFT_PLAN = plan_irfft(COMPLEX_CONV_INTERMEDIATE_A, 10);

global CONV_DT = Array{Float64, 2}(0,0);
global CONV_IMGPAD = Array{Float64, 2}(0,0);
global CONV_SUM = Array{Float64, 2}(0,0);
global CONV_SUM2 = Array{Float64, 2}(0,0);
global LOCAL_SUM = Array{Float64, 2}(0,0);
global LOCAL_SUM2 = Array{Float64, 2}(0,0);

end

init_Convolve();
=#

function convolve{T, n}(A::SubArray{T, n},B::Array{T, n},dims)
    common_size=tuple(map(max,size(A),size(B))...)
    pA=zeros(Complex{T},common_size)
    pB=zeros(Complex{T},common_size)
    rangesA=[1:x for x in size(A)]
    rangesB=[1:x for x in size(B)]
    pA[rangesA...]=A
    pB[rangesB...]=B
    
    real(ifft!(fft!(pA,dims).*fft!(pB,dims),dims))
end

function convolve_ComplexFloat64(A,B,dims)
    common_size=tuple(map(max,size(A),size(B))...)
    pA=zeros(Complex{Float64},common_size)
    pB=zeros(Complex{Float64},common_size)
    rangesA=[1:x for x in size(A)]
    rangesB=[1:x for x in size(B)]
    pA[rangesA...]=A
    pB[rangesB...]=B
    
    real(ifft!(fft!(pA,dims).*fft!(pB,dims),dims))
end

# got rid of dims argument
# using real to complex fft and ifft
function convolve_Float64(A,B)
    common_size=tuple(map(max,size(A),size(B))...)
    pA=zeros(Float64,common_size)
    pB=zeros(Float64,common_size)
    rangesA=[1:x for x in size(A)]
    rangesB=[1:x for x in size(B)]
    pA[rangesA...]=A
    pB[rangesB...]=B
    irfft(rfft(pA).*rfft(pB),common_size[1])
end

function convolve_Float64_planned(A,B; flip = true, factorable = nothing, crop = :valid, padding = :none, ce = get_convolveenv(A, B; crop = crop, padding = padding))
	clean!(ce)

    	rangesA=[1:x for x in size(A)]
    	rangesB= flip ? [x:-1:1 for x in size(B)] : [1:x for x in size(B)] 
	if padding == :mean
	  @fastmath @inbounds ce.intermediate_A[:] = mean(A)
	end
    	@inbounds ce.intermediate_A[rangesA...]=A
    	@inbounds ce.intermediate_B[rangesB...]=B

    	@fastmath A_mul_B!(ce.complex_intermediate_A, ce.rfft_plan, ce.intermediate_A);
    	@fastmath A_mul_B!(ce.complex_intermediate_B, ce.rfft_plan, ce.intermediate_B);
    	@fastmath elwise_mul!(ce.complex_intermediate_A, ce.complex_intermediate_B);
    	@fastmath A_mul_B!(ce.intermediate_A, ce.irfft_plan, ce.complex_intermediate_A)

	@inbounds ce.result[:] = view(ce.intermediate_A, ce.result_ranges...)

	return ce.result
end

function elwise_complex_mul!(A, B)
      @simd for i in 1:length(A)
	@fastmath @inbounds A[i] = A[i] * B[i]
      end
      return A
end
function elwise_mul!(A, B)
      @simd for i in 1:length(A)
	@fastmath @inbounds A[i] = A[i] * B[i]
      end
      return A
end
function elwise_add!(A, B)
      @simd for i in 1:length(A)
	@fastmath @inbounds A[i] = A[i] + B[i]
      end
      return A
end
function elwise_sub!(A, B)
      @simd for i in 1:length(A)
	@fastmath @inbounds A[i] = A[i] - B[i]
      end
      return A
end
function elwise_div!(A, B)
      @simd for i in 1:length(A)
	@fastmath @inbounds A[i] = A[i] / B[i]
      end
      return A
end
function calculate_dt!(arr)
  	m = mean(arr)
	@simd for i in 1:length(arr)
	 @fastmath @inbounds arr[i] = arr[i] - m
	end
end

function make_real_and_copy!(C, A, offset)
  size_A = size(A);
      @simd for i in 1:length(C)
	sub = ind2sub(C,i);
	@fastmath @inbounds newsub = sub[1] + offset[1], sub[2] + offset[2];
	@fastmath @inbounds newind = sub2ind(size_A, newsub...)
	@fastmath @inbounds C[i] = real(A[newind])
      end
      return C
end

function convolve_ComplexFloat64!(C,A,B,ranges)
    common_size=tuple(map(max,size(A),size(B))...)
    pA=zeros(Complex{Float64},common_size)
    pB=zeros(Complex{Float64},common_size)
    rangesA=[1:x for x in size(A)]
    rangesB=[1:x for x in size(B)]
    pA[rangesA...]=A
    pB[rangesB...]=B
    dims = 1:2
      
      fft!(pA,dims);
      fft!(pB,dims);
      elwise_complex_mul!(pA, pB)
      ifft!(pA, dims);
     offset = (first(ranges[1])-1), (first(ranges[2])-1)
     make_real_and_copy!(C, pA, offset)
   end

#=
function valid_convolve(A,B)
    ranges=[min(a,b):max(a,b) for (a,b) in zip(size(A),size(B))]
    #convolve_Float64(A,B)[ranges...]
    convolve_Float64(A,B)[ranges...]
end
=#

function convolve_ComplexFloat64_prealloc_flipped(A,B,ranges)
    common_size=tuple(map(max,size(A),size(B))...)
    if size(COMPLEX_CONV_INTERMEDIATE_A) != common_size
      global COMPLEX_CONV_INTERMEDIATE_A=zeros(Complex{Float64},common_size)
      global FFT_INPLACE_PLAN_A = plan_fft!(COMPLEX_CONV_INTERMEDIATE_A)
      global IFFT_INPLACE_PLAN_A = plan_ifft!(COMPLEX_CONV_INTERMEDIATE_A)
    else
      COMPLEX_CONV_INTERMEDIATE_A[:] = 0;
    end
    if size(COMPLEX_CONV_INTERMEDIATE_B) != common_size
      global COMPLEX_CONV_INTERMEDIATE_B=zeros(Complex{Float64},common_size)
      global FFT_INPLACE_PLAN_B = plan_fft!(COMPLEX_CONV_INTERMEDIATE_B)
    else
      COMPLEX_CONV_INTERMEDIATE_B[:] = 0;
    end
    rangesA=[1:x for x in size(A)]
    rangesB=[x:-1:1 for x in size(B)]
    COMPLEX_CONV_INTERMEDIATE_A[rangesA...]=A
    COMPLEX_CONV_INTERMEDIATE_B[rangesB...]=B
    dims = 1:2
      
     # fft!(COMPLEX_CONV_INTERMEDIATE_A,dims);
      #fft!(COMPLEX_CONV_INTERMEDIATE_B,dims);
      FFT_INPLACE_PLAN_A * COMPLEX_CONV_INTERMEDIATE_A
      FFT_INPLACE_PLAN_B * COMPLEX_CONV_INTERMEDIATE_B
      elwise_complex_mul!(COMPLEX_CONV_INTERMEDIATE_A, COMPLEX_CONV_INTERMEDIATE_B)
      #ifft!(COMPLEX_CONV_INTERMEDIATE_A, dims);
      IFFT_INPLACE_PLAN_A * COMPLEX_CONV_INTERMEDIATE_A
     offset = (first(ranges[1])-1), (first(ranges[2])-1)
     make_real_and_copy!(CONV_RESULT, COMPLEX_CONV_INTERMEDIATE_A, offset)
   end



function valid_convolve(A,B; flip = true)
    return convolve_Float64_planned(A, B; crop = :valid, flip = flip)
end

function same_convolve(A,B; flip = true)
    return convolve_Float64_planned(A, B; crop = :same, flip = flip)
end

function cumsum2{T,ndim}(A::Array{T,ndim})
    # cumulative sum in two dimensions
    if ndim!=2
        throw(ArgumentError("input must be two-dimensional array"))
    end
    B = similar(A);
    # first row and column are 1D cumulative sums
    B[:,1]=cumsum(A[:,1],1);
    B[1,:]=cumsum(A[1,:],2);      # B[1,1] is redundantly computed twice
    # compute rest of matrix from recursion
    (m,n)=size(A);
    if m>1 && n>1
        for j in 2:n
		for i in 2:m
            	@fastmath @inbounds B[i,j]=A[i,j]+B[i-1,j]+B[i,j-1]-B[i-1,j-1]
		end
        end
    end
    B
end

function optimize_normxcorr2(img)
    p=plan_rfft(img,flags=FFTW.MEASURE)
    q=plan_irfft(rfft(img),flags=FFTW.MEASURE,size(img,1))
    normxcorr2(img, img);
    return Void;
end

# extension:
# Params_session.jl: optimize_all_cores(params::Params)
#
function optimize_all_cores(img_d::Int64)
    img = rand(img_d, img_d)
	@sync begin
	  for p in procs()
	    @async begin
		remotecall_fetch(p, optimize_normxcorr2, img);
		end
	    end
	  end
    return Void;
end

function cumsum12!(A::Array{Float64,2}, conv_sum::Array{Float64,2}, conv_sum2::Array{Float64,2})
function calculate_cumsums!(sum, A)
    (m,n)=size(A);
    if m>1 && n>1
        for j in 2:n
		for i in 2:m
        	@fastmath @inbounds    sum[i,j]=A[i,j]+sum[i-1,j]+sum[i,j-1]-sum[i-1,j-1]
		end
        end
    end
end

    # cumulative sum in two dimensions
    # first row and column are 1D cumulative sums
    conv_sum[:,1]=cumsum(A[:,1],1);
    conv_sum[1,:]=cumsum(A[1,:],2);      # B[1,1] is redundantly computed twice
    # compute rest of matrix from recursion
    calculate_cumsums!(conv_sum, A);
	elwise_mul!(A, A)
    conv_sum2[:,1]=cumsum(A[:,1],1);
    conv_sum2[1,:]=cumsum(A[1,:],2);      # conv_sum[1,1] is redundantly computed twice
    # compute rest of matrix from recursion
    calculate_cumsums!(conv_sum2, A);
end


function calculate_local_sums(local_sum::Array{Float64,2}, conv_sum::Array{Float64,2}, local_sum2::Array{Float64,2}, conv_sum2::Array{Float64,2}, L1, L2, n1, n2)
  j0 = first(L2) - 1
  i0 = first(L1) - 1
  i_size, j_size = size(conv_sum)
  for j in L2
	  @simd for i in L1
    j_ind = j - j0
    i_ind = i - i0
    @fastmath @inbounds local_sum[i_ind,j_ind] = conv_sum[i, j] - conv_sum[(i-n1), j] - conv_sum[i, (j-n2)] + conv_sum[(i-n1), (j-n2)]
    end
   end 
  for j in L2
	  @simd for i in L1
    j_ind = j - j0
    i_ind = i - i0
	@fastmath @inbounds local_sum2[i_ind,j_ind] = conv_sum2[i, j] - conv_sum2[(i-n1), j] - conv_sum2[i, (j-n2)] + conv_sum2[(i-n1), (j-n2)]
	end
    end 
end

function calculate_local_variance!(sum2::Array{Float64,2}, sum::Array{Float64,2}, template::Array{Float64,2})
	den = prod(size(template))
	@simd for i in 1:length(sum2)
	  @fastmath @inbounds sum2[i] = sum2[i] - sum[i] * sum[i] / den;
	end
	  return sum2
end
function calculate_denominator!(localvariance::Array{Float64,2}, templatevariance::Float64)
	for i in 1:length(localvariance)
	 @fastmath @inbounds localvariance[i] = sqrt(localvariance[i] * templatevariance)
	end
	  return localvariance
end
function normxcorr2_preallocated(template,img; shape = "valid", highpass_sigma = 0)
    # "normalized cross correlation": slide template across img,
    # compute Pearson correlation coefficient for each template location
    # result has same size as MATLAB-style 'valid' convolution
    # efficient algorithm of J. P. Lewis
    # http://scribblethink.org/Work/nvisionInterface/nip.html
    #
    # return NaN for image patch or template with zero variance
    # Pearson correlation coefficient is undefined in this case
    # 
    # need some argument checking
    # e.g. sizes of template should be less than those of img
    # this works for arrays.  extend to Image defined in Holy's package?

    (n1,n2)=size(template);

    if shape == "full"
    (n1,n2)=size(template);
    (m01,m02)=size(img);

    k1 = size(template, 1) * 2 - 2 + size(img, 1);
    k2 = size(template, 2) * 2 - 2 + size(img, 2);
    img_new = fill(mean(img), k1, k2);
    img_new[n1:k1-n1+1, n2:k2-n2+1] = img;
    img = img_new;
    end

    ne = get_normxcorr2env(template, img);
    @inbounds ne.conv_dt[:] = template
    @fastmath @inbounds calculate_dt!(ne.conv_dt)
    @fastmath @inbounds numerator::Array{Float64} = valid_convolve(img, ne.conv_dt; flip = true)
    @fastmath @inbounds templatevariance::Float64 = sum(elwise_mul!(ne.conv_dt, ne.conv_dt))

    ##### local statistics of img
    # zero pad image in first row and column
    # so that cumulative sums will have zeros in the same place
    (m1,m2)=size(img);
    @inbounds ne.conv_imgpad[1, 1:end] = 0
    @inbounds ne.conv_imgpad[1:end, 1] = 0
    @fastmath @inbounds ne.conv_imgpad[2:end,2:end]=img;
    # define four combinations of Small and Large ranges
    if templatevariance==0
        return zeros(m1-n1+1,m2-n2+1)*NaN
    end
    @fastmath @inbounds LL=UnitRange{Int64}[1+(n1:m1),1+(n2:m2)]
#    @fastmath @inbounds SL=LL-[n1;0]; LS=LL-[0;n2]
#    @fastmath @inbounds SS=LL-[n1;n2]
    # sum of img and its square in template-sized neighborhoods
    @fastmath @inbounds cumsum12!(ne.conv_imgpad, ne.conv_sum, ne.conv_sum2)


	@fastmath @inbounds calculate_local_sums(ne.local_sum, ne.conv_sum, ne.local_sum2, ne.conv_sum2, LL[1], LL[2], n1, n2);
	# not actually local variance, but local variance * prod(size(img))
	@fastmath @inbounds localvariance = calculate_local_variance!(ne.local_sum2, ne.local_sum, template)

    # localvariance is zero for image patches that are constant
    # leading to undefined Pearson correlation coefficient
    # should only be negative due to roundoff error
    eps_scaled = Float64(eps_large * prod(size(template)))
    for i in 1:length(localvariance) 
    	@inbounds if localvariance[i] <= eps_scaled
	@fastmath @inbounds localvariance[i] = Alembic.eps 
        @fastmath @inbounds numerator[i] = 0.0
      end
    end
    @fastmath @inbounds denominator=calculate_denominator!(localvariance, templatevariance)
    @simd for i in 1:length(denominator) @inbounds if denominator[i] <= 0 
    	@inbounds denominator[i] = eps; 
	@inbounds numerator[i] = 0 
    	end end
    #@fastmath @inbounds numerator[denominator.<=0] = 0
    #@fastmath @inbounds denominator[denominator.<=0] = eps

    @fastmath @inbounds xc = elwise_div!(numerator, denominator);

#    @fastmath @inbounds xcd = Images.imfilter_gaussian(xc, sigma)
	return xc
end
