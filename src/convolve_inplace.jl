global CONV_RESULT = Array{Float64, 2}(10,10);
global COMPLEX_CONV_INTERMEDIATE_A = Array{Complex{Float64}, 2}(6,10);
global COMPLEX_CONV_INTERMEDIATE_B = Array{Complex{Float64}, 2}(6,10);

global CONV_INTERMEDIATE_A = Array{Float64, 2}(10,10);
global CONV_INTERMEDIATE_B = Array{Float64, 2}(10,10);

global IMG_CROPPED_FACTORABLE = Array{Float64, 2}(10,10);
global TEMPLATE_CROPPED_FACTORABLE = Array{Float64, 2}(10,10);

global FFT_INPLACE_PLAN_A = plan_fft!(COMPLEX_CONV_INTERMEDIATE_A)
global FFT_INPLACE_PLAN_B = plan_fft!(COMPLEX_CONV_INTERMEDIATE_B)
global IFFT_INPLACE_PLAN_A = plan_ifft!(COMPLEX_CONV_INTERMEDIATE_A)

global RFFT_PLAN = plan_rfft(CONV_INTERMEDIATE_A);
global IRFFT_PLAN = plan_irfft(COMPLEX_CONV_INTERMEDIATE_A, 10);

global CONV_DT = Array{Float64, 2}();
global CONV_IMGPAD = Array{Float64, 2}();
global CONV_SUM = Array{Float64, 2}();
global CONV_SUM2 = Array{Float64, 2}();
global LOCAL_SUM = Array{Float64, 2}();
global LOCAL_SUM2 = Array{Float64, 2}();

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

function convolve_Float64_planned(A,B, ranges; factorable = nothing)
    common_size=tuple(map(max,size(A),size(B))...)
    #pad to the smallest multiple of the factorables larger than the common size in each dimensions for speed
    if factorable != nothing
	common_size = (nextprod(factorable, common_size[1]), nextprod(factorable, common_size[2]))
    end
 
    if size(CONV_INTERMEDIATE_A) != common_size
      global CONV_INTERMEDIATE_A=zeros(Float64,common_size)
      global CONV_INTERMEDIATE_B=zeros(Float64,common_size)
      global COMPLEX_CONV_INTERMEDIATE_A=zeros(Complex{Float64}, div(common_size[1], 2)+1, common_size[2])
      global COMPLEX_CONV_INTERMEDIATE_B=zeros(Complex{Float64}, div(common_size[1], 2)+1, common_size[2])
      global RFFT_PLAN = plan_rfft(CONV_INTERMEDIATE_A)
      global IRFFT_PLAN = plan_irfft(COMPLEX_CONV_INTERMEDIATE_A, common_size[1])
    else
      CONV_INTERMEDIATE_A[:] = 0;
      CONV_INTERMEDIATE_B[:] = 0;
      COMPLEX_CONV_INTERMEDIATE_A[:] = 0;
      COMPLEX_CONV_INTERMEDIATE_B[:] = 0;
    end

    rangesA=[1:x for x in size(A)]
    rangesB=[x:-1:1 for x in size(B)]
    @inbounds CONV_INTERMEDIATE_A[rangesA...]=A
    @inbounds CONV_INTERMEDIATE_B[rangesB...]=B

    @fastmath A_mul_B!(COMPLEX_CONV_INTERMEDIATE_A, RFFT_PLAN, CONV_INTERMEDIATE_A);
    @fastmath A_mul_B!(COMPLEX_CONV_INTERMEDIATE_B, RFFT_PLAN, CONV_INTERMEDIATE_B);
    @fastmath elwise_mul!(COMPLEX_CONV_INTERMEDIATE_A, COMPLEX_CONV_INTERMEDIATE_B);
    @fastmath A_mul_B!(CONV_INTERMEDIATE_A, IRFFT_PLAN, COMPLEX_CONV_INTERMEDIATE_A)
    
   # (IRFFT_PLAN * elwise_mul!(RFFT_PLAN * CONV_INTERMEDIATE_A, RFFT_PLAN * CONV_INTERMEDIATE_B))[ranges...]
    CONV_RESULT[:] = CONV_INTERMEDIATE_A[ranges...]
   # irfft(rfft(pA).*rfft(pB),common_size[1])
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


function valid_convolve(A,B)
    ranges=[min(a,b):max(a,b) for (a,b) in zip(size(A),size(B))]
    #convolve_Float64(A,B)[ranges...]
    convolve_Float64(A,B)[ranges...]
end

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



function valid_convolve_flipped(A,B; factorable = nothing)
    ranges=[min(a,b):max(a,b) for (a,b) in zip(size(A),size(B))]
    if size(CONV_RESULT) != (length(ranges[1]), length(ranges[2]))
	global CONV_RESULT = Array{Float64}(length(ranges[1]), length(ranges[2]));
      end
    #=convolve_ComplexFloat64_prealloc_flipped(A,B, ranges)=#
    convolve_Float64_planned(A, B, ranges; factorable = factorable)
    return CONV_RESULT
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
        for j in 2:n, i in 2:m
            B[i,j]=A[i,j]+B[i-1,j]+B[i,j-1]-B[i-1,j-1]
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

function cumsum12!{T,ndim}(A::Array{T,ndim})
function calculate_cumsums!(SUM, A)
    (m,n)=size(A);
    if m>1 && n>1
        for j in 2:n, i in 2:m
        @fastmath @inbounds    SUM[i,j]=A[i,j]+SUM[i-1,j]+SUM[i,j-1]-SUM[i-1,j-1]
        end
    end
end

    # cumulative sum in two dimensions
    if ndim!=2
        throw(ArgumentError("input must be two-dimensional array"))
    end
    if size(CONV_SUM) != size(A)
      global CONV_SUM = similar(A);
    end
    if size(CONV_SUM2) != size(A)
      global CONV_SUM2 = similar(A);
    end
    # first row and column are 1D cumulative sums
    CONV_SUM[:,1]=cumsum(A[:,1],1);
    CONV_SUM[1,:]=cumsum(A[1,:],2);      # B[1,1] is redundantly computed twice
    # compute rest of matrix from recursion
    calculate_cumsums!(CONV_SUM, A);
	elwise_mul!(A, A)
    CONV_SUM2[:,1]=cumsum(A[:,1],1);
    CONV_SUM2[1,:]=cumsum(A[1,:],2);      # CONV_SUM[1,1] is redundantly computed twice
    # compute rest of matrix from recursion
    calculate_cumsums!(CONV_SUM2, A);
end


function calculate_local_sums(LOCAL_SUM, CONV_SUM, LOCAL_SUM2, CONV_SUM2, L1, L2, n1, n2)
  j0 = first(L2) - 1
  i0 = first(L1) - 1
  i_size, j_size = size(CONV_SUM)
  for j in L2,  i in L1
    j_ind = j - j0
    i_ind = i - i0
    @fastmath @inbounds LOCAL_SUM[i_ind,j_ind] = CONV_SUM[i, j] - CONV_SUM[(i-n1), j] - CONV_SUM[i, (j-n2)] + CONV_SUM[(i-n1), (j-n2)]
   end 
  for j in L2,  i in L1
    j_ind = j - j0
    i_ind = i - i0
	@fastmath @inbounds LOCAL_SUM2[i_ind,j_ind] = CONV_SUM2[i, j] - CONV_SUM2[(i-n1), j] - CONV_SUM2[i, (j-n2)] + CONV_SUM2[(i-n1), (j-n2)]
    end 
end

function calculate_local_variance!(SUM2, SUM, template)
	den = prod(size(template))
	for i in 1:length(SUM2)
	  SUM2[i] = SUM2[i] - SUM[i] * SUM[i] / den;
	end
	  return SUM2
end
function calculate_denominator!(localvariance, templatevariance)
	for i in 1:length(localvariance)
	  localvariance[i] = sqrt(localvariance[i] * templatevariance)
	end
	  return localvariance
end
function normxcorr2_preallocated(template,img; shape = "valid")
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

    factors = primes(10)

    if shape == "full"
    (n1,n2)=size(template);
    (m01,m02)=size(img);

    k1 = size(template, 1) * 2 - 2 + size(img, 1);
    k2 = size(template, 2) * 2 - 2 + size(img, 2);
    img_new = fill(mean(img), k1, k2);
    img_new[n1:k1-n1+1, n2:k2-n2+1] = img;
    img = img_new;
    end
    # sufficient to subtract mean from just one variable
    if size(CONV_DT) != size(template)
      global CONV_DT = Array(Float64, size(template)...);
    end
    @inbounds CONV_DT[:] = template[:]
    @fastmath @inbounds calculate_dt!(CONV_DT)
    @fastmath @inbounds numerator=valid_convolve_flipped(img,CONV_DT; factorable = factors)
    @fastmath @inbounds templatevariance=sum(elwise_mul!(CONV_DT, CONV_DT))

    
    ##### local statistics of img
    # zero pad image in first row and column
    # so that cumulative sums will have zeros in the same place
    (m1,m2)=size(img);
    if size(CONV_IMGPAD) != (m1 + 1, m2 + 1)
      global CONV_IMGPAD = zeros(Float64, m1+1, m2+1)
    end
    @fastmath @inbounds CONV_IMGPAD[2:end,2:end]=img;
    # define four combinations of Small and Large ranges
    (n1,n2)=size(template);
    if templatevariance==0
        return zeros(m1-n1+1,m2-n2+1)*NaN
    end
    @fastmath @inbounds LL=UnitRange{Int64}[1+(n1:m1),1+(n2:m2)]
#    @fastmath @inbounds SL=LL-[n1;0]; LS=LL-[0;n2]
#    @fastmath @inbounds SS=LL-[n1;n2]
    # sum of img and its square in template-sized neighborhoods
    @fastmath @inbounds cumsum12!(CONV_IMGPAD)
#    @fastmath @inbounds s=CONV_SUM
    if size(LOCAL_SUM) != (m1-n1+1, m2-n2+1)
      global LOCAL_SUM = Array(Float64, m1-n1+1, m2-n2+1);
      global LOCAL_SUM2 = Array(Float64, m1-n1+1, m2-n2+1);
    end

	@fastmath @inbounds calculate_local_sums(LOCAL_SUM, CONV_SUM, LOCAL_SUM2, CONV_SUM2, LL[1], LL[2], n1, n2);
	@fastmath @inbounds localvariance = calculate_local_variance!(LOCAL_SUM2, LOCAL_SUM, template)

    # localvariance is zero for image patches that are constant
    # leading to undefined Pearson correlation coefficient
    # should only be negative due to roundoff error
    #localvariance[localvariance.<=0] *= NaN
    @fastmath @inbounds denominator=calculate_denominator!(localvariance, templatevariance)
    @fastmath @inbounds numerator[denominator.<=0] = 0
    @fastmath @inbounds denominator[denominator.<=0] = eps

    @fastmath @inbounds return elwise_div!(numerator, denominator);
end
