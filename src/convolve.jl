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

function valid_convolve(A,B)
    ranges=[min(a,b):max(a,b) for (a,b) in zip(size(A),size(B))]
    convolve_Float64(A,B)[ranges...]
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
    return Void;
end

# extension:
# Params_session.jl: optimize_all_cores(params::Params)
#
function optimize_all_cores(img_d::Int64)
    img = rand(img_d, img_d)
	@sync begin
	  for p in 1:num_procs
	    @async begin
		remotecall_fetch(p, optimize_normxcorr2, img);
		end
	    end
	  end
    return Void;
end
    
function normxcorr2(template,img)
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

    # sufficient to subtract mean from just one variable
    dt=template-mean(template)
    templatevariance=sum(dt.^2)
    numerator=valid_convolve(img,dt[end:-1:1,end:-1:1])
    
    ##### local statistics of img
    # zero pad image in first row and column
    # so that cumulative sums will have zeros in the same place
    (m1,m2)=size(img);
    imgpad=zeros(m1+1,m2+1); imgpad[2:end,2:end]=img;
    # define four combinations of Small and Large ranges
    (n1,n2)=size(template);
    if templatevariance==0
        return zeros(m1-n1+1,m2-n2+1)*NaN
    end
    LL=Any[1+(n1:m1),1+(n2:m2)]
    SL=LL-[n1;0]; LS=LL-[0;n2]
    SS=LL-[n1;n2]
    # sum of img and its square in template-sized neighborhoods
    s=cumsum2(imgpad)
    localsum=s[LL...]-s[SL...]-s[LS...]+s[SS...]
    s2=cumsum2(imgpad.^2)
    localsum2=s2[LL...]-s2[SL...]-s2[LS...]+s2[SS...]
    localvariance=localsum2-localsum.^2/prod(size(template))
    # localvariance is zero for image patches that are constant
    # leading to undefined Pearson correlation coefficient
    # should only be negative due to roundoff error
    localvariance[localvariance.<=0] *= NaN
    denominator=sqrt(localvariance*templatevariance)

    numerator./denominator
end
    

