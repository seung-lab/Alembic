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

