using NearestNeighbors
using Distances

#Given features (e.g. timestamp and mean intensity) and values (e.g. mean intensity),
#compute bias corrected values by subtracting from each point the 
#mean of the values of its nearest neighbors in feature space
#The distance_threshold argument sets a cutoff for which neighbours
#to consider.

#features should have shape (nf, np), where nf is the
#number of features and np is the number of points
function f{T}(features::Array{T,2},values::Array,distance_threshold)
	n = size(A,2)
	kdtree=NNTree(features,NNTree.Chebyshev)

	target_values = zeros(size(values))
	for i in 1:n
		nns = map(x->x[1],filter(x->x[2]<distance_threshold,knn(kdtree, A[:,i], 30)))
		target_values[i]=values[i]-mean(values[nns])
	end
	return target_values
end
