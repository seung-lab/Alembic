# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150805
#
# Convert edge-node incidence matrix into triangle array (triples of
# node indices defining all triangles in the mesh). The triangle array
# will be used when applying piecewise affine transforms through the
# PiecewiseAffineTransform package.

# Incidence refers to an edge-node incidence, where every row represents
# an edge, and every column a node. If a node participates in an edge, it
# is marked by a 1.

using Base.Test

function find_nonzero_indices(a)
  eachindex(a)'[a .> 0]
end

function incidence2pairs(D)
  # Pare incidence matrix down to Nx2 vector of node index pairs
  D = abs(D)
  mapslices(find_nonzero_indices, D, 2)
end

function incidence2dict(D)
  # Create dictionary of sets of nodes connected to indexing node
  D = abs(D)
  node_dict = Dict()

  for j = 1:size(D,2)
    i = 1
    while D[i,j] == 0
      i += 1
    end
    if !(i in keys(node_dict))
      node_dict[i] = Set{Int64}()
    end
    k = i+1
    while D[k,j] == 0
      k += 1
    end
    push!(node_dict[i], k)
  end
  return node_dict
end

function dict2triangles(node_dict)
  # Convert node dictionary to list of 1x3 node index vectors for triangles
  triangles = Array(Int64, 0, 3)
  for a in sort(collect(keys(node_dict)))
    setA = node_dict[a]
    for b in sort(collect(setA))
      if !(b in keys(node_dict))
        continue
      end
      setB = node_dict[b]
      setC = intersect(setA, setB)
      for c in sort(collect(setC))
        triangles = vcat(triangles, [a b c])
      end
    end
  end
  return triangles
end

function incidence2triangles(D)
  node_dict = incidence2dict(D)
  dict2triangles(node_dict)
end

function test_incidence2triangles()
  D = [1 -1 0 0 0;
     1 0 -1 0 0;
     0 1 -1 0 0;
     0 1 0 -1 0;
     0 1 0 0 -1;
     0 0 0 1 -1]';
  triangles = [1 2 3;
        2 4 5];
  tri = incidence2triangles(D)
  @test triangles == tri

  D = [1 -1 0 0 0;
     1 0 -1 0 0;
     1 0 0 -1 0;
     1 0 0 0 -1;
     0 1 -1 0 0;
     0 1 0 -1 0;
     0 1 0 0 -1;
     0 0 1 -1 0;
     0 0 1 0 -1;
     0 0 0 1 -1]';
  triangles = [1 2 3;
        1 2 4;
        1 2 5;
        1 3 4;
        1 3 5;
        1 4 5;
        2 3 4;
        2 3 5;
        2 4 5;
        3 4 5];
  tri = incidence2triangles(D)
  @test triangles == tri

  D = [1 -1 0 0 0;
     0 0 1 -1 0;
     0 0 1 0 -1;
     1 0 -1 0 0;
     1 0 0 -1 0;
     1 0 0 0 -1;
     0 1 -1 0 0;
     0 1 0 -1 0;
     0 1 0 0 -1;
     0 0 0 1 -1]';
  triangles = [1 2 3;
        1 2 4;
        1 2 5;
        1 3 4;
        1 3 5;
        1 4 5;
        2 3 4;
        2 3 5;
        2 4 5;
        3 4 5];
  tri = incidence2triangles(D)
  @test triangles == tri
end

function test_incidence2dict()
  D = [1 -1 0 1 0;
     0 0 1 -1 0]';
  node_dict = Dict(1 => Set(2), 3 => Set(4))
  nd = incidence2dict(D)
  @test node_dict == nd

  D = [1 -1 0 1 0;
     1 -1 1 1 0;
     1 0 1 1 0]';
  node_dict = Dict(1 => Set([2, 3]))
  nd = incidence2dict(D)
  @test node_dict == nd 
end
