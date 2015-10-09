
"""
`find_rigid` - Calculate left-hand rigid transform from two N x d point sets

Inputs:
    points1, points2:
        2D arrays with each row being coords of a point of dimension d-1 (homogeneous coords) or d.
    homogeneous_coord = true:
        Whether the inputs are in homogeneous coords.
Returns:
    T: 
        Rigid transformation T that produces  points2 = points1 * T
        (if points1 and points2 are in homogeneous coords)
"""
function find_rigid(points1, points2; homogeneous_coord = true)
  n, dim = size(points1)
  p1_bar = mean(points1, 1)
  p2_bar = mean(points2, 1)
  p1 = points1 .- p1_bar
  p2 = points2 .- p2_bar
  C = p2.' * p1 / n
  U,s,V = svd(C)
  # Rotation R in least squares sense:  points2 - p2_bar = (points1 - p1_bar) * R
  R = ( U * diagm(vcat(ones(dim-1), det(U*V.'))) * V.' ).'
  trans = p2_bar - p1_bar * R
  if homogeneous_coord
    R[end, 1:end-1] = trans[1:end-1]
    return R
  else
    T = [R zeros(dim); trans 1]
    return T
  end
end

"""
moving_pts * tform = fixed_pts
"""
function find_affine(moving_pts, fixed_pts)
  return moving_pts \ fixed_pts
end