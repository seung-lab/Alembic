"""
`AFFINE_ALIGN_IMAGES`

Returns the affine transformation tform, such that p_in_B = p_in_A * tform, where p is 
point coordinates in row vector.
"""
function affine_align_images(moving_img::Array{}, fixed_img::Array{}, 
                                params=PARAMS_PREALIGNMENT; return_points=false)
  # points here are in column vector convention
  if params["scaling_factor"] != 1.0
    println("Finding points at ", params["scaling_factor"], "x")
    s = [params["scaling_factor"] 0 0; 0 params["scaling_factor"] 0; 0 0 1]
    moving_img, moving_offset = imwarp(moving_img, s)
    fixed_img, fixed_offset = imwarp(fixed_img, s)
  end
  points = generate_match_points(moving_img, fixed_img, params)
  moving_pointslist, fixed_pointslist = get_block_matches(moving_img, fixed_img, points, params)

  n_matches = length(moving_pointslist)
  println("n_matches: ", n_matches)

  moving_points = points_to_Nx3_matrix(moving_pointslist)
  fixed_points = points_to_Nx3_matrix(fixed_pointslist)

  tform = find_affine(moving_points, fixed_points)
  #tform = find_rigid(moving_points, fixed_points)
  residualIn2 = moving_points*tform - fixed_points
  rmsIn2 = mean( sum(residualIn2.^2, 1) )^0.5
  fixed_pointsin1 = fixed_points*inv(tform)
  residualIn1 = fixed_pointsin1 - moving_points
  rmsIn1 = mean( sum(residualIn1.^2, 1) )^0.5

  tolerance_ratio = 0.001
  rmsThres = tolerance_ratio * mean([size(moving_img)..., size(fixed_img)...])
  rmsTotal = mean([rmsIn1, rmsIn2].^2)^0.5

  if n_matches < 0.5 * length(points)
    println("WARNING [affine_align_sections]: # of matches is small. n_matches: ", n_matches)
  end
  if  rmsTotal > rmsThres
    println("WARNING [affine_align_sections]: high residual. RMS error: ", rmsTotal)
  end

  tform = adjust_affine_for_scaling(tform, params["scaling_factor"])

  if return_points
    return tform, moving_pointslist, fixed_pointslist, residualIn1, residualIn2, rmsIn1, rmsIn2, rmsTotal
  else
    return tform
  end
end

"""
Calculate coordinates of 'mesh' points to be blockmatched
"""
function generate_match_points(moving_img::Array{}, fixed_img::Array{}, 
                                                  params=PARAMS_PREALIGNMENT)
  border_ratio = 0.1
  mesh_length = params["mesh_length"] * params["scaling_factor"]
  grid_size = minimum(floor(Int, collect(size(moving_img))/mesh_length))
  block_radius = params["block_size"] * params["scaling_factor"]
  search_radius = params["search_r"] * params["scaling_factor"]
  overlap = [min(size(moving_img), size(fixed_img))...]
  #overlap = size(moving_img)

  border = round(Int, border_ratio * minimum(overlap))
  combined_radius = search_radius+block_radius
  border = border > combined_radius ? border : combined_radius
  println("border excluded: ", border)

  step = floor(Int, (overlap - 2*border - 1) ./ (grid_size))
  points = [[x; y] for x = 1+border:step[1]:overlap[1]-border, 
                                y = 1+border:step[2]:overlap[2]-border]
  points = points[:]
  println("n points: ", size(points))
  return points
end

"""
Cycle through points list, pass images from cross-correlation, and store results
"""
function get_block_matches(moving_img::Array{}, fixed_img::Array{}, points, 
                          params=PARAMS_PREALIGNMENT)
  moving_points = []
  fixed_points = []
  block_radius = floor(Int64, params["block_size"] * params["scaling_factor"])
  search_radius = floor(Int64, params["search_r"] * params["scaling_factor"])
  for pt = points
    offset, r, xc = block_match_at_point(moving_img, pt, fixed_img, pt, 
                                            block_radius, search_radius)
    println(pt, offset, r)
    if r >= params["min_r"]
      push!(moving_points, collect(pt))
      push!(fixed_points, collect(pt+offset))
    end
  end
  return moving_points, fixed_points
end

"""
Scale an affine transform computed on scaled images
"""
function adjust_affine_for_scaling(tform, scale)
  s = [scale 0 0; 0 scale 0; 0 0 1]
  return s * tform * inv(s)
end

"""
`FIND_AFFINE` - Calculate left-hand affine transform from two Nx3 point sets

Inputs:
    2D arrays with each column being a point in homogeneous coords.
Returns:
    Affine transformation that produces moving_points*tform = fixed_points


Affine transform of a global position:

* homogeneous coordinates [x, y, 1]  [ax + by + c, dx + ey + f, 1]
* or equivalently [x, y, 1]  [x, y, 1] * tform

        where `tform` = [a d 0;  
                         b e 0;  
                         c f 1]
"""
function find_affine(moving_points, fixed_points)
  return moving_points \ fixed_points
end

function block_match_at_point(A, pointInA, B, pointInB, block_radius, search_r)
# A is the template
# Returns:
# offset: offset from pointInB to pointInA's actual match in B.
# r_max:  correlation value at the match point
# xc:   raw cross-correlation map

  b = block_radius
  B_radius = block_radius + search_r

  A_lower = pointInA - b
  A_upper = pointInA + b
  B_lower = pointInB - B_radius
  B_upper = pointInB + B_radius

  xc = normxcorr2(sub(A, [l:u for (l,u) in zip(A_lower, A_upper)]...),
          sub(B, [l:u for (l,u) in zip(B_lower, B_upper)]...))

  r_max, ind = findmax(xc); 
  i_max, j_max = ind2sub(xc, ind)

  return [i_max-1-search_r; j_max-1-search_r], r_max, xc
  
end

function recompute_transform(meshset::MeshSet; method = "affine")
  params = meshset.params
  moving_pointslist = meshset.meshes[1].nodes
  fixed_pointslist = meshset.matches[1].dst_points
  moving_points = points_to_Nx3_matrix(moving_pointslist)
  fixed_points = points_to_Nx3_matrix(fixed_pointslist)
  if method == "affine"
    tform = find_affine(moving_points, fixed_points)
  elseif method == "rigid"
    tform = find_rigid(moving_points, fixed_points)
  end
  tform = adjust_affine_for_scaling(tform, params["scaling_factor"])
  return tform
end

function affine_align_sections(moving_img_filename::String, 
                                                fixed_img_filename::String, 
                                                params=PARAMS_PREALIGNMENT)
  meshset = MeshSet(params)
  meshset.N = 2
  moving_mesh = Mesh(moving_img_filename)
  fixed_mesh = Mesh(fixed_img_filename)
  meshset.meshes = [fixed_mesh, moving_mesh]

  fixed_img = get_ufixed8_image(meshset.meshes[1])
  moving_img = get_ufixed8_image(meshset.meshes[2])

  offsets = parse_offsets(prealigned_offsets_path)
  if is_aligned(meshset.meshes[1].index)
    offsets = ALIGNED_OFFSETS
  end
  idx = findfirst(offsets[:,1], fixed_img_filename)
  fixed_offset = offsets[idx, 3:4]
  # @time moving_img = rescopeimage(moving_img, fixed_offset, GLOBAL_BB)

  # Align
  tform, moving_points, fixed_points, residualIn1, residualIn2, rmsIn1, rmsIn2, rmsTotal = 
    affine_align_images(moving_img, fixed_img, params; return_points=true)
  println(tform)

  # Update matches object in meshset
  meshset.meshes[2].nodes = moving_points
  src_points_indices = 1:length(moving_points)
  matches = Match(meshset.meshes[2].index, meshset.meshes[1].index, 
    length(moving_points), src_points_indices, fixed_points, [], [], [])
  meshset.matches = [matches]
  println("Saving JLD for ", moving_img_filename, 
                                  " to ", fixed_img_filename)
  save(meshset)  

  warped_index = meshset.meshes[2].index
  
  # Write overlay thumbnail for review
  scale = 0.125
  s = [scale 0 0; 0 scale 0; 0 0 1]
  @time imgA, A_offset = imwarp(fixed_img, s)
  @time imgB, B_offset = imwarp(moving_img, tform*s)

  pt_scaling = meshset.params["scaling_factor"]
  src_nodes = hcat(fixed_points/pt_scaling...)[1:2, :]
  dst_nodes = (points_to_Nx3_matrix(moving_points/pt_scaling)*tform)[:, 1:2]'
  # src_offset = [GLOBAL_BB.i, GLOBAL_BB.j]/pt_scaling
  # dst_offset = [GLOBAL_BB.i, GLOBAL_BB.j]/pt_scaling
  # src_nodes = (src_nodes .- src_offset)*scale
  # dst_nodes = (dst_nodes .- dst_offset)*scale
  src_nodes *= scale
  dst_nodes *= scale

  O, O_bb = imfuse(imgA, A_offset, imgB, B_offset)
  imgc, img2 = view(O, pixelspacing=[1,1])
  vectors = [src_nodes; dst_nodes]
  println(size(vectors))
  an_pts, an_vectors = draw_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1), 10)
  c = draw_indices(imgc, img2, src_nodes, 48.0, [-40,-40])

  thumbnail_fn = string(join(warped_index[1:2], ","), "_prealigned_thumbnail.png")
  println("Writing ", thumbnail_fn)
  Cairo.write_to_png(imgc.c.back, joinpath(PREALIGNED_DIR, "review", thumbnail_fn))
  destroy(toplevel(imgc))

  # Render image
  warped_fn = string(join(warped_index[1:2], ","), "_prealigned.tif")
  println("Rendering ", warped_fn)
  @time warped_img, warped_offset = imwarp(moving_img, tform)
  println("Writing ", warped_fn)
  @time imwrite(warped_img, joinpath(PREALIGNED_DIR, warped_fn))

  log_path = joinpath(PREALIGNED_DIR, "prealigned_offsets.txt")  
  update_offset_log!(log_path, warped_fn, warped_offset, size(warped_img))
end

function compute_propogated_transform(index::Index; method = "affine")
  index = (index[1:2]..., 0, 0)
  filenames = filter(x -> x[end-2:end] == "jld", readdir(PREALIGNED_DIR))
  indices = [(parse_name(x), x) for x in filenames]
  sort!(indices)
  if indices[1][1] == zeros(Int, 4)
    error("Could not parse JLD filename to index: ", indices[1][2])
  end
  for k in 2:length(indices)
    if indices[k][1] > index
      break
    end
    if !is_adjacent(indices[k-1][1], indices[k][1], true)
      error("Missing section between ", indices[k-1][1], " and ", indices[k][1])
    end
  end
  T = diagm(ones(3))
  for k in 1:length(indices)
    ind, fn = indices[k]
    if ind > index
      break
    end
    meshset = load(joinpath(PREALIGNED_DIR, fn))
    A = recompute_transform(meshset; method = method)
    # row vector homogeneous point convention
    T *= A
  end
  return T
end

function prealignment()
    @time prealign_directory()
    @time render_prealignment_for_directory()
end
