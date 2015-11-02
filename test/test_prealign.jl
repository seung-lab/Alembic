function test()
  getimage(path) = convert(Array{Float64, 2}, convert(Array, Images.load(path)))
  moving_section = getimage("../sections/S2-W001_fixed_section4_0.175.tif")
  #fixed_section = getimage("./sections/S2-W001_fixed_section4_0.175.tif")
  fixed_section = getimage("../sections/S2-W001_fixed_section5_0.175.tif")
  println(size(moving_section))
  println(size(fixed_section))
  trans, moving_points, fixed_points, res1, res2 = affine_align_sections(moving_section, fixed_section, PARAMS_PREALIGNMENT; return_points=true)
  moving_points = points_to_3xN_matrix(moving_points)
  fixed_points = points_to_3xN_matrix(fixed_points)

  println(trans)
  #out_img, offset = imwarp(fixed_section, trans)
  #imwrite(moving_section, joinpath(".","test_outputs", string("moving_section", ".tif")))
  #imwrite(fixed_section, joinpath(".","test_outputs", string("fixed_section_", offset[1], "_", offset[2], ".tif")))
  #imwrite(out_img, joinpath(".","test_outputs", string("warped_", offset[1], "_", offset[2], ".tif")))
  p22 = fixed_points + res2
  p11 = moving_points + res1
  moving_points = moving_points[1:2,:]
  fixed_points = fixed_points[1:2,:]
  draw_vectors(fixed_section, vcat(fixed_points[1:2,:], p22[1:2,:]))
  #ccp.write_image_from_points(moving_points[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
  draw_points(moving_section, moving_points)
  draw_points(fixed_section, fixed_points)
  #draw_points(out_img, moving_points)
end


function test2()
  getimage(path) = convert(Array{Ufixed8, 2}, convert(Array, Images.load(path)))
  getimage = get_ufixed8_image
  #moving_section = getimage("../output_images/(1,1)_montage.tif")
  #fixed_section = getimage("../output_images/(1,2)_montage.tif")
  m = "/usr/people/smu/seungmount/research/Julimaps/datasets/piriform/2_montaged/1,77_montaged.tif"
  f = "/usr/people/smu/seungmount/research/Julimaps/datasets/piriform/2_montaged/1,78_montaged.tif"
  #m = "/usr/people/smu/seungmount/research/Julimaps/datasets/piriform/2_montaged/1,86_montaged.tif"
  #f = "/usr/people/smu/seungmount/research/Julimaps/datasets/piriform/2_montaged/1,87_montaged.tif"
  moving_section = getimage(m)
  fixed_section = getimage(f)
  println(size(moving_section))
  println(size(fixed_section))
  trans, moving_points, fixed_points, res1, res2 = affine_align_images(moving_section, fixed_section, PARAMS_PREALIGNMENT; return_points=true)
  res1 = res1.'
  res2 = res2.'
  moving_points = points_to_3xN_matrix(moving_points)
  fixed_points = points_to_3xN_matrix(fixed_points)

  println(trans)
  downsample = 4
  moving_section = moving_section[1:downsample:end, 1:downsample:end]
  fixed_section = fixed_section[1:downsample:end, 1:downsample:end]

  p11 = moving_points + res1
  p22 = fixed_points + res2

  p11 = p11[1:2,:]
  p22 = p22[1:2,:]
  moving_points = moving_points[1:2,:]
  fixed_points = fixed_points[1:2,:]

  p11 = p11/downsample
  p22 = p22/downsample
  moving_points = moving_points/downsample
  fixed_points = fixed_points/downsample
  
  trans = adjust_affine_for_scaling(trans, downsample)
  out_img, offset = imwarp(fixed_section, inv(trans))
  println(offset)
  fused, fused_offset = imfuse(moving_section, [0,0], out_img, offset)
  println(fused_offset)

  block_radius = 150
  scalebar = [1; 1; 2*block_radius/downsample; 2*block_radius/downsample]
  draw_vectors(make_isotropic(moving_section), hcat(vcat(moving_points, p11), scalebar))
  draw_vectors(make_isotropic(fixed_section), hcat(vcat(fixed_points, p22), scalebar+100))
  draw_vectors(make_isotropic(fused), hcat(vcat(moving_points.-fused_offset, p11.-fused_offset), scalebar)) # todo: check offset
    #ccp.write_image_from_points(moving_points[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
  decomp_affine(trans)
end

function test3()
  moving_section = "../output_images/(1,1)_montage.tif"
  fixed_section = "../output_images/(1,2)_montage.tif"
  A, meshset = affine_align_sections(moving_section, fixed_section)
  return meshset
end

function test_tile_align()
  # test on tiles
  params = copy(PARAMS_PREALIGNMENT)
  params["scaling_factor"] = 1
  params["mesh_length"] = 5000
  params["block_size"] = 300
  params["search_r"] = 1200
  params["min_r"] = 0.1
  sec = 77
  for row in 1:1:4
    for col in 1:1:4
      index = (1, sec, row, col)
      moving_section = get_ufixed8_image(index)
      index = (1, sec+1, row, col)
      fixed_section = get_ufixed8_image(index)
      println(size(moving_section))
      println(size(fixed_section))
      trans, = affine_align_images(moving_section, fixed_section, PARAMS_PREALIGNMENT; return_points=true)
      decomp_affine(trans)
    end
  end
end