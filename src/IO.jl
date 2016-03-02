function save(path::String, img::Array)
      f = h5open(path, "w")
      @time f["img", "chunk", (1000,1000)] = img
      close(f)
    end

# extensions:
# Mesh.jl	get_image(mesh::Mesh)
# filesystem.jl	get_image() 
function get_image(path::String, dtype = UInt8)
	ext = splitext(path)[2];
  	if ext == ".tif"
  		img = data(FileIO.load(path))
  		img = img[:, :, 1]'
   		#img.properties["timedim"] = 0
  		return convert(Array{dtype, 2}, round(convert(Array, img)*255))
	elseif ext == ".h5"
 		return convert(Array{dtype, 2}, h5read(path, "img"))
	end
end

function get_uint8_image(path::String)
  img = Images.load(path)
  return reinterpret(UInt8, data(img)[:,:,1])'
end

function get_image(index, dtype = UInt8)
  return get_image(get_path(index), dtype)
end

function ufixed8_to_uint8(img)
  reinterpret(UInt8, -img)
end

function get_h5_path(index::Index)
  return string(get_path(index)[1:end-4], ".h5")
end

function get_h5_slice(path::String, slice)
  # return reinterpret(Ufixed8, h5read(path, "img", slice))
  return h5read(path, "img", slice)
end

function get_h5_image(path::String)
  return h5read(path, "img")
end

# extensions:
# Mesh.jl get_float_image(mesh::Mesh)
function get_float_image(path::String)
  img = Images.load(path)
  img.properties["timedim"] = 0
  return convert(Array{Float64, 2}, convert(Array, img[:,:,1]))
end

function get_ufixed8_image(path::String)
  return convert(Array{Ufixed8}, data(Images.load(path))[:,:,1])'
end

function load_affine(path::String)
  affinePath = joinpath(AFFINE_DIR, string(path, ".csv"))
  return readcsv(path)
end

function parse_rough_align(info_path::String)
  file = readdlm(info_path)
  session = cell(size(file, 1), 4); # name, index, dx, dy
  for i in 1:size(file, 1)
  m = match(r"(Tile\S+).tif", file[i, 1])
  session[i, 1] = m[1]
  session[i, 2] = m[1]
  session[i, 3] = parse_name(array[i, 1])
  session[i, 4] = file[i, 2]
  session[i, 5] = file[i, 3]
  end
  return session
end

# extensions:
# MeshSet.jl load_section_images(Ms::MeshSet)
function load_section_images(session, section_num)
  indices = find(i -> session[i,2][2] == section_num, 1:size(session, 1))
  max_tile_size = 0
  num_tiles = length(indices)
  paths = Array{String, 1}(num_tiles)
  for i in 1:num_tiles
    name = session[i, 1]
    paths[i] = get_path(name)
    image = get_image(paths[i])
    max_size = max(size(image, 1), size(image, 2))
    if max_tile_size < max_size max_tile_size = max_size; end
  end
  imageArray = SharedArray(UInt8, max_tile_size, max_tile_size, num_tiles)

  for k in 0:num_procs:num_tiles
    @sync @parallel for l in 1:num_procs
    i = k+l
    if i > num_tiles return; end
    image = get_image(paths[i])
    imageArray[1:size(image, 1), 1:size(image, 2), i] = image
    end
  end

  return imageArray
end

function load_overview(session, section_num)
  """
  Load overview image
  """
end

