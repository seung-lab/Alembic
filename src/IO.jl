

function get_name(index::Index)
    if is_overview(index)
        if cur_dataset == "zebrafish"
            return string("MontageOverviewImage_W00", index[1], "_sec", index[2])
        else
            return string("MontageOverviewImage_S2-W00", index[1], "_sec", index[2])
        end
    elseif is_montaged(index)
    return string(index[1], ",", index[2], "_montaged")
    elseif is_prealigned(index)
    return string(index[1], ",", index[2], "_prealigned")
    elseif is_aligned(index)
    return string(index[1], ",", index[2], "_aligned")
    else
    return string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
    end
end

# function get_path()
# methods: 
#     
# extensions:
# Mesh.jl: get_path(mesh::Mesh)
function get_path(index, ext = ".h5")
    name = get_name(index)
    if is_overview(index)
        if cur_dataset == "zebrafish"
            section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
        else
            section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        end
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ext))
    elseif is_montaged(index)
        path = joinpath(MONTAGED_DIR, string(name, ext))
    elseif is_prealigned(index)
        path = joinpath(PREALIGNED_DIR, string(name, ext))
    elseif is_aligned(index)
        path = joinpath(ALIGNED_DIR, string(name, ext))
    else
        if cur_dataset == "zebrafish"
            section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
        else
            section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        end
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ext))
    end
    println(path)
    return path
end

function get_path(name::String)
    return get_path(parse_name(name))
end

function get_image(index::Index)
    return get_image(get_path(index))
end

function get_thumbnail_path(index::Index)
  fn = string(join(index[1:2], ","), "_thumb.jpg")
  println(fn)
  return joinpath(MONTAGED_DIR, "review", fn) 
end

function get_thumbnail_path(indexA::Index, indexB::Index)
  fn = string(join(indexA[1:2],","), "-", join(indexB[1:2],","), "_thumb.png")
  println(fn)
  if is_prealigned(indexA)
    return joinpath(PREALIGNED_DIR, "review", fn) 
  else
    return joinpath(ALIGNED_DIR, "review", fn) 
  end
end

# extensions:
# Mesh.jl	get_image(mesh::Mesh)
# filesystem.jl	get_image() 
function get_image(path::String, dtype = UInt8)
	ext = splitext(path)[2];
  	if ext == ".tif"
		img = imread(path)
		img = img[:, :, 1]
 		img.properties["timedim"] = 0
  		return convert(Array{dtype, 2}, round(convert(Array, img)*255))
	elseif ext == ".h5"
 		return convert(Array{dtype, 2}, h5read(path, "img"))
	end
end

function get_uint8_image(path::String)
  img = imread(path)
  return reinterpret(UInt8, data(img)[:,:,1])
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
  return convert(Array{Ufixed8}, h5read(path, "img", slice))
end

function get_h5_image(path::String)
  return convert(Array{Ufixed8}, h5read(path, "img"))
end

# extensions:
# Mesh.jl get_float_image(mesh::Mesh)
function get_float_image(path::String)
  img = imread(path)
  img.properties["timedim"] = 0
  return convert(Array{Float64, 2}, convert(Array, img[:,:,1]))
end

function get_ufixed8_image(path::String)
  return convert(Array{Ufixed8}, data(imread(path))[:,:,1])'
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

