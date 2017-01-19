# size in bytes
global const IMG_CACHE_SIZE = 20 * 2^30 # n * gibibytes
global const IMG_ELTYPE = UInt8

if myid() == 1
	global const IMG_CACHE_DICT = Dict{Tuple{AbstractString, Float64}, SharedArray}()
	global const IMG_CACHE_LIST = Array{Any, 1}();
end

global const IO_PROC = nprocs();
if nprocs() > 2
global const WORKER_PROCS = setdiff(procs(), [1, IO_PROC]);
elseif nprocs() == 2
global const WORKER_PROCS = setdiff(procs(), [1]);
else
global const WORKER_PROCS = [1];
end

#=
function save(path::AbstractString, img::Array)
	f = h5open(path, "w")
	@time f["img", "chunk", (1000,1000)] = img
        close(f)
end=#

function save(path::AbstractString, data; chunksize = 1000)
	println("Saving $(typeof(data)) to ", path)
  	ext = splitext(path)[2];
	if ext == ".h5"
		f = h5open(path, "w"); @time f["img", "chunk", (chunksize, chunksize)] = (typeof(data) <: SharedArray ? data.s : data); close(f)
	elseif ext == ".jls"
	open(path, "w") do file serialize(file, data) end
      	elseif ext == ".jld"
	open(path, "w") do file write(file, "data", data) end
      	end
end

function load(path::AbstractString)
	println("Loading data from ", path)
	if !isfile(path) return nothing end
	ext = splitext(path)[2];
	if ext == ".h5"
    data = get_image_disk(path)
	elseif ext == ".jls"
  	data = open(deserialize, path)
    # data = open(path)
	elseif ext == ".jld"
  	data = load(path, "data")
  elseif ext == ".png"
    data = load_mask(path)
  elseif ext == ".jpg"
    data = load_mask(path)
  elseif ext == ".txt"
    data = readdlm(path)
  elseif ext == ".json"
    data = JSON.parsefile(path; dicttype=Dict, use_mmap=true)
	end
  println("Loaded $(typeof(data)) from ", path)
	return data
end

function save(data; kwargs...)
	save(get_path(data), data; kwargs...)
end

function load(args...)
	return load(get_path(args...));
end

# extensions:
# Mesh.jl	get_image(mesh::Mesh)
# filesystem.jl	get_image() 
function get_image_disk(path::AbstractString, dtype = IMG_ELTYPE; shared = false)
	ext = splitext(path)[2];
  if ext == ".tif" || ext == ".png"
    img = data(FileIO.load(path))
    # img = img[:, :, 1]' # incredibly slow in Julia v0.5.0
    #img.properties["timedim"] = 0
    if shared
      img = convert(Array{dtype, 2}, round(convert(Array, img)*255))'
      shared_img = SharedArray(dtype, size(img)...);
      @inbounds shared_img.s[:,:] = img[:,:]
      return shared_img;
    else
      if dtype == UInt8
        return convert(Array{dtype, 2}, round(convert(Array, img)*255))'
      else
        return convert(Array{dtype, 2}, img)'
      end
    end
	elseif ext == ".h5"
    if shared
      img = h5read(path, "img")
      if typeof(img) != Array{dtype, 2} img = convert(Array{dtype, 2}, img) end
      @everywhere gc();
      shared_img = SharedArray(dtype, size(img)...);
      @inbounds shared_img.s[:,:] = img[:,:]
      return shared_img;
    else
      return convert(Array{dtype, 2}, h5read(path, "img"))
    end
	end
end

function prefetch(index, scale=1.0, dtype = IMG_ELTYPE)
        # clean_cache()  
	if haskey(IMG_CACHE_DICT, (get_image_path(index), scale)) index = NO_INDEX end;
	println("$(get_path(index)) is being prefetched at scale $scale...")
	return remotecall(IO_PROC, get_image_disk_async, index, scale, dtype);
end

function get_image_disk_async(index, scale=1.0, dtype = IMG_ELTYPE) 
        if index == NO_INDEX
	  return nothing;
	end
        path = get_path(index);
	img = get_image_disk(path, dtype);
	scaled_img = imscale(img, scale)[1]
	ref = RemoteRef();
	put!(ref, scaled_img);
	remotecall_fetch(1, load_image, path, scale, ref, dtype)
	return nothing;
end

function load_image(path::AbstractString, scale, img::Array, dtype = IMG_ELTYPE)
	    push!(IMG_CACHE_LIST, (path, scale))
	    IMG_CACHE_DICT[(path, scale)] = img;
end

# function load_image(path::AbstractString, scale, imgref::RemoteRef, dtype = IMG_ELTYPE)
#             img = take!(imgref);
# 	    close(imgref);
# 	    load_image(path, scale, img, dtype);
# 	    img = 0;
# 	   fetch = remotecall(IO_PROC, gc); gc(); wait(fetch)
# end

function reset_cache()
  IMG_CACHE_DICT = Dict{Tuple{AbstractString, Float64}, SharedArray}()
  IMG_CACHE_LIST = Array{Any, 1}()
  @time @everywhere gc();
  @time @everywhere gc();
end

function clean_cache()
	if sum(map(Int64, map(length, values(IMG_CACHE_DICT)))) > IMG_CACHE_SIZE && !(length(IMG_CACHE_DICT) < 2)
	while sum(map(Int64, map(length, values(IMG_CACHE_DICT)))) > IMG_CACHE_SIZE * 0.50 && !(length(IMG_CACHE_DICT) < 2)
		todelete = shift!(IMG_CACHE_LIST);
		#IMG_CACHE_DICT[todelete] = zeros(IMG_ELTYPE,0,0)
		delete!(IMG_CACHE_DICT, todelete)
	end
	print("cache garbage collection:")
	@time @everywhere gc();
      end

	cur_cache_size = sum(map(Int64, map(length, values(IMG_CACHE_DICT))));

	println("current cache usage: $cur_cache_size / $IMG_CACHE_SIZE (bytes), $(round(Int64, cur_cache_size/IMG_CACHE_SIZE * 100))%")
end

function get_image(path::AbstractString, scale=1.0, dtype = IMG_ELTYPE)
#=  	if myid() != IO_PROC
	  return remotecall_fetch(IO_PROC, get_image, path, scale, dtype);
	end =#

  	if haskey(IMG_CACHE_DICT, (path, scale))
	  println("$path is in cache at scale $scale - loading from cache...")
	  # @everywhere gc();
	  return IMG_CACHE_DICT[(path, scale)]
	end

	clean_cache();

	if !haskey(IMG_CACHE_DICT, (path, 1.0))
	    println("$path is not in cache at full scale - loading into cache...")
	    #img = get_image_disk(path, dtype);
	#    shared_img = SharedArray(dtype, size(img)...);
	#    shared_img[:,:] = img[:,:];

            #print("Getting image from disk:") 	

	    push!(IMG_CACHE_LIST, (path, 1.0))
	    #IMG_CACHE_DICT[(path, 1.0)] = img;
	    @time IMG_CACHE_DICT[(path, 1.0)] = get_image_disk(path, dtype; shared = true)
            #print("image share and store to cache:")
	    #@time IMG_CACHE_DICT[(path, 1.0)] = img
	    #img = 0;
	    #img = 0;
	#    @everywhere gc();
	end

	if scale != 1.0
	  println("$path is in cache at full scale - downsampling to scale $scale...")
	  scaled_img = imscale(sdata(IMG_CACHE_DICT[(path, 1.0)]), scale)[1]
    	  #shared_img_scaled = SharedArray(dtype, size(scaled_img)...);
	  #shared_img_scaled[:,:] = scaled_img[:,:];

	  push!(IMG_CACHE_LIST, (path, scale))
	  #@time @everywhere gc();
	  #IMG_CACHE_DICT[(path, scale)] = scaled_img;
	  push!(IMG_CACHE_DICT, (path, scale), scaled_img)
	  scaled_img = 0;
	  scaled_img = 0;
        end

	    #@everywhere gc();

	return IMG_CACHE_DICT[(path, scale)];
end

function get_image(index, scale = 1.0, dtype = IMG_ELTYPE)
  return get_image(get_path(index), scale, dtype)
end

function ufixed8_to_uint8(img)
  reinterpret(UInt8, -img)
end

function get_slice(index::Index, bb::ImageRegistration.BoundingBox, scale=1.0; is_global=true, thumb=false)
  if is_global
    if thumb
      offset = [0,0]
    else
      offset = get_offset(index)
    end
    bb = translate_bb(bb, -offset+[1,1])
  end
  return get_slice(index, bb_to_slice(bb), scale, thumb=thumb)
end

function get_slice(index::Index, slice::Tuple{UnitRange{Int64},UnitRange{Int64}}, scale=1.0; thumb=false)
  path = thumb ? get_path("thumbnail", index) : get_path(index)
  return get_slice(path, slice, scale)
end

function get_slice(path::AbstractString, slice, scale=1.0)
  dtype = UInt8
  output_bb = slice_to_bb(slice)
  scaled_output_bb = snap_bb(scale_bb(output_bb, scale))
  output = zeros(dtype, ImageRegistration.get_size(scaled_output_bb)...)
  o = [1,1]

  fid = h5open(path, "r")
  dset = fid["img"]
  data_size = size(dset)
  image_bb = ImageRegistration.BoundingBox(o..., data_size...)

  if intersects(output_bb, image_bb)
    shared_bb = image_bb - output_bb
    image_slice = bb_to_slice(shared_bb)
    img = h5read(path, "img", image_slice)

    if scale != 1.0
      img, _ = imscale(img, scale)
    end

    output_offset = ImageRegistration.get_offset(output_bb)
    output_roi = translate_bb(shared_bb, -output_offset+o)
    output_roi = translate_bb(scale_bb(translate_bb(output_roi, -o), scale), o)
    output_slice = bb_to_slice(snap_bb(output_roi))
    output[output_slice...] = img
  end
  return output
end

function load_mask(index::Index; clean=true)
    path = get_path("mask", index)
    return load_mask(path, clean=clean)
end

function load_mask(path; clean=true)
    img = Images.load(path).data'
    mask = reinterpret(UInt8, img)
    mask = permutedims(mask[1,:,:], [2, 3, 1])[:,:,1]
    if clean
        clean_mask!(mask)
    end
    return mask
end

function load_affine(path::AbstractString)
  affinePath = joinpath(AFFINE_DIR_PATH, string(path, ".csv"))
  return readcsv(path)
end

function parse_rough_align(info_path::AbstractString)
  file = readdlm(info_path)
  session = cell(size(file, 1), 4); # name, index, dx, dy
  for i in 1:size(file, 1)
  m = Base.match(r"(Tile\S+).tif", file[i, 1])
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
  paths = Array{AbstractString, 1}(num_tiles)
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

function expunge_tiles(indices)
  for index in indices
    assert(is_premontaged(index))
    src_path = get_path(index)
    dst_path = get_path("expunge", index)
    println("Expunging $index")
    println("src_path: $src_path")
    println("dst_path: $src_path")
    assert(isfile(src_path) && !isfile(dst_path))
    mv(src_path, dst_path)
  end
  # metadata = get_metadata(index)
  purge_from_registry!(indices)
end

function expunge_tile(index::Index)
  assert(is_premontaged(index))
  src_path = get_path(index)
  dst_path = get_path("expunge", index)
  println("Expunging $index")
  println("src_path: $src_path")
  println("dst_path: $src_path")
  assert(isfile(src_path) && !isfile(dst_path))
  mv(src_path, dst_path)
  # metadata = get_metadata(index)
  purge_from_registry!(index)
end

function expunge_section(index::Index)
  tiles = get_index_range(premontaged(index), premontaged(index))
  for tile in tiles
    expunge_tile(tile)
  end
  purges = []
  if is_finished(index)
    purges = [aligned(index), prealigned(index), montaged(index)]
  elseif is_aligned(index)
    purges = [aligned(index), prealigned(index), montaged(index)]
  elseif is_prealigned(index)
    purges = [prealigned(index), montaged(index)]
  elseif is_montaged(index)
    purges = [montaged(index)]
  end
  for purge in purges
    purge_from_registry!(purge)
  end
end

function resurrect_tile(index::Index)
  assert(is_premontaged(index))
  fn = get_name(index)
  path = joinpath(EXPUNGED_DIR, fn)
  new_path = get_path(index)
  println("Resurrecting $fn")
  println("Moving back to $new_path")
  assert(isfile(path) && !isfile(new_path))
  mv(path, new_path)
end

function is_expunged(index::Index)
  assert(is_premontaged(index))
  fn = get_filename(index)
  expunged_path = joinpath(EXPUNGED_DIR, fn)
  included_path = get_path(index)
  return assert(isfile(expunged_path) && !isfile(included_path))
end

function make_stack(firstindex::Index, lastindex::Index, slice=(1:255, 1:255); scale=1.0, thumb=false)
  # dtype = h5read(get_path(firstindex), "dtype")
  dtype = UInt8
  global_bb = ImageRegistration.slice_to_bb(slice)
  stack_offset = ImageRegistration.get_offset(global_bb)
  stack_size = ImageRegistration.get_size(global_bb)

  thumbnail_scale = 0.02
  if thumb
    thumb_path = get_path("thumbnail", firstindex)
    thumbnail_scale = h5read(thumb_path, "scale")
  end

  stack_bb = sz_to_bb(stack_size)
  sbb = scale_bb(stack_bb, scale)
  scaled_size = sbb.h, sbb.w
  # scaled_size = round(Int64, stack_size[1]*scale+1), round(Int64, stack_size[2]*scale+1)

  imgs = zeros(dtype, scaled_size..., length(get_index_range(firstindex, lastindex)))

  for (i, index) in enumerate(get_index_range(firstindex, lastindex))
    print(string(join(index[1:2], ",") ,"|"))
    #img = zeros(dtype, scaled_size...)
    offset = thumb ? [0,0] : get_offset(index)
    sz = get_image_size(index)
    bb = snap_bb(ImageRegistration.BoundingBox(offset..., sz...))
    if intersects(bb, global_bb)
      shared_bb = global_bb - bb
      stack_roi = snap_bb(scale_bb(translate_bb(shared_bb, -stack_offset+[1,1]), scale))
      #stack_roi = snap_bb(translate_bb(scale_bb(translate_bb(shared_bb, -stack_offset), scale), [1,1]))
      img_slice = bb_to_slice(stack_roi)
      imgs[img_slice..., i] = get_slice(index, shared_bb, scale, is_global=true, thumb=thumb) 
    end
#    push!(imgs, img)
  end
#  return cat(3, imgs...)
  return imgs
end

function make_slice(center, radius)
  x, y = center
  return (x-radius+1):(x+radius), (y-radius+1):(y+radius)
end

function save_stack(firstindex::Index, lastindex::Index, center, radius; scale=1.0)
  slice = make_slice(center, radius)
  stack = make_stack(firstindex, lastindex, slice, scale=scale)
  return save_stack(stack, firstindex, lastindex, slice, scale=scale)
end

function save_stack(firstindex::Index, lastindex::Index, slice=(1:200, 1:200); scale=1.0)
  stack = make_stack(firstindex, lastindex, slice, scale=scale)
  return save_stack(stack, firstindex, lastindex, slice, scale=scale)
end

function save_stack(stack::Array{UInt8,3}, firstindex::Index, lastindex::Index, slice=(1:200, 1:200); scale=1.0)
  # if sum(stack) == 0 return end
  # perm = [3,2,1]
  # stack = permutedims(stack, perm)
  orientation = "zyx"
  dataset = DATASET
  origin = [0,0]
  x_slice = [slice[1][1], slice[1][end]] + origin
  y_slice = [slice[2][1], slice[2][end]] + origin
  z_slice = [find_in_registry(firstindex), find_in_registry(lastindex)]
  phasename = is_prealigned(firstindex) ? "prealigned" : "aligned"
  filename = string(phasename, "_", join([join(x_slice, "-"), join(y_slice, "-"), join(z_slice,"-")], "_"), ".h5")
  filepath = joinpath(FINISHED_DIR_PATH, filename)
  println("\nSaving stack to ", filepath)
  f = h5open(filepath, "w")
  # Omni can't handle chunked channel data
  # chunksize = min(512, min(size(stack)...))
  # @time f["main", "chunk", (chunksize,chunksize,chunksize)] = stack
  f["img"] = stack
  f["orientation"] = orientation
  f["origin"] = origin
  f["x_slice"] = x_slice
  f["y_slice"] = y_slice
  f["z_slice"] = z_slice
  f["z_start_index"] = [firstindex...]
  f["z_end_index"] = [lastindex...] 
  f["scale"] = scale
  f["by"] = ENV["USER"]
  f["machine"] = gethostname()
  f["timestamp"] = string(now())
  f["dataset"] = DATASET
  close(f)
end

function build_range(origin, cube_dim, overlap, grid_dim, i)
  return origin[i] : cube_dim[i]-overlap[i] : origin[i] + (cube_dim[i]-overlap[i]) * grid_dim[i]-1
end

"""
Save out a grid of cubes starting from the orgin

All dimensions in i,j,k format
"""
function save_cubes(cube_origin::Tuple{Int64, Int64, Int64}, cube_dims::Tuple{Int64, Int64, Int64}, overlap::Tuple{Int64, Int64, Int64}, grid_dims::Tuple{Int64, Int64, Int64})
  indices = get_indices(aligned(0,0))
  start_i = build_range(cube_origin, cube_dims, overlap, grid_dims, 1)
  start_j = build_range(cube_origin, cube_dims, overlap, grid_dims, 2)
  start_k = build_range(cube_origin, cube_dims, overlap, grid_dims, 3)
  x,y,z = cube_dims
  n = 1
  for i in start_i
    for j in start_j
      for k in start_k
        firstindex = indices[k]
        lastindex = indices[k+z-1]
        slice = (i:i+x-1, j:j+y-1)
        println("processing cube $n of $(length(start_i) * length(start_j) * length(start_k)) at ", join([firstindex, lastindex, slice], ", "))
        n += 1
        save_stack(firstindex, lastindex, slice)
      end
    end
  end
end

"""
Return bool indicating if image has mask
"""
function image_has_mask(path::AbstractString)
  if ishdf5(path)
    f = h5open(path, "r")
    return HDF5.exists(attrs(f), "masks")
  else
    return false
  end
end

"""
Given image size, create Array{Bool, 2} for all of image

Key
  false:  padding
  true:   tissue

Note: HDF5.jl does not currently support H5T_BITFIELD, so the BitArray is 
reinterpretted into an Array{UInt8,2}. In order to reinterpret it back
into a proper BitArray, the size of the image needs to be given.

Might consider implementing the mask as a set of polygons, using fillpoly!
in ImageRegistration.jl to generate the image mask when needed.
"""
function create_mask(sz::Array{Int,2})
  return ones(Bool, sz)
end

"""
Given HDF5 path, return Array{Bool,2} representing image padding mask
"""
function get_mask(path::AbstractString)
  f = h5open(path, "r")
  sz = size(f["img"])
  if "mask" in names(f)
    return convert(Array{Bool,2}, UInt8_to_BitArray(f["mask"], sz))
  else
    return create_mask(sz)
  end
end

function BitArray_to_UInt8(B)
  return reinterpret(UInt8, B.chunks)
end

function UInt8_to_BitArray(b, sz)
  B = BitArray(sz)
  B.chunks = reinterpret(UInt64, b)
  return B
end

function compile_save_stacks(indices)
  stack = []
  for index in indices
    path = get_path(finished(index))
    if !isfile(path)
      println(index)
    end
  end
  for index in indices
    println(index)
    path = get_path(finished(index))
    img = h5read(path, "img", (1:5000,1:20))
    push!(stack, img)
  end
  s = cat(3, stack...)
  path = get_path(finished(indices[1]))
  f = h5open(path)
  x_slice = colon(f["x_slice"][:]...)
  y_slice = colon(f["y_slice"][:]...)
  close(f)
  save_stack(s, indices[1], indices[end], (x_slice, y_slice), scale=1)
end