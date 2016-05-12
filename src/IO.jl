global const IMG_ELTYPE = UInt8
#global const IMG_SUP_SIZE = (75000, 75000)

# size in bytes
global const IMG_CACHE_SIZE = 8 * 2^30 # n * gibibytes

if myid() == 1
	global const IMG_CACHE_DICT = Dict{Any, SharedArray}()
	global const IMG_CACHE_LIST = Array{Any, 1}();
end

global const IO_PROC = nprocs();
#global const WORKER_PROCS = setdiff(procs(), IO_PROC);
if nprocs() > 2
global const WORKER_PROCS = setdiff(procs(), [1, IO_PROC]);
else 
global const WORKER_PROCS = setdiff(procs(), [1]);
end

#global const SHARED_SRC_IMAGE = SharedArray(IMG_ELTYPE, SUP_SIZE)
#global const SHARED_DST_IMAGE = SharedArray(IMG_ELTYPE, SUP_SIZE)


function save(path::String, img::Array)
      f = h5open(path, "w")
      @time f["img", "chunk", (1000,1000)] = img
      close(f)
    end

# extensions:
# Mesh.jl	get_image(mesh::Mesh)
# filesystem.jl	get_image() 
function get_image_disk(path::String, dtype = IMG_ELTYPE)
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

function prefetch(index, scale=1.0, dtype = IMG_ELTYPE)
        clean_cache()  
	if haskey(IMG_CACHE_DICT, (get_path(index), scale)) index = NO_INDEX end;
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

function load_image(path::String, scale, img::Array, dtype = IMG_ELTYPE)
	    push!(IMG_CACHE_LIST, (path, scale))
	    IMG_CACHE_DICT[(path, scale)] = img;
end

function load_image(path::String, scale, imgref::RemoteRef, dtype = IMG_ELTYPE)
            img = take!(imgref);
	    close(imgref);
	    load_image(path, scale, img, dtype);
	    img = 0;
	   fetch = remotecall(IO_PROC, gc); gc(); wait(fetch)
end

function clean_cache()
	if sum(map(Int64, map(length, values(IMG_CACHE_DICT)))) > IMG_CACHE_SIZE && !(length(IMG_CACHE_DICT) < 2)
	while sum(map(Int64, map(length, values(IMG_CACHE_DICT)))) > IMG_CACHE_SIZE * 0.75 && !(length(IMG_CACHE_DICT) < 2)
		todelete = shift!(IMG_CACHE_LIST);
		IMG_CACHE_DICT[todelete] = SharedArray(IMG_ELTYPE, 0, 0);
		delete!(IMG_CACHE_DICT, todelete)
	end
	@everywhere gc();
      end

	cur_cache_size = sum(map(Int64, map(length, values(IMG_CACHE_DICT))));

	println("current cache usage: $cur_cache_size / $IMG_CACHE_SIZE (bytes), $(round(Int64, cur_cache_size/IMG_CACHE_SIZE * 100))%")
end

function get_image(path::String, scale=1.0, dtype = IMG_ELTYPE)
#=  	if myid() != IO_PROC
	  return remotecall_fetch(IO_PROC, get_image, path, scale, dtype);
	end =#

  	if haskey(IMG_CACHE_DICT, (path, scale))
	  println("$path is in cache at scale $scale - loading from cache...")
	  return IMG_CACHE_DICT[(path, scale)]
	end

	clean_cache();

	if !haskey(IMG_CACHE_DICT, (path, 1.0))
	    println("$path is not in cache at full scale - loading into cache...")
	    img = get_image_disk(path, dtype);
	 #   shared_img = SharedArray(dtype, size(img)...);
	 #   shared_img[:,:] = img[:,:];

	    push!(IMG_CACHE_LIST, (path, 1.0))
	    IMG_CACHE_DICT[(path, 1.0)] = img;
	    img = 0;
	end

	if scale != 1.0
	  println("$path is in cache at full scale - downsampling to scale $scale...")
	  scaled_img = imscale(sdata(IMG_CACHE_DICT[(path, 1.0)]), scale)[1]
    	  #shared_img_scaled = SharedArray(dtype, size(scaled_img)...);
	  #shared_img_scaled[:,:] = scaled_img[:,:];

	  push!(IMG_CACHE_LIST, (path, scale))
	  IMG_CACHE_DICT[(path, scale)] = scaled_img;
	  scaled_img = 0;
        end

	@everywhere gc();
	return IMG_CACHE_DICT[(path, scale)];
end

function get_image(index, scale = 1.0, dtype = IMG_ELTYPE)
  return get_image(get_path(index), scale, dtype)
end

function ufixed8_to_uint8(img)
  reinterpret(UInt8, -img)
end

function get_slice(path::String, slice)
  # return reinterpret(Ufixed8, h5read(path, "img", slice))
  if splitext(path)[2] != ".h5" return get_image(path)[slice...] end
  return h5read(path, "img", slice)
end

function load_affine(path::String)
  affinePath = joinpath(AFFINE_DIR, string(path, ".csv"))
  return readcsv(path)
end

function parse_rough_align(info_path::String)
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

function expunge_tile(index::Index)
  assert(is_premontaged(index))
  fn = get_filename(index)
  path = get_path(index)
  new_path = joinpath(EXPUNGED_DIR, fn)
  println("Expunging $fn")
  println("Moving to $new_path")
  assert(isfile(path) && !isfile(new_path))
  mv(path, new_path)
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
  fn = get_filename(index)
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

function make_stack_from_finished(firstindex::Index, lastindex::Index, slice=(1:200, 1:200))
    imgs = []
    bb = nothing
    for index in get_index_range(firstindex, lastindex)
        index = finished(index)
        print(string(join(index[1:2], ",") ,"|"))
        img = get_slice(get_path(index), slice)
        if bb == nothing
            bb = h5read(get_path(index), "bb")
        else
            current_bb = h5read(get_path(index), "bb")
            if current_bb != bb
                error("FINISHED IMAGE, $index, NOT IN SAME BOUNDING BOX: $bb")
            end
        end
        push!(imgs, img)
    end
    return cat(3, imgs...), bb
end

function make_stack(firstindex::Index, lastindex::Index, slice=(1:255, 1:255))
  # dtype = h5read(get_path(firstindex), "dtype")
  dtype = UInt8
  stack_offset = [slice[1][1], slice[2][1]] - [1,1]
  stack_size = map(length, slice)
  global_bb = BoundingBox(stack_offset..., stack_size...)
  imgs = []
  for index in get_index_range(firstindex, lastindex)
    print(string(join(index[1:2], ",") ,"|"))
    img = zeros(dtype, stack_size...)
    offset = get_offset(index)
    sz = get_image_size(index)
    bb = BoundingBox(offset..., sz...)
    if intersects(bb, global_bb)
      shared_bb = global_bb - bb
      img_roi = translate_bb(shared_bb, -offset)
      stack_roi = translate_bb(shared_bb, -stack_offset)
      h5_slice = bb_to_slice(img_roi)
      img_slice = bb_to_slice(stack_roi)
      # println(h5_slice, " " , img_slice)
      img[img_slice...] = get_slice(get_path(index), h5_slice) 
    end
    push!(imgs, img)
  end
  return cat(3, imgs...)
end

function save_stack(firstindex::Index, lastindex::Index, slice=(1:200, 1:200))
  stack = make_stack(firstindex, lastindex, slice)
  return save_stack(stack, firstindex, lastindex, slice)
end

function save_stack(stack::Array{UInt8,3}, firstindex::Index, lastindex::Index, slice=(1:200, 1:200))
  scale = 1.0
  # perm = [3,2,1]
  # stack = permutedims(stack, perm)
  orientation = "zyx"
  dataset = cur_dataset
  origin = [0,0]
  x_slice = [slice[1][1], slice[1][end]] + origin
  y_slice = [slice[2][1], slice[2][end]] + origin
  z_slice = [find_in_registry(aligned(firstindex)), find_in_registry(aligned(lastindex))]
  filename = string(cur_dataset, "_", join([join(x_slice, "-"), join(y_slice, "-"), join(z_slice,"-")], "_"), ".h5")
  filepath = joinpath(FINISHED_DIR, filename)
  println("\nSaving stack to ", filepath)
  f = h5open(filepath, "w")
  chunksize = min(512, min(size(stack)...))
  @time f["main", "chunk", (chunksize,chunksize,chunksize)] = stack
  # f = Dict()
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
  f["dataset"] = cur_dataset
  close(f)
end

function build_range(origin, cube_dim, overlap, grid_dim, i)
  return origin[i] : cube_dim[i]-overlap[i] : origin[i] + (cube_dim[i]-overlap[i]) * grid_dim[i]-1
end

"""
Save out a grid of cubes starting from the orgin

All dimensions in i,j,k format
"""
function save_cubes(origin::Tuple{Int64, Int64, Int64}, cube_dims::Tuple{Int64, Int64, Int64}, overlap::Tuple{Int64, Int64, Int64}, grid_dims::Tuple{Int64, Int64, Int64})
  indices = get_indices(aligned(1,2))
  start_i = build_range(origin, cube_dims, overlap, grid_dims, 1)
  start_j = build_range(origin, cube_dims, overlap, grid_dims, 2)
  start_k = build_range(origin, cube_dims, overlap, grid_dims, 3)
  x,y,z = cube_dims
  n = 1
  for i in start_i
    for j in start_j
      for k in start_k
        firstindex = indices[k]
        lastindex = indices[k+z-1]
        slice = (i:i+x-1, j:j+y-1)
        println(join([n, firstindex, lastindex, slice], ", "))
        n += 1
        save_stack(firstindex, lastindex, slice)
      end
    end
  end
end
