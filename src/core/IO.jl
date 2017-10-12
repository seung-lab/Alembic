#global const IMG_ELTYPE = UInt8
global IMG_ELTYPE = UInt8

global const IO_PROC = nprocs();
if nprocs() > 2
global const WORKER_PROCS = setdiff(procs(), [1, IO_PROC]);
elseif nprocs() == 2
global const WORKER_PROCS = setdiff(procs(), [1]);
else
global const WORKER_PROCS = [1];
end


# generates Janelia-type tilespec from the registry
function generate_tilespec_from_registry(index; dataset = "", stack = "", stack_type = "", parent_stack = "")
  tileid = "$(index[1]),$(index[2])_$stack"
  tilespec = Dict{Symbol, Any}()
  tilespec[:tileId] = tileid
  tilespec[:z] = find_in_registry(index)
  tilespec[:width] = get_image_size(index)[2]
  tilespec[:height] = get_image_size(index)[1]
  tilespec[:minX] = get_offset(index)[2]
  tilespec[:minY] = get_offset(index)[1]
  # we need the meta field w
  tilespec[:meta] = Dict{Symbol, Any}()
  tilespec[:meta][:dataset] = dataset
  tilespec[:meta][:stack] = stack
  tilespec[:meta][:stack_type] = stack_type
  tilespec[:meta][:parent_stack] = parent_stack
  tilespec[:meta][:resolution] = (30,30,40)
  tilespec[:channels] = Dict{Symbol, Any}()
  tilespec[:channels][:nccnet] = Dict{Symbol, Any}()
  tilespec[:channels][:nccnet][:imageUrl] = "file://$(get_path(bucket = BUCKET, dataset = dataset, stack = stack, obj_name = "nccnet", tileid = tileid))"
  #tilespec[:render] = Dict{Symbol, Any}()
  #tilespec[:render][:transform] = "file://$(get_path(bucket = BUCKET, dataset = dataset, stack = stack, obj_name = "cumulative_tform", tileid = tileid))"
  tilespec[:mipmapLevels] = Dict{Symbol, Any}()
  tilespec[:mipmapLevels][Symbol(0)] = Dict{Symbol, Any}()
  tilespec[:mipmapLevels][Symbol(0)][:imageUrl] = "file://$(get_path(bucket = BUCKET, dataset = dataset, stack = stack, obj_name = "image", tileid = tileid))"
  #tilespec[:mipmapLevels][:0][:imageUrl] = "file://$BUCKET/$(get_path(dataset = dataset, stack = stack, tileid = tileid))"
  return tilespec
end

function generate_stackspec(dataset = "", stack = "")
  stackspec = Dict{Symbol, Any}()
  tileid = "$(index[1]),$(index[2])_$stack"
  tilespec = Dict{Symbol, Any}()
  tilespec[:tileId] = tileid
  tilespec[:z] = find_in_registry(index)
  tilespec[:width] = get_image_size(index)[2]
  tilespec[:height] = get_image_size(index)[1]
  tilespec[:minX] = get_offset(index)[2]
  tilespec[:minY] = get_offset(index)[1]
  # we need the meta field w
  tilespec[:meta] = Dict{Symbol, Any}()
  tilespec[:meta][:dataset] = dataset
  tilespec[:meta][:stack] = stack
  tilespec[:meta][:parent_stack] = parent_stack
  tilespec[:meta][:resolution] = (30,30,40)
  tilespec[:mipmapLevels] = Dict{Symbol, Any}()
  tilespec[:mipmapLevels][Symbol(0)] = Dict{Symbol, Any}()
  tilespec[:mipmapLevels][Symbol(0)][:imageUrl] = "file://$(get_path(bucket = BUCKET, dataset = dataset, stack = stack, tileid = tileid))"
  #tilespec[:mipmapLevels][:0][:imageUrl] = "file://$BUCKET/$(get_path(dataset = dataset, stack = stack, tileid = tileid))"
  return tilespec
end

function get_path(obj_name::AbstractString)
  return PARAMS[:dirs][Symbol(obj_name)]
end

function get_scale()
  return PARAMS[:match][:blockmatch_scale]
end

function get_z(index):
  return index
end

function save(path::AbstractString, data; chunksize = 1000)
  println("Saving $(typeof(data)) to ", path)
    ext = splitext(path)[2];
  if ext == ".h5"
    f = h5open(path, "w"); @time f["img", "chunk", (chunksize, chunksize)] = (typeof(data) <: SharedArray ? data.s : data); close(f)
  elseif ext == ".jls"
  open(path, "w") do file serialize(file, data) end
        elseif ext == ".jld"
  open(path, "w") do file write(file, "data", data) end
        elseif ext == ".json"
  open(path, "w") do file JSON.print(file, data) end
        end
end

function load(path::AbstractString; kwargs...)
  println("Loading data from ", path)
  if !isfile(path) return nothing end
  ext = splitext(path)[2];
  if ext == ".h5"
    data = get_image_disk(path; kwargs...)
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
    #data = JSON.parsefile(path; dicttype=Dict{Symbol, Any}, use_mmap=true)
    data = JSON.parsefile(path; dicttype=Dict{Symbol, Any}, use_mmap=false)
  end
  println("Loaded $(typeof(data)) from ", path)
#=  if typeof(data) == Match
  data.correspondence_properties = Array{Dict{Any, Any},1}(); end
  gc(); =#
  return data
end

function save(data; kwargs...)
  save(get_path(data), data; kwargs...)
end

function load(args...; kwargs...)
  return load(get_path(args...); kwargs...);
end

function get_offset(obj_name, mip=get_scale())
  return offset(CloudVolumeWrapper(get_path(obj_name), mip=mip))
end

function get_image(index, mip=get_scale())
  cv = CloudVolumeWrapper(get_path("src_image"), mip=mip)
  offset = offset(cv)[1:2]
  sz = size(cv)[1:2]
  z = get_z(index)
  return cv[map(range, zip(offset, sz)...)..., z]
end

function get_slice(index::FourTupleIndex, bb::ImageRegistration.BoundingBox, scale=1.0; is_global=true, thumb=false)
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

function get_slice(index::FourTupleIndex, slice::Tuple{UnitRange{Int64},UnitRange{Int64}}, scale=1.0; thumb=false)
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

function load_mask(index::FourTupleIndex; clean=true)
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
  imageArray = SharedArray{UInt8}(max_tile_size, max_tile_size, num_tiles)

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

function make_slice(center, radius)
  x, y = center
  return (x-radius+1):(x+radius), (y-radius+1):(y+radius)
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