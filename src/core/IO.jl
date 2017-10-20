global const IO_PROC = nprocs();
if nprocs() > 2
global const WORKER_PROCS = setdiff(procs(), [1, IO_PROC]);
elseif nprocs() == 2
global const WORKER_PROCS = setdiff(procs(), [1]);
else
global const WORKER_PROCS = [1];
end
global const STORAGE_OBJECTS = ["mesh", "match", "meshset"]

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

function get_mip()
  return PARAMS[:match][:mip]
end

function get_scale()
  return 1/(2^get_mip())
end

function get_z(index::Number)
  return Int(floor(index))
end

function is_subsection(index::Number)
  return get_subsection(index) != 0
end

function get_subsection(index::Number)
  return Int(round((index % get_z(index)) * SPLIT_MESH_BASIS))
end

function get_z_range()
  return PARAMS[:match][:z_start]:PARAMS[:match][:z_stop]
end

function get_name(obj)
  return string(get_index(obj))
end

function get_preceding(index::Number, n=1)
  return get_z(index)-n
end

function get_succeeding(index::Number, n=1)
  return get_z(index)+n
end

function use_cache()
  return PARAMS[:dirs][:cache]
end

function is_preceding(indexA::Number, indexB::Number, within=1)
  for i in 1:within 
    if indexA == get_preceding(indexB, i)
      return true
    end
  end
  return false;
end

## Non-image related (Storage)

function load(fn, obj_name::AbstractString)
  println("Loading $obj_name for $fn")
  if obj_name in STORAGE_OBJECTS
    s = StorageWrapper(get_path(obj_name))
    return s[fn]
  else
    println("Not a Storage object, use `get_image`")
  end
end

function save(fn, obj)
  obj_name = lowercase(string(typeof(obj)))
  if obj_name in STORAGE_OBJECTS
    s = StorageWrapper(get_path(obj_name))
    s[fn] = obj
  else
    println("Not a Storage object")
  end
end

## Image related (CloudVolume)

function get_cloudvolume(obj_name::AbstractString)
  path = get_path(obj_name)
  mip = get_mip()
  cache = use_cache()
  return CloudVolumeWrapper(path, mip=mip, cache=cache, 
                                  bounded=false, fill_missing=true)
end

function get_image_size(obj_name::AbstractString)
  cv = get_cloudvolume(obj_name)
  return size(cv)[1:2]
end

function get_offset(obj_name::AbstractString)
  cv = get_cloudvolume(obj_name)
  return offset(cv)[1:2]
end

"""
Get 3-tuple of ranges representing index in cloudvolume
"""
function get_image_slice(index::Number, obj_name)
  cv = get_cloudvolume(obj_name)
  o = get_offset(obj_name)
  s = get_image_size(obj_name)
  z = get_z(index)
  xy_slice = map(range, (o, s)...)
  return tuple(xy_slice..., z:z)
end

function get_image(index::Number, obj_name::AbstractString)
  cv = get_cloudvolume(obj_name)
  slice = get_image_slice(index, obj_name)
  return get_image(cv, slice)
end

function get_image(index::Number, obj_name::AbstractString, slice)
  cv = get_cloudvolume(obj_name)
  return get_image(cv, slice)
end

function get_image(cv::CloudVolumeWrapper, slice)
  # return OffsetArray(cv[slice...], slice[1:2])
  img = cv[slice...]
  shared_img = SharedArray{UInt8}(size(img))
  shared_img[:] = img[:]
  return shared_img
end

"""
Snap to interval-aligned value, given start & interval values (using method)
"""
function snap(x, x_start, x_inter, method=floor)
   return method(Int, (x-x_start)/x_inter)*x_inter+x_start
end

"""
Snap slice to be chunk-aligned
"""
function chunk_align(obj_name::AbstractString, slice)
  cv = get_cloudvolume(obj_name)
  o = offset(cv)
  c = chunks(cv)
  x_start = snap(slice[1][1], o[1], c[1], floor)
  x_stop = snap(slice[1][end], o[1], c[1], ceil)-1
  y_start = snap(slice[2][1], o[2], c[2], floor)
  y_stop = snap(slice[2][end], o[2], c[2], ceil)-1
  z_start = snap(slice[3][1], o[3], c[3], floor)
  z_stop = snap(slice[3][end], o[3], c[3], ceil)
  return x_start:x_stop, y_start:y_stop, z_start:z_stop
end

"""
Rescope image to be chunk-aligned with cloudvolume data
"""
function chunk_align(obj_name::AbstractString, img, src_slice)
  dst_slice = chunk_align(obj_name, src_slice)
  return rescope(img, src_slice, dst_slice), dst_slice
end

function save_image(index::Number, obj_name::AbstractString, src_img, src_slice)
  cv = get_cloudvolume(obj_name)
  dst_slice = chunk_align(obj_name, src_slice)
  dst_img = rescope(src_img, src_slice, dst_slice)
  return save_image(cv, dst_slice, dst_img)
end

"""
Save 2D image to slice in CloudVolume
"""
function save_image(cv::CloudVolumeWrapper, slice, img)
  cv[slice...] = reshape(img, (size(img)..., 1, 1))
end
