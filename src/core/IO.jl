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

function get_mip()
  return PARAMS[:match][:mip]
end

function get_scale()
  return 1/(2^get_mip())
end

function get_z(index)
  return Int(floor(index))
end

function is_subsection(index)
  return get_subsection(index) != 0
end

function get_subsection(index)
  return Int(round((index % get_z(index)) * SPLIT_MESH_BASIS))
end

function get_z_range()
  return PARAMS[:match][:z_start]:PARAMS[:match][:z_stop]
end

function get_name(obj)
  return string(get_index(obj))
end

function get_preceding(index, n=1)
  return get_z(index)-n
end

function use_cache()
  return PARAMS[:dirs][:cache]
end

function is_preceding(indexA, indexB, within=1)
  for i in 1:within 
    if indexA == get_preceding(indexB, i)
      return true
    end
  end
  return false;
end

## Non-image related (Storage)

function load(fn, obj_name)
  println("Loading $obj_name for $fn")
  storage_objects = ["mesh", "match"]
  if obj_name in storage_objects
    s = StorageWrapper(get_path(obj_name))
    return s[fn]
  else
    println("Not a Storage object, use `get_image`")
  end
end

function save(fn, obj)
  obj_name = lowercase(string(typeof(obj)))
  if obj_name in ["mesh", "match"]
    s = StorageWrapper(get_path(obj_name))
    s[fn] = obj
  else
    println("Not a Storage object")
  end
end

## Image related (CloudVolume)

function get_cloudvolume(obj_name)
  path = get_path(obj_name)
  mip = get_mip()
  cache = use_cache()
  return CloudVolumeWrapper(path, mip=mip, cache=cache)
end

function get_image_size(obj_name)
  cv = get_cloudvolume(obj_name)
  return size(cv)[1:2]
end

function get_offset(obj_name)
  cv = get_cloudvolume(obj_name)
  return offset(cv)[1:2]
end

function get_image(index, obj_name::AbstractString)
  cv = get_cloudvolume(obj_name)
  o = get_offset(obj_name)
  s = get_image_size(obj_name)
  z = get_z(index)
  xy_slice = map(range, (o, s)...)
  slice = tuple(xy_slice..., z)
  return get_image(cv, slice)
end

function get_image(index, obj_name::AbstractString, slice)
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

function save_image(index, obj_name::AbstractString, img)
  cv = get_cloudvolume(obj_name)
  o = get_offset(obj_name)
  s = get_image_size(obj_name)
  z = get_z(index)
  xy_slice = map(range, (o, s)...)
  slice = tuple(xy_slice..., z)
  return get_image(cv, slice, img)
end

function save_image(cv::CloudVolumeWrapper, slice, img)
  cv[slice...] = img
end
