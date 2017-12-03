global const IO_PROC = nprocs();
if nprocs() > 2
global const WORKER_PROCS = setdiff(procs(), [1, IO_PROC]);
elseif nprocs() == 2
global const WORKER_PROCS = setdiff(procs(), [1]);
else
global const WORKER_PROCS = [1];
end
global const IMG_CACHE_SIZE = 20 * 2^30
global const IMG_ELTYPE = UInt8

if myid() == 1
  # Image Cache indexed by index, obj_name, miplevel
  global const IMG_CACHE_DICT = Dict{Tuple{Number, AbstractString, Int64}, SharedArray{IMG_ELTYPE, 2}}()
  global const IMG_CACHE_LIST = Array{Any, 1}()
end


  
function reset_cache()
  IMG_CACHE_DICT = Dict{Tuple{AbstractString, Float64}, SharedArray{IMG_ELTYPE, 2}}()
  IMG_CACHE_LIST = Array{Any, 1}()
  @time @everywhere gc();
end

function clean_cache()
	if sum(map(Int64, map(length, values(IMG_CACHE_DICT)))) > IMG_CACHE_SIZE && !(length(IMG_CACHE_DICT) < 2)
	while sum(map(Int64, map(length, values(IMG_CACHE_DICT)))) > IMG_CACHE_SIZE * 0.50 && !(length(IMG_CACHE_DICT) < 2)
		todelete = shift!(IMG_CACHE_LIST);
		IMG_CACHE_DICT[todelete] = zeros(IMG_ELTYPE,0,0)
		delete!(IMG_CACHE_DICT, todelete)
	end
	print("cache garbage collection:")
	@time gc();
      end

	cur_cache_size = sum(map(Int64, map(length, values(IMG_CACHE_DICT))));
	println("current cache usage: $cur_cache_size / $IMG_CACHE_SIZE (bytes), $(round(Int64, cur_cache_size/IMG_CACHE_SIZE * 100))%")
end

function get_path(obj_name::AbstractString)
  return joinpath(PARAMS[:dirs][:bucket], PARAMS[:dirs][:dataset], 
                                            PARAMS[:dirs][Symbol(obj_name)])
end

"""
Get mip level specified in params
"""
function get_mip(k = :match)
  return PARAMS[k][:mip]
end

function get_scale(k = :match)
  return 1/(2^get_mip(k))
end

function get_z(index::Number)
  return Int(floor(abs(index)))*sign(index)
end

function is_subsection(index::Number)
  return get_subsection(index) != 0
end

function get_subsection(index::Number)
  return Int(round((index - get_z(index)) * SPLIT_MESH_BASIS))
end

function get_z_range()
  return PARAMS[:mesh][:z_start]:PARAMS[:mesh][:z_stop]
end

function get_name()
  z_range = get_z_range()
  return string("($(z_range[1]), $(z_range[end]))")
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

function use_defect_mask()
  return PARAMS[:dirs][:mask] != ""
end

function use_roi_mask()
  return PARAMS[:dirs][:roi] != ""
end

function globalize!{T}(pts::Points{T}, offset::Array{Int,1})
  @inbounds o1 = offset[1];
  @inbounds o2 = offset[2];

  @simd for i in 1:size(pts, 2)
    @fastmath @inbounds pts[1,i] += o1;
    @fastmath @inbounds pts[2,i] += o2;
  end 
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

function load(obj_name::AbstractString, fn::AbstractString)
  println("Loading $obj_name for $fn")
  s = StorageWrapper(get_path(obj_name))
  return s[fn]
end

function save(obj, obj_name::AbstractString, fn::AbstractString=get_name(obj))
  println("Saving $fn to $(get_path(obj_name))")
  s = StorageWrapper(get_path(obj_name))
  s[fn] = obj;
end

## Image related (CloudVolume)

function get_cloudvolume(obj_name::AbstractString; mip::Int64 = get_mip(:match))
  path = get_path(obj_name)
  cache = use_cache()
  return CloudVolumeWrapper(path, mip=mip, cache=cache, 
                                  bounded=false, fill_missing=true)
end

function get_image_size(obj_name::AbstractString; mip::Int64 = get_mip(:match))
  cv = get_cloudvolume(obj_name, mip=mip)
  return size(cv)[1:2]
end

function get_offset(obj_name::AbstractString; mip::Int64 = get_mip(:match))
  cv = get_cloudvolume(obj_name, mip=mip)
  return offset(cv)[1:2]
end

"""
Get 3-tuple of ranges representing index in cloudvolume
"""
function get_image_slice(index::Number, obj_name; mip::Int64 = get_mip(:match))
  cv = get_cloudvolume(obj_name, mip=mip)
  o = get_offset(obj_name, mip=mip)
  s = get_image_size(obj_name, mip=mip)
  z = get_z(index)
  xy_slice = map(range, (o, s)...)
  return tuple(xy_slice..., z)
end

function get_image(index::Number, obj_name::AbstractString; mip::Int64 = get_mip(:match))
  if haskey(IMG_CACHE_DICT, (index, obj_name, mip))
    println("$index, $obj_name, at miplevel $mip is in the image cache.")
  else
    println("$index, $obj_name, at miplevel $mip is not in the image cache. Downloading...")
    @time begin
    push!(IMG_CACHE_LIST, (index, obj_name, mip))

    cv = get_cloudvolume(obj_name, mip=mip)
    slice = get_image_slice(index, obj_name, mip=mip)

    IMG_CACHE_DICT[(index, obj_name, mip)] = SharedArray(get_image(cv, slice))
    end
  end
    return IMG_CACHE_DICT[(index, obj_name, mip)]::SharedArray{IMG_ELTYPE, 2}
end

function get_image(index::Number, obj_name::AbstractString, slice; mip::Int64 = get_mip(:match))
  cv = get_cloudvolume(obj_name, mip=mip)
  return get_image(cv, slice)
end

function get_image(cv::CloudVolumeWrapper, slice)
  # return OffsetArray(cv[slice...], slice[1:2])
  img = cv[slice...]
  return img
  #=
  shared_img = SharedArray{UInt8}(size(img))
  shared_img[:] = img[:]
  return shared_img
  =#
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
function chunk_align(obj_name::AbstractString, slice; mip::Int64=get_mip(:render))
  cv = get_cloudvolume(obj_name, mip=get_mip(:render))
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
function chunk_align(obj_name::AbstractString, img, src_slice; mip::Int64=get_mip(:render))
  dst_slice = chunk_align(obj_name, src_slice, mip=get_mip(:render))
  return rescope(img, src_slice, dst_slice), dst_slice
end

function save_image(index::Number, obj_name::AbstractString, src_img, src_slice; mip::Int64=get_mip(:render))
  println("Saving $index @ mip=$mip to $(get_path(obj_name))")
  cv = get_cloudvolume(obj_name, mip=mip)
  dst_slice = chunk_align(obj_name, src_slice)
  if dst_slice != src_slice
    src_img = rescope(src_img, src_slice, dst_slice)
  end
  return save_image(cv, dst_slice, src_img)
end

"""
Save 2D image to slice in CloudVolume
"""
function save_image(cv::CloudVolumeWrapper, slice, img)
  cv[slice...] = reshape(img, (size(img)..., 1, 1))
end

function get_local_root()
  return homedir()
end

function get_local_dir(dir="match")
  return joinpath(".alembic", get_path(dir)[6:end]) # hardcoded for gs:// header
end

function get_local_path(dir="match")
  return joinpath(get_local_root(), get_local_dir(dir))
end

function check_local_dir(dir="match")
  root = get_local_root()
  path = get_local_dir(dir)
  folders = split(path, '/')
  for i in 1:length(folders)
    temp_path = joinpath(root, folders[1:i]...)
    if !isdir(temp_path)
      mkdir(temp_path)
    end
  end
end

function flush_cv_cache(obj_name; mip::Int64=get_mip(:render))
  cv = get_cloudvolume(obj_name, mip=mip)
  flush(cv)
end

# function make_roi_mask(dir,)
