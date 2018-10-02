# Immutable type that contains the preallocated space needed for match
immutable MatchEnv
  # sizes of the full size patches
  size_src::NTuple{2, Int64}
  size_dst::NTuple{2, Int64}

  # blockmatch scale
  scale::Float64
  # bandpass values
  bandpass

  # bandpass filter for xcorr
  kernel::Array{Float64, 2}

  # tforms for scaling
  tform_src::Array{Float64, 2}
  tform_dst::Array{Float64, 2}

  # locations for the full size patches
  src_patch_full::Array{Float64, 2}
  dst_patch_full::Array{Float64, 2}

  # locations for the scaled patches
  src_patch::Array{Float64, 2}
  dst_patch::Array{Float64, 2}
end

global MATCH_ENVS = Dict{Symbol, MatchEnv}()

function clear_matchenvs()
    for key in keys(MATCH_ENVS)
        MATCH_ENVS[key] = MatchEnv(0:0,0:0)
    end
    gc();
    global MATCH_ENVS = Dict{Symbol, MatchEnv}()
end

function register_matchenv(me::MatchEnv)
    MATCH_ENVS[Symbol(me.size_src, me.size_dst, me.scale, me.bandpass)] = me
end

function get_matchenv(src_range, dst_range; scale = 1.0, bandpass = (0,0))
  if haskey(MATCH_ENVS, Symbol(map(length, src_range), map(length, dst_range), scale, bandpass)) 
    return MATCH_ENVS[Symbol(map(length, src_range), map(length, dst_range), scale, bandpass)]
  end

  me = MatchEnv(src_range, dst_range; scale = scale, bandpass = bandpass)
  register_matchenv(me)

  return me
end

function MatchEnv(src_range, dst_range; scale = 1.0, bandpass = (0,0))
  size_src = map(length, src_range)
  size_dst = map(length, dst_range)

  src_patch_full = zeros(Float64, size_src...)
  dst_patch_full = zeros(Float64, size_dst...)

  if scale == 1.0
    src_patch = zeros(Float64, size_src...)
    dst_patch = zeros(Float64, size_dst...)
    tform_src = zeros(Float64, 0, 0)
    tform_dst = zeros(Float64, 0, 0)
  else
    tform = [scale 0 0; 0 scale 0; 0 0 1]

    bb_src = BoundingBox{Float64}(0,0, size_src...)
    wbb_src = tform_bb(bb_src, tform)
    tbb_src = snap_bb(wbb_src)
    src_patch = zeros(Float64, tbb_src.h, tbb_src.w)
    tform_src = [tbb_src.h/size_src[1] - eps 0 0; 0 tbb_src.w/size_src[2] - eps 0; 0 0 1]

    bb_dst = BoundingBox{Float64}(0,0, size_dst...)
    wbb_dst = tform_bb(bb_dst, tform)
    tbb_dst = snap_bb(wbb_dst)
    dst_patch = zeros(Float64, tbb_dst.h, tbb_dst.w)
    tform_dst = [tbb_dst.h/size_dst[1] - eps 0 0; 0 tbb_dst.w/size_dst[2] - eps 0; 0 0 1]
  end
   
  if bandpass == (0, 0)
    kernel = zeros(Float64, 0, 0)
  else
    kernel = make_bandpass_kernel(bandpass...) 
  end

    return MatchEnv(
  size_src::NTuple{2, Int64},
  size_dst::NTuple{2, Int64},

  scale::Float64,
  bandpass,

  kernel::Array{Float64, 2},

  tform_src::Array{Float64, 2},
  tform_dst::Array{Float64, 2},

  src_patch_full::Array{Float64, 2},
  dst_patch_full::Array{Float64, 2},

  src_patch::Array{Float64, 2},
  dst_patch::Array{Float64, 2}
    )
end

function clean!(me::MatchEnv)
  @inbounds begin
    me.src_patch_full[:] = Float64(0)
    me.dst_patch_full[:] = Float64(0)
    me.src_patch[:] = Float64(0)
    me.dst_patch[:] = Float64(0)
  end
end



function init_Match()

global SRC_PATCH_FULL = Array{Float64, 2}(1,1);
global SRC_PATCH_FULL_H = Array{Float64, 2}(1,1);
global SRC_PATCH = Array{Float64, 2}(1,1);
#global SRC_PATCH_G_H = Array{Float64, 2}(1,1);
#global SRC_PATCH_G_L = Array{Float64, 2}(1,1);
global DST_PATCH_FULL = Array{Float64, 2}(1,1);
global DST_PATCH_FULL_H = Array{Float64, 2}(1,1);
global DST_PATCH = Array{Float64, 2}(1,1);
#global DST_PATCH_G_H = Array{Float64, 2}(1,1);
#global DST_PATCH_G_L = Array{Float64, 2}(1,1);

end

init_Match();

struct Match{T} <: AbstractMatch
  src_index          # source mesh index
  dst_index          # destination mesh index
  src_points::Points{T}     # source point coordinate in image
  dst_points::Points{T}        # destination point coordinate in image
  correspondence_properties    # in the same array order as points, contains properties of correspondences

  filters                       # Array of filters 
  properties::Dict{Symbol, Any}
end


function Base.deepcopy(m::Match; src_index=m.src_index, dst_index=m.dst_index,
                        src_points=m.src_points, dst_points=m.dst_points,
                        correspondence_properties=m.correspondence_properties,
                        filters=m.filters, properties=m.properties)
  return Match(src_index, dst_index, deepcopy(src_points), deepcopy(dst_points), 
                        deepcopy(correspondence_properties), deepcopy(filters), 
                        deepcopy(properties))
end

## IO
function get_correspondences_df(match::Match; globalized=true, manual_only=true)
    pre, post, filter = get_correspondences(match, globalized=globalized, manual_only=manual_only)
    pre_z, post_z = get_index(match)
    n = length(filter)
    df = DataFrame(vcat(pre, ones(n)'*pre_z, post, ones(n)'*post_z, filter')')
    names!(df, [:pre_x, :pre_y, :pre_z, :post_x, :post_y, :post_z, :filter])
    return df
end

function to_dataframe(match::Match; manual_only=true)
    return hcat(get_correspondences_df(match, manual_only=manual_only), match.correspondence_properties)
end

function to_csv(match::Match; manual_only=true)
    check_local_dir(:match)
    fn = joinpath(get_local_path(:match), get_name(match))
    println("Writing correspondences to $fn")
    # writetable(fn, to_dataframe(match, manual_only=manual_only), separator=',')
    CSV.write(fn, to_dataframe(match, manual_only=manual_only))
end

function load_manual_filters!(match::Match)
    fn = joinpath(get_local_path(:match), get_name(match))
    println("Loading manual inspection filters from $fn")
    # df = readtable(fn)
    df = CSV.read(fn)
    match.filters[:manual] = .~convert(Array{Bool,1}, df[:filter])
end

### index related
function get_src_index(match::Match)
  return match.src_index
end
function get_dst_index(match::Match)
  return match.dst_index
end
function get_index(match::Match)
    return get_src_and_dst_indices(match::Match)
end
function get_src_and_dst_indices(match::Match)
  return match.src_index, match.dst_index
end

### counting
function count_correspondences(match::Match) return size(match.src_points, 2);  end
function count_filtered_correspondences(match::Match; manual_only=false) return sum(get_filtered_indices(match, manual_only=manual_only)); end
function count_rejected_correspondences(match::Match; manual_only=false) return sum(get_rejected_indices(match, manual_only=manual_only)); end

function count_filters(match::Match) return size(match.filters, 2); end
function count_filtered_properties(match::Match, property_name, compare, threshold)
    return sum(compare(get_correspondence_properties(match::Match, property_name; filtered = true), threshold))
end

### correspondences related
function get_filtered_indices(match::Match; manual_only=false)
    return map(!, get_rejected_indices(match; manual_only=manual_only))
end
function get_rejected_indices(match::Match; manual_only=false)
    ret = fill(false, count_correspondences(match));
    if manual_only
      if :manual in names(match.filters)
        ret = match.filters[:manual]
      end
    else
      for j in 1:count_filters(match)
        dfcol::Array{Bool, 1} = match.filters.columns[j]
        for i in 1:length(dfcol)
          @inbounds ret[i] = ret[i] || dfcol[i]
        end
      end
    end
    return ret
end
#=
function (match::Match)
    return match.correspondence_properties[get_filtered_indices(match), :]
end
=#

function get_correspondences(match::Match; globalized::Bool=false, filtered::Bool=false, use_post::Bool=false, src_mesh = nothing, dst_mesh = nothing, manual_only=false)
    if filtered
      src_pts = match.src_points[:, get_filtered_indices(match, manual_only=manual_only)]; dst_pts = match.dst_points[:, get_filtered_indices(match, manual_only=manual_only)];
    else
      src_pts = copy(match.src_points); dst_pts = copy(match.dst_points);
    end
    if use_post
      src_pts = deform(src_pts, src_mesh);
      dst_pts = deform(dst_pts, dst_mesh);
    end
    if globalized
      globalize!(src_pts, get_offset(:match_image)); 
      globalize!(dst_pts, get_offset(:match_image));
    end
    filtered ? ret = (src_pts, dst_pts) : ret = (src_pts, dst_pts, get_filtered_indices(match, manual_only=manual_only));
    return ret;
end

### property handling
function get_dfs(dict::Dict, keytofetch)
    for key in keys(dict)
        if key == keytofetch
          return Base.get(dict, key, nothing)
        elseif typeof(Base.get(dict, key, nothing)) <: Dict
          if get_dfs(Base.get(dict, key, nothing), keytofetch) != nothing 
            return get_dfs(Base.get(dict, key, nothing), keytofetch) 
          end
        end
    end
    return nothing
end

function get_correspondence_properties(match::Match, property_name; filtered = false)
    if length(match.correspondence_properties) > 0
        ret = match.correspondence_properties[property_name];
        if length(ret) != 0
            ret = Array{typeof(ret[1])}(ret);
        end
        return filtered ? ret[get_filtered_indices(match)] : ret
    else
        return []
    end
end

function get_correspondence_properties(match::Match, fn::Function, args...; filtered = false)
    return filtered ? fn(match, args...)[get_filtered_indices(match)] : fn(match, args...)
end

### reviewing
function set_reviewed!(match::Match)
    match.properties[:review][:author] = author();
    return;
end

function is_reviewed(match::Match)
    return match.properties[:review][:author] == null_author()
end

function get_author(match::Match)
    author = match.properties[:review][:author]
    return author
end

function is_flagged(match::Match, filter_name = nothing)
    if filter_name == nothing return match.properties[:review][:flagged];
    else
    return haskey(match.properties[:review][:flags], filter_name); end
end

function flag!(match::Match, crit = nothing)
    match.properties[:review][:flagged] = true;
    if crit != nothing
    match.properties[:review][:flags][crit[1]] = crit[2:end];
      end
end

"""
Flag a match based on property criteria
"""
function flag!(match::Match, property_name, compare, threshold)
    attributes = get_correspondence_properties(match, property_name; filtered = true)
    inds_to_filter = find(i -> compare(i, threshold), attributes)
    if length(inds_to_filter) > 0
        flag!(match, (property_name, compare, threshold))
    end
end

function unflag!(match::Match, filter_name = nothing)
    if filter_name == nothing
    match.properties[:review][:flagged] = false;
    for key in keys(match.properties[:review][:flags])
      pop!(match.properties[:review][:flags], key);
    end
      else
      if haskey(match.properties[:review][:flags], filter_name)
      pop!(match.properties[:review][:flags], filter_name);
    end
      if length(match.properties[:review][:flags]) == 0
      match.properties[:review][:flagged] = false;
    end
      end
end

### helper methods
function get_ranges(pt, src_index, src_offset, src_image_size, dst_index, dst_offset, dst_image_size)
    get_ranges(pt, src_index, dst_index, PARAMS[:match][:block_r], PARAMS[:match][:search_r]);
end

function get_ranges(pt, src_index, src_offset, src_image_size, dst_index, dst_offset, dst_image_size, block_r::Int64, search_r::Int64)
    # convert to local coordinates in both src / dst images, and then round up to an integer
    src_pt = ceil.(Int64, pt);
    # if global_offsets
    rel_offset = src_offset - dst_offset;
    # else
    #   rel_offset = src_offset
    # end
    dst_pt = src_pt + rel_offset
    dst_pt = ceil.(Int64, dst_pt);

    block_range = -block_r:block_r;
    search_range = -(block_r+search_r):(block_r+search_r);

    src_range_full = src_pt[1] + block_range, src_pt[2] + block_range;
    dst_range_full = dst_pt[1] + search_range, dst_pt[2] + search_range;

    range_in_src = intersect(src_range_full[1], 1:src_image_size[1]), intersect(src_range_full[2], 1:src_image_size[2]);
    if length(range_in_src[1]) % 2 != 1
      range_in_src = range_in_src[1][1:end-1], range_in_src[2];
    end
    if length(range_in_src[2]) % 2 != 1
      range_in_src = range_in_src[1], range_in_src[2][1:end-1];
    end
    range_in_dst = intersect(dst_range_full[1], 1:dst_image_size[1]), intersect(dst_range_full[2], 1:dst_image_size[2]);

    src_pt_locs = findfirst(range_in_src[1] .== src_pt[1]), findfirst(range_in_src[2] .== src_pt[2]);
    dst_pt_locs = findfirst(range_in_dst[1] .== dst_pt[1]), findfirst(range_in_dst[2] .== dst_pt[2]);
    dst_pt_locs_full = findfirst(dst_range_full[1] .== dst_pt[1]), findfirst(dst_range_full[2] .== dst_pt[2]);

    if src_pt_locs[1] == 0 || src_pt_locs[2] == 0 || dst_pt_locs[1] == 0 || dst_pt_locs[2] == 0 
    return nothing
    end


    return src_index, range_in_src, [src_pt_locs[1], src_pt_locs[2]], dst_index, range_in_dst, dst_range_full, [dst_pt_locs[1], dst_pt_locs[2]], [dst_pt_locs_full[1], dst_pt_locs_full[2]], rel_offset;
end

function meanpad!(img, val)
  @fastmath avg = mean(img)
  @fastmath zero_entries = img .== val
  @inbounds img[zero_entries] = avg
  #=
  running_sum = 0.0;
  count = 0;
  zero_entries = Array{Int64, 1}();
    @simd for i in 1:length(img)
      @inbounds px = img[i];
      if px != 0
        @fastmath running_sum += px;
        @fastmath count += 1;
      end
      if px == 0
        push!(zero_entries, i)
      end
    end
    if length(zero_entries) == 0 return end
    @fastmath nonzero_mean = running_sum / count
    if eltype(img) != Float64
        @fastmath nonzero_mean = round(eltype(img), nonzero_mean)
    end
    for i in zero_entries
      @inbounds px = img[i];
      if px == 0
        @fastmath @inbounds img[i] = nonzero_mean;
#=      if i%2 != 0
        @fastmath @inbounds img[i] = nonzero_mean;
      else
        @fastmath @inbounds img[i] = nonzero_mean + eps::Float64;
        end =#
      end
    end
    =#
end

# if from_disk src_image / dst_image are indices
# function prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, scale, highpass_sigma; from_disk = false)
function prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, bandpass_sigmas; me = get_matchenv(src_range, dst_range_full, bandpass = bandpass_sigmas), meanpad = true, meanpad_val=0)
    clean!(me)
    indices_within_range = findin(dst_range_full[1], dst_range[1]), findin(dst_range_full[2], dst_range[2])
    @fastmath @inbounds me.src_patch_full[:] = view(src_image, src_range...)
    @fastmath @inbounds me.dst_patch_full[indices_within_range...] = view(dst_image, dst_range...)

    if meanpad
        meanpad!(me.src_patch_full, meanpad_val)
        meanpad!(me.dst_patch_full, meanpad_val)
    end

    if (bandpass_sigmas != (0,0)) & (bandpass_sigmas != [0,0])
        @inbounds me.src_patch_full[:] = convolve_Float64_planned(me.src_patch_full, me.kernel; crop = :same, padding = :mean)
        @inbounds me.dst_patch_full[:] = convolve_Float64_planned(me.dst_patch_full, me.kernel; crop = :same, padding = :mean)
    end

    return me.src_patch_full, me.dst_patch_full
end

function make_bandpass_kernel(lowpass_sigma, highpass_sigma)
    kernel_l = Images.Kernel.gaussian(lowpass_sigma)
    kernel_h = Images.Kernel.gaussian(highpass_sigma)
    oset = kernel_l.offsets[1]+1

    for j in oset:-oset
      for i in oset:-oset
        kernel_h[i, j] -= kernel_l[i,j]
      end
    end
    kernel = -kernel_h.parent
    return kernel
end

"""
Template match two image patches to produce point pair correspondence
"""
function get_match(pt, ranges, src_image, dst_image, scale = 1.0, bandpass_sigmas = (0, 0), full = false, meanpad = true, meanpad_val = 0)
    src_index, src_range, src_pt_loc, dst_index, dst_range, dst_range_full, dst_pt_loc, dst_pt_loc_full, rel_offset = ranges;
#=  if sum(src_image[src_range[1], first(src_range[2])]) == 0 && sum(src_image[src_range[1], last(src_range[2])]) == 0 &&
            sum(src_image[first(src_range[1]), src_range[2]]) == 0 && sum(src_image[last(src_range[1]), src_range[2]]) == 0 return nothing end
    if sum(dst_image[dst_range[1], first(dst_range[2])]) == 0 && sum(dst_image[dst_range[1], last(dst_range[2])]) == 0 &&
            sum(dst_image[first(dst_range[1]), dst_range[2]]) == 0 && sum(dst_image[last(dst_range[1]), dst_range[2]]) == 0 return nothing end=#

    #see if any of the edges in the template are fully padded
    if sum(src_image[src_range[1], first(src_range[2])]) == 0 return nothing end
    if sum(src_image[src_range[1], last(src_range[2])]) == 0 return nothing end
    if sum(src_image[first(src_range[1]), src_range[2]]) == 0 return nothing end
    if sum(src_image[last(src_range[1]), src_range[2]]) == 0 return nothing end


    
    if sum(src_image[src_range[1], round(Int64,median(src_range[2]))]) == 0 return nothing end
    if sum(src_image[round(Int64,median(src_range[1])), src_range[2]]) == 0 return nothing end
    dst_quart_range_i = linspace(dst_range[1][1], dst_range[1][end], 5)
    dst_quart_range_j = linspace(dst_range[2][1], dst_range[2][end], 5)
    cent_sum = 0;
    cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]):round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[2])]);
    cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]):round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[4])]);
    cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]), round(Int64, dst_quart_range_j[2]):round(Int64, dst_quart_range_j[4])]);
    cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[2]):round(Int64, dst_quart_range_j[4])]);
    if cent_sum == 0 return nothing end
    

#   if sum(dst_image[dst_range[1], round(Int64,linspace(dst_range[2][1], dst_range[2][end], 5)[2])]) == 0 return nothing end
#   if sum(dst_image[dst_range[1], round(Int64,linspace(dst_range[2][1], dst_range[2][end], 5)[4])]) == 0 return nothing end
#   if sum(dst_image[round(Int64,linspace(dst_range[1][1], dst_range[1][end], 5)[2]), dst_range[2]]) == 0 return nothing end
#   if sum(dst_image[round(Int64,linspace(dst_range[1][1], dst_range[1][end], 5)[4]), dst_range[2]]) == 0 return nothing end

    correspondence_properties = DataFrame();
    correspondence_properties[:ranges_src_pt_loc] = [src_pt_loc];
    correspondence_properties[:ranges_src_range] = src_range;
    correspondence_properties[:ranges_dst_pt_loc] = [dst_pt_loc];
    correspondence_properties[:ranges_dst_range] = dst_range;
    correspondence_properties[:ranges_dst_range_full] = dst_range_full;
    correspondence_properties[:ranges_scale] = scale;

#=
    src_range = ceil(Int64, scale * first(src_range[1])) : ceil(Int64, scale * last(src_range[1])), ceil(Int64, scale * first(src_range[2])) : ceil(Int64, scale * last(src_range[2]))
    dst_range = ceil(Int64, scale * first(dst_range[1])) : ceil(Int64, scale * last(dst_range[1])), ceil(Int64, scale * first(dst_range[2])) : ceil(Int64, scale * last(dst_range[2]))
    dst_range_full = ceil(Int64, scale * first(dst_range_full[1])) : ceil(Int64, scale * last(dst_range_full[1])), ceil(Int64, scale * first(dst_range_full[2])) : ceil(Int64, scale * last(dst_range_full[2]))
    src_pt_loc = ceil(Int64, scale * src_pt_loc);
    dst_pt_loc = ceil(Int64, scale * dst_pt_loc);
    dst_pt_loc_full = ceil(Int64, scale * dst_pt_loc_full);


    # padding
    if dst_range != dst_range_full
        indices_within_range = findin(dst_range_full[1], dst_range[1]), findin(dst_range_full[2], dst_range[2])
        intersect_img = dst_image[dst_range...];
        avg = mean(intersect_img);
        if isnan(avg) return nothing; end
        avg = round(eltype(dst_image), avg);
        padded_img = fill(avg, length(dst_range_full[1]), length(dst_range_full[2]));
        padded_img[indices_within_range...] = intersect_img;
        xc = normxcorr2_preallocated(src_image[src_range[1], src_range[2]], padded_img);
    else
    xc = normxcorr2_preallocated(src_image[src_range[1], src_range[2]], dst_image[dst_range[1], dst_range[2]]);
    end
    =#
    #tic()
    src_patch, dst_patch = prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, bandpass_sigmas, meanpad = meanpad, meanpad_val = meanpad_val)
    inputs_have_zeros = (0 in src_patch) | (0 in dst_patch)

    #if (pp == nothing) return nothing end;
    #   prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full, scale, highpass_sigma)
    xc::Array{Float64,2} = normxcorr2_preallocated(src_patch, dst_patch; shape = full ? "full" : "valid");

    r_max = maximum(xc)
    if isnan(r_max) 
        return nothing 
    end;

    function compute_dv(xc, r_val)
        #   if r_max > 1.0 println("rounding error") end
        ind = findfirst(r_max .== xc)
        i_max, j_max = ind2sub(xc, ind)
        if i_max == 0 
            i_max = size(xc, 1)
        end

        i_max_int = i_max
        j_max_int = j_max

        if 2 < i_max < size(xc, 1) - 1 && 2 < j_max < size(xc, 2) - 1
            xc_w = sum(xc[i_max-2:i_max+2, j_max-2:j_max+2])
            xc_i = sum(xc[i_max+2, j_max-2:j_max+2]) * 2 + sum(xc[i_max+1, j_max-2:j_max+2]) - sum(xc[i_max-1, j_max-2:j_max+2]) - sum(xc[i_max-2, j_max-2:j_max+2]) * 2
            xc_j = sum(xc[i_max-1:i_max+1, j_max+1]) - sum(xc[i_max-1:i_max+1, j_max-1])
            i_max += xc_i / xc_w 
            j_max += xc_j / xc_w 
            if isnan(i_max) || isnan(j_max) 
              println(xc_i)
              println(xc_j)
              println(xc_w)
              return nothing
            end
        end


        @fastmath @inbounds begin
        if scale != 1.0
            #length of the dst_range_full is always odd, so only need to care about the oddity / evenness of the source to decide the size of the xc in full
            if full
              xc_i_len = length(dst_range_full[1]) + length(src_range[1]) - 1
              xc_j_len = length(dst_range_full[2]) + length(src_range[2]) - 1
            else
                if isodd(length(src_range[1])) 
                  xc_i_len = length(dst_range_full[1]) - length(src_range[1]) + 1
                else
                  xc_i_len = length(dst_range_full[1]) - length(src_range[1]) 
                end
                if isodd(length(src_range[2])) 
                  xc_j_len = length(dst_range_full[2]) - length(src_range[2]) + 1
                else
                  xc_j_len = length(dst_range_full[2]) - length(src_range[2]) 
                end
            end
            xc_i_locs = linspace(1, xc_i_len, size(xc, 1))
            xc_j_locs = linspace(1, xc_j_len, size(xc, 2))
            i_max = (i_max - 1) * step(xc_i_locs) + 1
            j_max = (j_max - 1) * step(xc_j_locs) + 1
        end
        if full
            di = Float64(i_max + src_pt_loc[1] - dst_pt_loc_full[1] - length(src_range[1]))
            dj = Float64(j_max + src_pt_loc[2] - dst_pt_loc_full[2] - length(src_range[2]))
        else
            di = Float64(i_max - 1 + src_pt_loc[1] - dst_pt_loc_full[1])
            dj = Float64(j_max - 1 + src_pt_loc[2] - dst_pt_loc_full[2])
        end
        end
        return di, dj, i_max_int, j_max_int
    end

    di, dj, i_max, j_max = compute_dv(xc, r_max)

    function recompute_r_without_zeros(src_patch, dst_patch, di, dj)
        dst_bbox = BoundingBox(1, 1, size(dst_patch)...)
        src_offset = dst_pt_loc_full[1] - src_pt_loc[1] + di,
                        dst_pt_loc_full[2] - src_pt_loc[2] + dj
        src_bbox = BoundingBox(src_offset..., size(src_patch)...)
        dst_slice = bb_to_slice(dst_bbox - src_bbox)
        
        src_bbox = BoundingBox(1, 1, size(src_patch)...)
        dst_offset = src_pt_loc[1] - dst_pt_loc_full[1] - di,
                        src_pt_loc[2] - dst_pt_loc_full[2] - dj
        dst_bbox = BoundingBox(dst_offset..., size(dst_patch)...)
        src_slice = bb_to_slice(dst_bbox - src_bbox)
        src_sub = src_patch[src_slice...]
        dst_sub = dst_patch[dst_slice...]
        mask = (src_sub .!= 0) & (dst_sub .!= 0)
        return cor(src_sub[mask], dst_sub[mask])
    end

    if inputs_have_zeros
        r_max = recompute_r_without_zeros(src_patch, dst_patch, di, dj)
    end

    correspondence_properties[:patches_src_normalized_dyn_range] = (maximum(view(src_image, src_range...)) - minimum(view(src_image, src_range...))) / typemax(eltype(src_image));
    correspondence_properties[:patches_src_kurtosis] = kurtosis(view(src_image, src_range...));
    correspondence_properties[:patches_dst_kurtosis] = kurtosis(view(dst_image, dst_range...));
    correspondence_properties[:xcorr_r_max] = r_max;

    function compute_r_diff(xc_orig, rad, r_max_original, i_max, j_max)
        xc = copy(xc_orig)
        range_i = intersect(1:size(xc, 1), (-rad:rad) + i_max)
        range_j = intersect(1:size(xc, 2), (-rad:rad) + j_max)
        @inbounds xc[range_i, range_j] = -Inf
        r_max_new = maximum(xc)
        if inputs_have_zeros
            di, dj, i_max, j_max = compute_dv(xc_orig, r_max_new)
            r_max_new = recompute_r_without_zeros(src_patch, dst_patch, di, dj)
        end
        return r_max_original - r_max_new
    end

    for rad in [5, 10, 15]
        correspondence_properties[Symbol(string("xcorr_delta_", rad))] = compute_r_diff(xc, rad, r_max, i_max, j_max)
    end
    for beta in [0.5, 0.75, 0.95]
        correspondence_properties[Symbol(string("xcorr_sigma_", beta))] = sigma(xc, beta) / scale;
    end
    correspondence_properties[:vects_dv] = [[di, dj]];
    correspondence_properties[:vects_norm] = norm([di, dj]);
    #correspondence_properties[:posts] = Dict{Any, Any}();

    return vcat(pt + rel_offset + [di, dj], correspondence_properties);
end

function Base.filter!(match::Match, priority, function_name, compare, threshold, vars...)
    # attributes = get_correspondence_properties(match, property_name)
    attributes = eval(function_name)(match, vars...)
    compare = eval(compare)
    filter_col = fill(false, count_correspondences(match))
    @inbounds for i in 1:length(attributes)
      filter_col[i] = compare(attributes[i], threshold)
    end
    match.filters[Symbol(tuple(priority, function_name, compare, threshold, vars...))] = filter_col;
    return sum(filter_col)
    #=
    if attributes == nothing return 0; end
    inds_to_filter = find(i -> compare(i, threshold), attributes);
    type_name = function_name
    if length(vars) > 0
        type_name = vars[1]
    end 
    push!(match.filters, Dict{Any, Any}(
                :author => author(),
                :type     => type_name,
                :threshold => threshold,
                :rejected  => inds_to_filter,
                :function => function_name
                  ));
    #println("$(length(inds_to_filter)) / $(count_correspondences(match)) rejected.");
    return length(inds_to_filter);
    =#
end

function Base.filter!(match::Match, filter::Tuple)
    return filter!(match, filter...)
end

function get_residual_norms_post(match, ms)
    src_pts_after, dst_pts_after, filtered = get_globalized_correspondences_post(ms, findfirst(match_in_ms -> match_in_ms.src_index == match.src_index && match_in_ms.dst_index == match.dst_index, ms.matches));
    norms = zeros(eltype(dst_pts_after), size(dst_pts_after, 2))
    @fastmath @inbounds begin
      @simd for i in 1:size(dst_pts_after, 2)
      d1 = dst_pts_after[1,i] - src_pts_after[1,i]
      d2 = dst_pts_after[2,i] - src_pts_after[2,i]
      norms[i] = sqrt(d1^2+d2^2)
    end
        end
    return norms
end

function check(match::Match, function_name, compare, threshold, vars...)
     return eval(compare)(eval(function_name)(match, vars...), threshold)
end

function check!(match::Match, crits) 
    unflag!(match);
    for crit in crits
        if check(match, crit...) 
            flag!(match, crit);
        end
    end
    return is_flagged(match);
end

### ADD MANUAL FILTER
# function filter_manual!(match::Match, inds_to_filter; filtertype="manual")
#   push!(match.filters, Dict{Symbol, Any}(
#               :author   => author(),
#               :type     => filtertype,
#               :rejected  => inds_to_filter
#                 ));
#   return;
# end

function clear_filters!(match::Match; filtertype=nothing)
    # match.filters = match.filters[setdiff(1:length(match.filters), find(filter -> filter==filtertype, MeshSet(match.filters)))]
    if filtertype == nothing    
        delete!(match.filters, 1:size(match.filters,2))
    end
end

function undo_filter!(match::Match)
    if length(match.filters) > 0
        pop!(match.filters);
    end
end

"""
Create matches between two meshes

Args:
    * src_mesh (source image)
    * dst_mesh (destination image)

Output:
    * List of Match objects for each pair of src & dst mask components

Description:
    * Load mesh images & associated masks (indexed image)
    * Calculate the number of mask components for both src & dst
    * Create Match object for each pair of mask components
"""
function get_matches{T}(src_mesh::Mesh{T}, dst_mesh::Mesh{T})
    src_index = get_index(src_mesh); 
    dst_index = get_index(dst_mesh);
    src_image = get_image(src_index, :match_image, mip=get_mip(:match_image));  
    dst_image = get_image(dst_index, :match_image, mip=get_mip(:match_image));
    if use_roi_mask()
        src_roi = get_image(src_index, :roi_mask, mip=get_mip(:match_image), input_mip=get_mip(:roi_mask));
        dst_roi = get_image(dst_index, :roi_mask, mip=get_mip(:match_image), input_mip=get_mip(:roi_mask));
        roi_value = get_mask_value(:roi_mask)
        unsafe_mask_image!(src_image, src_roi, roi_value, src_image, keep_id=true)
        unsafe_mask_image!(dst_image, dst_roi, roi_value, dst_image, keep_id=true)
    end
    matches = Array{Match,1}()
    if use_defect_split()
      src_image_sub = deepcopy(src_image)
      dst_image_sub = deepcopy(dst_image)
      println("Getting masks for $(src_index) & $(dst_index):")
      src_mask = get_image(src_index, :defect_split, mip=get_mip(:match_image), input_mip=get_mip(:defect_split));
      dst_mask = get_image(dst_index, :defect_split, mip=get_mip(:match_image), input_mip=get_mip(:defect_split));
      println("Compiling mask components")
      src_subsections = Array{IMG_ELTYPE, 1}(unique(src_mask))
      dst_subsections = Array{IMG_ELTYPE,1}(unique(dst_mask))
      defect_value = get_mask_value(:defect_split)
      src_subsections = filter!(x -> x != defect_value, src_subsections)
      dst_subsections = filter!(x -> x != defect_value, dst_subsections)
      src_mesh.properties[:subsections] = src_subsections / SPLIT_MESH_BASIS
      dst_mesh.properties[:subsections] = dst_subsections / SPLIT_MESH_BASIS
      println("src_subsections: $src_subsections")
      println("dst_subsections: $dst_subsections")
      for ss in src_subsections
        for ds in dst_subsections
          src_id = src_index + ss/SPLIT_MESH_BASIS
          dst_id = dst_index + ds/SPLIT_MESH_BASIS
          unsafe_copy_image!(src_image, src_image_sub)
          unsafe_mask_image!(src_image, src_mask, ss, src_image_sub, keep_id=true)
          unsafe_copy_image!(dst_image, dst_image_sub)
          unsafe_mask_image!(dst_image, dst_mask, ds, dst_image_sub, keep_id=true)
          push!(matches, Match(src_mesh, dst_mesh, 
                                          src_id, dst_id, 
                                          src_image_sub, 
                                          dst_image_sub))
        end
      end
      src_image_sub = 0
      src_image_sub = 0
      dst_image_sub = 0
      dst_image_sub = 0
    else
      push!(matches, Match(src_mesh, dst_mesh, src_index, dst_index, 
                                                  src_image, dst_image))
    end
    gc();
    gc();
    return matches
end

function Match{T}(src_mesh::Mesh{T}, dst_mesh::Mesh{T}, 
                src_index=get_index(src_mesh), dst_index=get_index(dst_mesh), 
                src_image=get_image(src_mesh), dst_image=get_image(dst_mesh))
    println("Matching $(src_index) -> $(dst_index):")
    image_size = get_image_size(src_mesh);
    image_offset = get_offset(src_mesh);

    print("computing ranges:")
    @time ranges = map(get_ranges, columnviews(src_mesh.src_nodes), repeated(src_index), repeated(image_offset), repeated(image_size), repeated(dst_index), repeated(image_offset), repeated(image_size), repeated(PARAMS[:match][:block_r]), repeated(PARAMS[:match][:search_r]));
    ranged_inds = find(i -> i != nothing, ranges);
    ranges = ranges[ranged_inds];
    print("    ")

    println("$(length(ranged_inds)) / $(size(src_mesh.src_nodes, 2)) nodes to check.")

    if length(ranged_inds) != 0
        ranges = Array{typeof(ranges[1]), 1}(ranges);
    end


    print("computing matches:")
    print("    ")
    use_full_conv = false
    meanpad = true
    @time dst_allpoints = pmap(get_match, columnviews(src_mesh.src_nodes[:,ranged_inds]), 
                                    ranges, 
                                    repeated(src_image), 
                                    repeated(dst_image), 
                                    repeated(get_scale(:match_image)), 
                                    repeated(PARAMS[:match][:bandpass_sigmas]), 
                                    repeated(use_full_conv), 
                                    repeated(meanpad), 
                                    repeated(get_mask_value(:match_image))) 

    matched_inds = find(i -> i != nothing, dst_allpoints);

    #=for ap in dst_allpoints
      if typeof(ap) == RemoteException
        println(ap)
      end
    end=#

    src_points = copy(src_mesh.src_nodes[:,ranged_inds[matched_inds]]);
    dst_points = similar(src_points)
    for j in 1:size(dst_points, 2)
      @inbounds dst_points[1,j] = dst_allpoints[matched_inds[j]][1]
      @inbounds dst_points[2,j] = dst_allpoints[matched_inds[j]][2]
    end
    correspondence_properties = vcat([dst_allpoints[ind][3] for ind in matched_inds]...)

    filters = DataFrame();
    properties = Dict{Symbol, Any}(
        :params => PARAMS,
        :review => Dict{Symbol, Any}(
                :flagged => false,
                :flags => Dict{Symbol, Any}(),
                :author => null_author()
                ) 
            );

#   @everywhere init_Match();
    @everywhere gc();

    return Match{T}(src_index, dst_index, src_points, dst_points, correspondence_properties, filters, properties);
end
