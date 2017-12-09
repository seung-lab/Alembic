function get_correspondence_patches(m::Match, ind)
  src_z, dst_z = get_index(m);
  props = m.correspondence_properties[ind, :]

  # scale = props[:ranges_scale][1];
  bandpass_sigmas = m.properties[:params][:match][:bandpass_sigmas]
  src_pt = props[:ranges_src_pt_loc][1];
  dst_pt = props[:ranges_dst_pt_loc][1];
  # dst_range_full = props[:ranges_dst_range_full][1]

  function adjust_range(r, offset)
    return r[1] - offset[1], r[2] - offset[2]
  end

  mip = m.properties[:params][:match][:mip]
  offset = get_offset("src_image", mip=mip)
  src_range = props[:ranges_src_range][1]
  dst_range = props[:ranges_dst_range][1]

  src_range_global = tuple(adjust_range(src_range, offset)..., src_z)
  dst_range_global = tuple(adjust_range(dst_range, offset)..., dst_z)

  src_patch = get_image("src_image", src_range_global, mip=mip)
  dst_patch = get_image("src_image", dst_range_global, mip=mip)
  src_match = get_image("match_image", src_range_global, mip=mip)
  dst_match = get_image("match_image", dst_range_global, mip=mip)

  src_dummy_range = (1:size(src_match,1), 1:size(src_match,2))
  dst_dummy_range = (1:size(dst_match,1), 1:size(dst_match,2))

  # src_patch, dst_patch = prepare_patches(src_patch, dst_patch, 
  #                                         src_dummy_range, dst_dummy_range, 
  #                                         dst_dummy_range, bandpass_sigmas, 
  #                                         meanpad=true)
  src, dst = prepare_patches(src_match, dst_match, 
                                          src_dummy_range, dst_dummy_range, 
                                          dst_dummy_range, bandpass_sigmas, 
                                          meanpad=true)
  xc = normxcorr2_preallocated(src, dst; shape = "valid");
  dv = ceil(Int64, props[:vects_dv][1])
  src_patch = convert(Array{UInt8,2}, src_patch) 
  dst_patch = convert(Array{UInt8,2}, dst_patch)
  src_match = convert(Array{UInt8,2}, src_match)
  dst_match = convert(Array{UInt8,2}, dst_match)
  xc = xc_to_uint8(xc)
  src_offset = -src_pt+dv
  dst_offset = -dst_pt
  xc_offset = -dst_pt+src_pt

  @time save_image("src_patch", src_patch, src_offset, 0, mip=mip, cdn_cache=false)
  @time save_image("src_match", src_match, src_offset, 0, mip=mip, cdn_cache=false)
  @time save_image("dst_patch", dst_patch, dst_offset, 0, mip=mip, cdn_cache=false)
  @time save_image("dst_match", dst_match, dst_offset, 0, mip=mip, cdn_cache=false)
  @time save_image("xc", xc, xc_offset, 0, mip=mip, cdn_cache=false)
  return src_patch, dst_patch, src_match, dst_match, xc, src_pt, dst_pt, dv
end

function get_url()
  url = "https://neuromancer-seung-import.appspot.com/#!"
  url = string(url, "{'layers':")
  url = string(url, "{'xc':{'type':'image'_'source':")
  url = string(url, "'precomputed://$(get_path("xc"))'}")
  url = string(url, "_'dst_patch':{'type':'image'_'source':")
  url = string(url, "'precomputed://$(get_path("dst_patch"))'_")
  url = string(url, "'shader':'void main() {\nfloat val = toNormalized(getDataValue());\n if (val == 0.0) {\n emitRGBA(vec4(0.0));\n }\n else{\n emitRGBA(vec4(0.0, val, 0.0, 1.0));\n }\n}'}")
  url = string(url, "_'dst_match':{'type':'image'_'source':")
  url = string(url, "'precomputed://$(get_path("dst_match"))'_")
  url = string(url, "'shader':'void main() {\nfloat val = toNormalized(getDataValue());\n if (val == 0.0) {\n emitRGBA(vec4(0.0));\n }\n else{\n emitRGBA(vec4(0.0, val, 0.0, 1.0));\n }\n}'}")
  url = string(url, "_'src_patch':{'type':'image'_'source':")
  url = string(url, "'precomputed://$(get_path("src_patch"))'_")
  url = string(url, "'shader':'void main() {\nfloat val = toNormalized(getDataValue());\n if (val == 0.0) {\n emitRGBA(vec4(0.0));\n }\n else{\n emitRGBA(vec4(val, 0.0, 0.0, 1.0));\n }\n}'}")
  url = string(url, "_'src_match':{'type':'image'_'source':")
  url = string(url, "'precomputed://$(get_path("src_match"))'_")
  url = string(url, "'shader':'void main() {\nfloat val = toNormalized(getDataValue());\n if (val == 0.0) {\n emitRGBA(vec4(0.0));\n }\n else{\n emitRGBA(vec4(val, 0.0, 0.0, 1.0));\n }\n}'}")
  url = string(url, "}_")
  url = string(url, "'navigation':{'pose':{'position':")
  url = string(url, "{'voxelSize':[4_4_40]_'voxelCoordinates':[1024_1024_0]}}_")
  url = string(url, "'zoomFactor':16}_'layout':'xy'}")
  return url
end

function xc_to_uint8(xc)
  b = xc / maximum(xc)
  b[b .> 1] = 1
  b[b .< 0] = 0
  b = b * 254 + 1
  b[isnan(b)] = 0
  return round(UInt8, b)
end

