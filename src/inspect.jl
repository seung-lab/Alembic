"""
Display matches index k as overlayed images with points
"""
function inspect_matches(meshset, k)
  matches = meshset.matches[k]
  indexA = matches.src_index
  indexB = matches.dst_index

  path = get_outline_filename("seam", indexB, indexA)
  if !isfile(path)
    path = get_outline_filename("seam", indexA, indexB)
  end
  img = h5read(path, "img")
  offset = h5read(path, "offset")
  scale = h5read(path, "scale")

  img = vcat(img, ones(UInt32, 400, size(img, 2)))

  params = meshset.properties["params"]["match"]
  params["scale"] = scale
  params["thumb_offset"] = offset
  params["src_index"] = matches.src_index
  params["dst_index"] = matches.dst_index
  params["src_offset"] = get_offset(meshset.meshes[find_mesh_index(meshset, indexA)])
  params["dst_offset"] = get_offset(meshset.meshes[find_mesh_index(meshset, indexB)])
  params["src_size"] = get_image_sizes(matches.src_index)
  params["dst_size"] = get_image_sizes(matches.dst_index)

  src_nodes, dst_nodes = get_globalized_correspondences(meshset, k)
  vectors = [hcat(src_nodes...); hcat(dst_nodes...)]
  src_nodes, dst_nodes = get_globalized_correspondences_post(meshset, k)
  vectors_t = [hcat(src_nodes...); hcat(dst_nodes...)]
  vectorsA = scale_matches(src_nodes, scale)
  vectorsB = scale_matches(dst_nodes, scale)
  vecs = offset_matches(vectorsA, vectorsB, offset)

  imview = view(img, pixelspacing=[1,1])
  big_vecs = change_vector_lengths([hcat(vecs[1]...); hcat(vecs[2]...)], 10)
  an_pts, an_vectors = show_vectors(imview..., big_vecs, RGB(0,0,1), RGB(1,0,1))
  return imview, vectors, vectors_t, copy(an_vectors.ann.data.lines), params
end