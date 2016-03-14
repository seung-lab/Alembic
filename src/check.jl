function check_montage_review_renders()
	jls_files = readdir(MONTAGED_DIR)
	seam_paths = readdir(joinpath(MONTAGED_DIR, "review"))
	seam_paths = [joinpath(MONTAGED_DIR, "review", fn) for fn in seam_paths]
	missed_seams = []
	path = joinpath(MONTAGED_DIR, "review", "missed_seams.txt")
	println("Saving missed seams:\n", path)
	for fn in jls_files
		if fn[end-2:end] == "jls"
			meshset = load(joinpath(MONTAGED_DIR, fn))
			bbs = find_boundingboxes(meshset)
			indices = [mesh.index for mesh in meshset.meshes]
			overlap_tuples = find_overlaps(bbs)
			for (k, (i,j)) in enumerate(overlap_tuples)
				seam_path = get_review_path(indices[i], indices[j])
				if !(seam_path in seam_paths)
					println("MISSED ", seam_path)
					push!(missed_seams, seam_path)
					writedlm(path, missed_seams)
				end
			end
		end
	end
	return missed_seams
end

function check_montage_review_copy()
end

function check_montage_jls_copy()
end

function check_tile_copy()
	omni_path = joinpath(homedir(), "seungmount/Omni/alignment/datasets/")
	research_path = joinpath(homedir(), "seungmount/")
	for (k, wafer_path) in WAFER_DIR_DICT
		omni_wafer = joinpath(omni_path, wafer_path)
		research_wafer = joinpath(research_path, wafer_path)
	end
end

function remove_tiffs_from_Omni_GABA()
end