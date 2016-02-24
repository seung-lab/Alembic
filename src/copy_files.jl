function move_files()
	for i in 124:168
		folder = "research/GABA/data/atlas/MasterUTSLdirectory/07122012S2/S2-W002/HighResImages_ROI1_W002_7nm_120apa/S2-W002_Sec$(i)_Montage"
		original_path = joinpath(homedir(), "seungmount", folder)
		new_path = joinpath(homedir(), "seungmount/Omni/alignment/datasets", folder)
		if isdir(original_path)
			tiles = readdir(original_path)
			for tile in tiles
				if tile[end-1:end] == "h5"
					original_file = joinpath(original_path, tile)
					new_file = joinpath(new_path, tile)
					println(original_file)
					println(new_file)
					cp(original_file, new_file, remove_destination=true)
				end
			end
		end
	end
end

"""
rename all review images to start with "review"
"""
function rename_review_images()
	review_names = [(MONTAGED_DIR, "seam"), (PREALIGNED_DIR, "thumb"), (ALIGNED_DIR, "thumb_imfuse")]
	for (dir, prefix) in review_names
		for fn in readdir(joinpath(dir, "review"))
			print(fn)
			if length(fn) > length(prefix)
				if fn[1:length(prefix)] == prefix
					new_fn = string("review", fn[length(prefix)+1:end])
					print("\t", new_fn)
					mv(joinpath(dir, "review", fn), 
							joinpath(dir, "review", new_fn); 
							remove_destination=true)
				end
			end
			print("\n")
		end
		print("\n")
	end
end
