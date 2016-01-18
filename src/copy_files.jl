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
