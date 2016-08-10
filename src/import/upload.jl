function upload_AIBS_google_z_index()
	# fn = "/media/tmacrina/667FB0797A5072D7/3D_align/mosaiced_images_160720_google_cloud_upload.csv" # FINISHED 160724 to allen-brain-institute
	# fn = "/media/tmacrina/667FB0797A5072D7/3D_align/mosaiced_images_160726_google_cloud_upload.csv" # FINISHED 160728 to allen-brain-institute
	fn = "/media/tmacrina/667FB0797A5072D7/3D_align/mosaiced_images_160729_google_cloud_upload_seung_import.csv" # FINISHED 160731 to seung-import
	meta = readdlm(fn, ',')
	src_dir = "/media/tmacrina/"
	# dst_dir = "pinky" # allen-brain-institute
	dst_dir = "pinky_raw" # neuromancer-seung-import
	for i in 1:size(meta,1)
		src_path = joinpath(src_dir, meta[i,2], "*")
		dst_path = joinpath("gs://", dst_dir, string(meta[i,1]), "")
		run(`gsutil -m cp -r $src_path $dst_path`)
	end
end

# Started 160803 16:22 on seungworkstation11
function upload_AIBS_google_ALL()
	dirs = ["05CD6FE95EABBD48",
			"2AFA162542D90D86",
			"3E68987D7020E00D",
			"416062E31672DE2C",
			"4BED39E032CF5004",
			"5C2D807B721C2F35",
			"61B81A1F318B9618",
			"667FB0797A5072D7"]
	subdirs_to_ignore = ["3D_align", "2D_alignment", "datasets", ".Trash-1753"]
	src_dir = "/media/tmacrina/"
	# dst_dir = "pinky" # allen-brain-institute
	dst_dir = "gs://pinky_raw" # neuromancer-seung-import
	for dir in dirs
		for subdir in readdir(joinpath(src_dir, dir))
			src_path = joinpath(src_dir, dir, subdir)
			if isdir(src_path) && !(subdir in subdirs_to_ignore)
				dst_path = joinpath(dst_dir, dir, subdir)
				run(`gsutil -m cp -r $src_path $dst_path`)
			end
		end
	end
end

# Started 160803 16:30 on seungworkstation11
function upload_AIBS_AWS_ALL()
	dirs = ["05CD6FE95EABBD48",
			"2AFA162542D90D86",
			"3E68987D7020E00D",
			"416062E31672DE2C",
			"4BED39E032CF5004",
			"5C2D807B721C2F35",
			"61B81A1F318B9618",
			"667FB0797A5072D7"]
	subdirs_to_ignore = ["3D_align", "2D_alignment", "datasets", ".Trash-1753"]
	src_dir = "/media/tmacrina/"
	# dst_dir = "pinky" # allen-brain-institute
	dst_dir = "s3://seunglab/datasets/pinky" # neuromancer-seung-import
	for dir in dirs
		for subdir in readdir(joinpath(src_dir, dir))
			src_path = joinpath(src_dir, dir, subdir)
			if isdir(src_path) && !(subdir in subdirs_to_ignore)
				dst_path = joinpath(dst_dir, dir, subdir)
				run(`aws s3 sync $src_path $dst_path`)
			end
		end
	end
end
