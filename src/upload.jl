function upload_AIBS_google()
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
