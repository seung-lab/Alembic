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
# Stopped 160806 incomplete - see EverNote for error
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
# Stopped 160806 incomplete - see EverNote for error
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

# Started 160819 17:20 on seungworkstation11
function restart_upload_AIBS_google_ALL()
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
		src_path = joinpath(src_dir, dir)
		dst_path = joinpath(dst_dir, dir)
		run(`gsutil -m rsync -r $src_path $dst_path`)
	end
end

# Started 160819 17:20 on seungworkstation11
function restart_upload_AIBS_AWS_ALL()
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
		src_path = joinpath(src_dir, dir)
		dst_path = joinpath(dst_dir, dir)
		run(`aws s3 sync $src_path $dst_path`)
	end
end

function upload_ground_truth_to_google()
	dirs = ["vol01"]
	subdirs = ["Array{Any}_segmentation", "raw"]
	src_root = joinpath(homedir(), "seungmount/Omni/TracerTasks/AIBS_practice_234251S6R_01_01_aligned_01/ground_truth/")
	dst_root = "gs://s6r_ground_truth"
	for dir in dirs
		for subdir in subdirs
			src_path = joinpath(src_root, dir, subdir)
			dst_path = joinpath(dst_root, dir)
			run(`gsutil -m cp -r $src_path $dst_path`)
		end
	end
end

# Started 160819 17:20 on seungworkstation11
function upload_some_tiles_to_AWS()
	dst_dir = "s3://seunglab/datasets/pinky/1_premontaged/" # neuromancer-seung-import
	indices = get_index_range(premontaged(1,4001), premontaged(1,4021))
	for index in indices
		src_path = get_path(index)
		# dst_path = joinpath(dst_dir, string(get_name(index), ".h5"))
		Base.run(`aws s3 cp $src_path $dst_dir`)
	end
	src_path = get_registry_path(indices[1])
	# dst_path = joinpath(dst_dir, "registry.txt")
	Base.run(`aws s3 cp $src_path $dst_dir`)
end

function clean_up_aibs_google_cloud_dir_paths(fn)
	a = readdlm(fn)
	a = [join(a[i,:]...) for i in 1:size(a,1)]
	b = [split(i, '/')[2] for i in a]
	c = [length(i)>2 for i in b]
	d = b[convert(BitArray, c)]
	e = [i[end-2:end]=="tif" for i in d]
	f = d[convert(BitArray, e)]
	g = map(length, d)
	h = d[g .> 3]
	l = [join(split(i, "_")[2:end], "_")[1:end-4] for i in h]
	m = [join(split(i, "_")[1:end-2], "_") for i in l]
	n = [split(i, "_")[1] == "243774" for i in m]
	o = m[convert(BitArray, n)]
	p = unique(o)
	writedlm(fn, p)
end