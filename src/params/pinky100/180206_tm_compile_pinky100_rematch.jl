rematched_list = readdlm(joinpath(homedir(), ".julia/v0.6/Alembic/src/params/pinky100/180204_tm_rematch.txt"), Int64)

rematch_dir = joinpath(homedir(), ".alembic/remtach")
mkdir(rematch_dir)

# Copy expanded list of all files from gcloud
for i in rematched_list
	pre_fn = joinpath(rematch_dir, "$(i)_pre")
	pre_cmd = `gsutil -m ls -l "gs://neuroglancer/pinky100_v0/father_of_alignment/match/($i.*"`
	post_fn = joinpath(rematch_dir, "$(i)_post")
	post_cmd = `gsutil -m ls -l "gs://neuroglancer/pinky100_v0/father_of_alignment/match/*, $i*"`
	Base.run(pipeline(pre_cmd, pre_fn))
	Base.run(pipeline(post_cmd, post_fn))
end

# Remove last line from all the files (file count & download rate)
all_files = readdir(rematch_dir)
for f in all_files
	fn = joinpath(rematch_dir, f)
	rm_cmd = `sed -i '$ d' $fn`
	Base.run(rm_cmd)
end

# Concatenate all the files
full_files = [joinpath(rematch_dir, f) for f in readdir(rematch_dir)]
cat_cmd = `cat $full_files`
cat_fn = joinpath(homedir(), ".alembic/180205_tm_gsutil_rematch_files.txt")
Base.run(pipeline(cat_cmd, cat_fn))

# Load combined file list
gs_match_list = readdlm(cat_fn)

# Compile list of objects to keep and remove
remove_list = gs_match_list[[!(parse(i[9:10]) in [4, 5]) for i in gs_match_list[:,2]], :]
keep_list = gs_match_list[[(parse(i[9:10]) in [4, 5]) for i in gs_match_list[:,2]], :]
unique([i[1:10] for i in remove_list[:,2]])
unique([i[1:10] for i in keep_list[:,2]])

# Save out the remove list
remove_file_list = [string(remove_list[i,3], " ", remove_list[i,4]) for i in 1:size(remove_list,1)]
rm_fn = joinpath(homedir(), ".alembic/180205_tm_gsutil_remove_files.txt")
writedlm(rm_fn, remove_file_list)