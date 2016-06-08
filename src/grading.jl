function show_supervoxel_counts(fn, attr)
	gt = h5read(fn, attr)

	id = []
	for i = 1:size(gt, 3)
		push!(id, unique(gt[:,:,i]))
	end

	num_id = map(length, id)

	delta_num_id = []
	for i = 2:length(id)
	    push!(delta_num_id, length(setdiff(id[i-1], id[i])))
	end

	x = []
	y = []
	for i = 1:length(id)
		append!(x, ones(length(id[i]))*i)
		append!(y, id[i])
	end

	fig = figure("supervoxel statistics", figsize=(20,20))
	PyPlot.clf()
	subplot(211)
	title(fn)
	p = plot(1:length(num_id[1:end-1]), num_id[1:end-1])
	# title("supervoxel count by section")
	xlabel("section no.")
	ylabel("supervoxel count")

	p = plot(1:length(delta_num_id[1:end-1]), delta_num_id[1:end-1])
	# title("supervoxel id changes by section")
	# xlabel("section no.")
	# ylabel("count of supervoxel id differences to the next section")

	subplot(212)
	scatter(x,y)
	# title("supervoxel attendance per section")
	xlabel("section no.")
	ylabel("supervoxel id")
end

function show_ground_truth()
	fn = joinpath(homedir(), "seungmount/Omni/TracerTasks/Piriform_cortex_Corrections/merge_lbl.h5")
	show_supervoxel_counts(fn, "label")
end

function show_watershed_no_threshold()
	fn = joinpath(homedir(), "seungmount/research/Julimaps/datasets/piriform/5_finished/stacks/3x3/piriform_17000-18023_13000-14023_484-739_SUPERVOXELS.h5")
	show_supervoxel_counts(fn, "main")
end

function show_watershed_thresholded()
	fn = joinpath(homedir(), "seungmount/research/Julimaps/datasets/piriform/5_finished/stacks/3x3/piriform_17000-18023_13000-14023_484-739_SUPERVOXELS_0.3.h5")
	show_supervoxel_counts(fn, "main")
end

function show_realigned_no_threshold()
	fn = joinpath(homedir(), "seungmount/research/Julimaps/datasets/piriform/5_finished/stacks/3x3/piriform_17000-18023_13000-14023_484-739_REALIGNED_SUPERVOXELS.h5")
	show_supervoxel_counts(fn, "main")
end
