
function preprocess_data(set, slice)
	img = load(prealigned(set,slice));
	if size(img)[1]>21000
		img = img[5000:20000, 5000:20000]
		save(get_path(prealigned(set, slice)), img)
	end
	return img
end

function preprocess_all(dataset)
	return [ preprocess_data(x[1],x[2]) for x in dataset]
end

function construct_alignment(set, slice, skip = 1, blur = (0,0))
	PARAMS_ALIGNMENT["match"]["bandpass_sigmas"] = blur
	ms = MeshSet(prealigned(set, slice), prealigned(set, slice+skip), solve = false);
	inspect(ms, 1)
	stats(ms)
	return ms
end

function construct_all(dataset)
	[ construct_alignment(x[1],x[2]) for x in dataset]
end

function compute_error_rate(set, slice, skip = 1, blur = (0,0))
	PARAMS_ALIGNMENT["match"]["bandpass_sigmas"] = blur;
	ms = MeshSet(prealigned(set, slice), prealigned(set, slice+2), solve = false);
	x = count_rejected_correspondences(ms.matches[1]);
	y = count_correspondences(ms.matches[1]);
	return (x,y)
end

function compare_with_net(dataset)
	blurs = [(0,0), (2.5, 12), (0,0)]	

	errors, overall = zeros(size(blurs)), zeros(size(blurs));
	for (set, slice) in dataset  				
		(x,y) = compute_error_rate(set, slice, 2)
		if x<200
			errors[1] += x
			overall[1] += y
			(x,y) = compute_error_rate(set, slice, 2, blurs[2])
			errors[2] += x
			overall[2] += y
			(x,y) = compute_error_rate(set, slice+1000, 2)
			errors[3] += x		
			overall[3] += y
		end
	end
	print("~~~~~~~~~~~~~~~~")
	print(errors./overall)
	return errors, overall
end


function find_best_bandpass(dataset)	
	blurs = [(0,0), (0,8), (0.5, 8), (1,8), (1.5, 8), (2, 8), (2.5, 8), (3, 8),
		 (0,12), (0.5, 12), (1,12), (1.5, 12), (2, 12), (2.5, 12), (3, 12),
		 (0,15), (0.5, 15), (1,15), (1.5, 15), (2, 15), (2.5, 15), (3, 15),
		 (0,20), (0.5, 20), (1,20), (1.5, 20), (2, 20), (2.5, 20), (3, 20),
		 (0,25), (0.5, 25), (1,25), (1.5, 25), (2, 25), (2.5, 25), (3, 25)]
	#blurs = [(0,0),(0.5,11), (0.5,12), (0.5, 0.13), (0.4,11)] #, (0.4,12), (0.4, 0.13), (0.6,11), (0.6,12), (0.6, 0.13)]
	
	errors, overall = zeros(size(blurs)), zeros(size(blurs));
	for (set, slice) in dataset  				
		(x,y) = compute_error_rate(set, slice)
		if x<500
			errors[1] += x
			overall[1] += y
			count = 2
			for blur in blurs
				if blur == (0,0)
					continue
				end
				(x,y) = compute_error_rate(set, slice, 2, blur)
				errors[count] += x
				overall[count] += y
				count += 1
			end
		end
	end
	print("~~~~~~~~~~~~~~~~")
	print(errors./overall)
	return errors, overall
end

function compute_r_variable(ms)
	error = zeros(100, 4)

	
	for i in collect(1:100)
		
		clear_filters!(ms.matches[1])
		#filter!(ms.matches[1], (1,:get_properties, >, 20-20*i/100, 0.95)) # 5-20
		filter!(ms.matches[1], (1,:get_properties, >, log(700*i), 0.75)) # 6-50
		#filter!(ms.matches[1], (1,:get_properties, >, 6*(i/100)+7, 0.50)) # 25-500
		#filter!(ms.matches[1], (4, :get_properties, <, i/100, "r_max"))
		# (4, :get_properties, <, i/100, "r_max") 
		r_filtered = get_rejected_indices(ms.matches[1])
		clear_filters!(ms.matches[1])
		filter!(ms.matches[1], (0,:get_properties,>,85,"norm"))
		dist_filtered = get_rejected_indices(ms.matches[1])

		if count_rejected_correspondences(ms.matches[1])>1000
			return error
		end
			
		a = setdiff(dist_filtered, r_filtered) # dist_but_not_r
		c = intersect(r_filtered, dist_filtered) # r_and_dist
		d = setdiff(r_filtered, c)             # Correct but in r
		b = count_correspondences(ms.matches[1]) - length(d) - length(c) - length(a) # Correct

		error[i, 1] = length(a); error[i, 2] = b; 
		error[i, 3] = length(c); error[i, 4] = length(d);
	end
	return error
end

function preprocess_meshes(dataset, skip = 1, blur= (0,0))
	PARAMS_ALIGNMENT["match"]["bandpass_sigmas"] = blur;
	for (set, slice) in dataset
		ms = MeshSet(prealigned(set, slice), prealigned(set, slice+skip), solve = false);
	end
end

using Plots
plotly()

function get_stats(mss)
	errors = zeros(100, 4)
	for ms in mss
		errors += compute_r_variable(ms)
	end
	return errors
end

function plot_rejection_error(errors, name="")
	x = linspace(0, 1, 100);
	y = zeros(100, 2)
	y[:,1] = (errors[:,3]+errors[:,4])./(errors[:,3]+errors[:,4]+errors[:,1]+errors[:,2])
	y[:,2] = (errors[:,1])./(errors[:,1]+errors[:,2])
	
	#trace1 = scatter(; x=x, y=y[:,2], name="rejections", mode="lines")
	#trace2 = scatter(; x=x, y=y[:,2], name="rate", mode="lines") 

	plot(y[:,1],y[:,2], linewidth=2, title="Rejection rate: "*name ) 
end

function plot_all(errors, name="")
	x = linspace(0,1,100)
	return plot(x,errors, linewidth=2, title="All: "*name ) 
end

function plot_bandpass_and_net_accross()
	slices = [(1, 2*x) for x in [1:47]]
	slices_net = [(1, 2*x+1000) for x in [1:47]]
	
	error_blur = get_stats(slices, 2)
	error_net = get_stats(slices_net, 2)
	
	plot_rejection_error(error_blur, "bandpass")
	plot_rejection_error(error_net, "net")
	plot_all(error_blur, "bandpass")
	plot_all(error_net, "net")
	
end

function load_meshsets(slices, skip = 1)
	return [load(MeshSet, (prealigned(i, j), prealigned(i, j+skip))) for (i, j) in slices]
end 

function compute_meshsets_error(mss, dist = 90, rvalue=0.15, kurtosis = 25, s_max_v = 3, s_mid_v = 15, s_low_v = 30)
	error = zeros(8)	
	for ms in mss
		if ms == nothing
			continue
		end
		s_max = filter_by_property(ms, (1,:get_properties, >, s_max_v, 0.95))
		s_mid = filter_by_property(ms, (2,:get_properties, >, s_mid_v, 0.75))
		s_low = filter_by_property(ms, (3,:get_properties, >, s_low_v, 0.50))
		r = filter_by_property(ms, (4, :get_properties,<, rvalue, "r_max"))
		k = filter_by_property(ms, (5,:get_properties, >, kurtosis, "src_kurtosis"))
		n = filter_by_property(ms, (0,:get_properties, >, dist,"norm"))

		clear_filters!(ms.matches[1])
		if size(n)[1]>1000
			continue
		end

		error[1] += size(n)[1]
		error[2] += size(r)[1]
		error[3] += size(k)[1]
		error[4] += size(s_max)[1]
		error[5] += size(s_mid)[1]
		error[6] += size(s_low)[1]
		error[7] += size(union(s_max, s_mid, s_low, r, n, k))[1]
		error[8] += count_correspondences(ms.matches[1])		
	end
	return error#./error[8]
end

function filter_by_property(ms, property)
	clear_filters!(ms.matches[1])
	filter!(ms.matches[1], property)
	ind = get_rejected_indices(ms.matches[1])
	return ind
end
	
