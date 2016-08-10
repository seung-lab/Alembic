function eval_filters(match::Match, filters, conjunction=false, meshset=nothing)

	inds_to_filter = Array{Any, 1}();
	thresholds = Array{Int64, 1}();

	for filter in filters
	if typeof(filter[1]) == Function
	attributes = get_properties(match, filter[1], meshset)
	else
	attributes = get_properties(match, filter[1]);
	end
	push!(inds_to_filter, find(i -> filter[2](i, filter[3]), attributes));
	push!(thresholds, filter[4]);
	end

	rejected_inds = get_rejected_indices(match);

	if conjunction == false
		if Base.|((thresholds .< map(length, inds_to_filter))...) filter_reject_match = true
		else filter_reject_match = false; end
	else
		if Base.&((thresholds .< map(length, inds_to_filter))...) filter_reject_match = true
		else filter_reject_match = false; end
	end

	inds_to_filter = union(inds_to_filter...)

	if length(rejected_inds) > 0 actual_reject_match = true;
	else actual_reject_match = false; end

	false_rejections = setdiff(inds_to_filter, rejected_inds)
	false_acceptances = setdiff(rejected_inds, inds_to_filter)
	common_rejections = intersect(rejected_inds, inds_to_filter)


	return length(false_rejections), length(false_acceptances), length(common_rejections), count_correspondences(match), filter_reject_match, actual_reject_match;
end

function eval_filters(ms::MeshSet, filters, range, conjunction=false)
	evals = map(eval_filters, ms.matches[range], repeated(filters), repeated(conjunction), repeated(ms))

	total_false_rej= 0;
	total_false_acc = 0;
	total_correct = 0;
	total_corresp = 0;

	match_false_rej = 0;
	match_false_acc = 0;
	match_correct = 0;

	for i in evals
	total_false_rej = total_false_rej + i[1];
	total_false_acc = total_false_acc + i[2];
	total_correct = total_correct + i[3];
	total_corresp = total_corresp + i[4];

	if i[5] == true && i[6] == true match_correct = match_correct + 1; end
	if i[5] == true && i[6] == false match_false_rej = match_false_rej + 1; end
	if i[5] == false && i[6] == true match_false_acc = match_false_acc + 1;
	println(findfirst(ind-> ind==i, evals), " was wrongly accepted by filters")
	end

	end

	total = total_false_acc + total_correct

	println("filtering by: $(filters...)")
#=
	println("false rejections: $(total_false_rej / total)")
	println("false non-rejections: $(total_false_acc / total)")
	println("correct rejections: $(total_correct / total)")
=#
	println("Per match:")
	println("precision: $(100 * match_correct / (match_false_rej + match_correct)) %")
	println("recall: $(100 * match_correct / (match_correct + match_false_acc)) %")
	println();

	println("total matches with issue: $(100 * (match_correct + match_false_acc) / (count_matches(ms))) %")
	println("workload reduced to: $(100 * (match_correct + match_false_rej) / (count_matches(ms))) %")
	println();

	println("Per correspondence:")
	println("precision: $(100 * total_correct / (total_false_rej + total_correct)) %")
	println("recall: $(100 * total_correct / (total_correct + total_false_acc)) %")
	return total_false_rej, total_false_acc, total_correct, total_corresp
end

function eval_filters_meshsets(mses, filters)
	evals = pmap(eval_filters, mses, repeated(filters))

	total_false_rej = 0;
	total_false_acc = 0;
	total_correct = 0;

	for i in evals
	total_false_rej = total_false_rej + i[1];
	total_false_acc = total_false_acc + i[2];
	total_correct = total_correct + i[3];
	total_corresp = total_corresp + i[4];
	end

	println("filtering by: $(filters...)")
	println("precision: $(100 * total_correct / (total_false_rej + total_correct)) %")
	println("recall: $(100 * total_correct / (total_correct + total_false_acc)) %")
	println("total correspondences: $total_corresp")
	println("total correct: $total_correct")
	println("total falsely rejected: $total_false_rej")
	println("total falsely accepted: $total_false_acc")
end















# if bias is a number in (0, 1) then it will fill the training data with that proportion being flagged correspondences
function make_training_data(ms::MeshSet, bias=0.0)

       X = Array{Float64, 2}(count_correspondences(ms), 4)
       Y = fill(-1, count_correspondences(ms))
       ranges = Array{UnitRange, 1}();
       current = 0

       for m in ms.matches
	range = (1:count_correspondences(m)) + current;

              X[range, 1] = get_properties(m, "norm")[:]';
              X[range, 2] = get_properties(m, "r_val")[:]';
              X[range, 3] = get_properties(m, "src_normalized_dyn_range")[:]';
       	   X[(1:count_correspondences(m)) + current, 4] = get_properties(m, "src_kurtosis")[:]';  
       	   Y[(collect(Int64, get_rejected_indices(m))) + current] = 1
	  push!(ranges, range)
       	  current = current + count_correspondences(m);
       end
	if bias != 0.0
		flagged_inds = find(i -> i == 1, Y);
		unflagged_inds = setdiff(1:length(Y), flagged_inds);
		X_f = X[flagged_inds, :];
		Y_f = Y[flagged_inds];
		X_uf = X[unflagged_inds, :];
		Y_uf = Y[unflagged_inds];
		length_uf = round(Int64, length(flagged_inds) * ((1 / bias) - 1), RoundUp);
		X = vcat(X_f, X_uf[1:length_uf, :]);
		Y = vcat(Y_f, Y_uf[1:length_uf]);
       		ranges = Array{UnitRange, 1}();
		push!(ranges, 1:length(Y))
	end
       return X, Y, ranges
end

function logreg_and_optimize(X, Y)
	logreg = BinomialLogReg(X, Y)
	model = optimize(logreg)
	return logreg, model;
end

function evaluate_filter(model, X, Y, ranges)
	Ybar = predict(model, X);

	println("Evaluating filter:")
	total_false_flagged = 0;
	total_false_passed = 0;
	total_correct_flagged = 0;
	total_corresp = length(Y);

	for ind in 1:length(Y)
		total_corresp = total_corresp + 1;
		if Y[ind] == Ybar[ind] && Y[ind] == 1 total_correct_flagged = total_correct_flagged + 1 end;
		if Y[ind] != Ybar[ind] && Y[ind] == -1 total_false_flagged = total_false_flagged + 1 end; 
		if Y[ind] != Ybar[ind] && Y[ind] == 1 total_false_passed = total_false_passed + 1 end; 
	end

	match_false_flagged = 0;
	match_false_passed = 0;
	match_correct_flagged = 0;
	match_total = 0;

	for range in ranges
		match_total = match_total + 1;
		if in(1, Y[range]) == in(1, Ybar[range]) && in(1, Y[range]) == true match_correct_flagged = match_correct_flagged + 1 end;
		if in(1, Y[range]) != in(1, Ybar[range]) && in(1, Y[range]) == false match_false_flagged = match_false_flagged + 1 end;
		if in(1, Y[range]) != in(1, Ybar[range]) && in(1, Y[range]) == true match_false_passed = match_false_passed + 1 end;
	end

	println("Per match:")
	println("precision: $(100 * match_correct_flagged / (match_false_flagged + match_correct_flagged)) %")
	println("recall: $(100 * match_correct_flagged / (match_correct_flagged + match_false_passed)) %")
	println();

	println("total matches with issue: $(100 * (match_correct_flagged + match_false_passed) / (match_total)) %")
	println("workload reduced to: $(100 * (match_correct_flagged + match_false_flagged) / (match_total)) %")
	println();

	println("Per correspondence:")
	println("precision: $(100 * total_correct_flagged / (total_false_flagged + total_correct_flagged)) %")
	println("recall: $(100 * total_correct_flagged / (total_correct_flagged + total_false_passed)) %")

	println("total correspondences with issue: $(100 * (total_correct_flagged + total_false_passed) / (total_corresp)) %")
	println("workload reduced to: $(100 * (total_correct_flagged + total_false_flagged) / (total_corresp)) %")
	println();
end
