function montage(firstindex::Index, lastindex::Index)
  for index in get_index_range(montaged(firstindex[1:2]...), montaged(lastindex[1:2]...))
    ms = MeshSet(index)
    if is_flagged(ms)
      render_montaged(ms; render_full=false, render_review=true, flagged_only=true)
    else
      render_montaged(ms; render_full=true, render_review=false)
    end
  end
end

function montage_migrate(firstindex::Index, lastindex::Index)
  for index in get_index_range(montaged(firstindex[1:2]...), montaged(lastindex[1:2]...))
    ms = load(index)
    migrate!(ms)
    check!(ms)
    save(ms)
    if is_flagged(ms)
      render_montaged(ms; render_full=true, render_review=true, flagged_only=true)
    else
      render_montaged(ms; render_full=true, render_review=false)
    end
  end
end


# function montage_section(wafer_num, n)
#   @time Ms, images = load_section(PREMONTAGED_OFFSETS, wafer_num, n);
#   @time add_all_matches!(Ms, images);
#   @time solve_meshset!(Ms);
#   save(Ms);
#   images = 0;
# end

# function montage_sections(wafer_num, k::UnitRange{Int64})
#   optimize_all_cores(PARAMS_MONTAGE);
#   for n in k
#     @time montage_section(wafer_num, n);
#   end
# end

function align_stack(first_wafer_num, first_sec_num, last_wafer_num, last_sec_num)
  println("Elastically aligning $first_wafer_num, $first_sec_num -> ... -> $last_wafer_num, $last_sec_num:")
  first_index = (first_wafer_num, first_sec_num, PREALIGNED_INDEX, PREALIGNED_INDEX);
  last_index = (last_wafer_num, last_sec_num, PREALIGNED_INDEX, PREALIGNED_INDEX);
  @time Ms = make_stack(first_index, last_index);
  for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index;
    @time b = Ms.meshes[i+1].index;
    @time add_pair_matches_reflexive!(Ms, a, b);
  end
  save(Ms)
  # solve_meshset!(Ms);
  # save(Ms);
  return Ms;
end

function prealign_section(src_wafer_num, src_sec_num)
  src_index = (src_wafer_num, src_sec_num, MONTAGED_INDEX, MONTAGED_INDEX);
  dst_index = find_preceding(src_index);
  if dst_index == NO_INDEX return Void; end
  println("Prealigning $src_index -> $dst_index:")
  @time images = affine_load_section_pair(src_index, dst_index)
  @time Ms = make_stack(dst_index, src_index);
  @time add_pair_matches_with_thumbnails!(Ms, src_index, dst_index, images);
  affine_solve_meshset!(Ms);
  save(Ms);
  return Ms;
end

function prealign_stack(first_wafer_num, first_sec_num, last_wafer_num, last_sec_num)
  println("Prealigning $first_wafer_num, $first_sec_num -> ... -> $last_wafer_num, $last_sec_num:")
end

#### ASSUMES HOMOGENEOUS TILE SIZE
#function premontage(wafer_range::UnitRange{Int64})
function premontage(wafer, start)
#  for wafer in wafer_range
    premontaged_path = get_path(premontaged(wafer, 1))
    dir,name = splitdir(premontaged_path)
    tiles = sort_dir(dir, ".h5");
    tiles = filter(x->contains(x,"Tile"), tiles)
    tile_indices = sort(map(parse_name, tiles))
    tile_indices = filter(x->x[1] == wafer, tile_indices)
    section_range = start:tile_indices[end][2]
  for sec in section_range
    premontaged_path = get_path(premontaged(wafer, sec))
    dir,name = splitdir(premontaged_path)
    println("$wafer, $sec")
    tiles = sort_dir(dir, ".h5");
    tiles = filter(x->contains(x,"Tile"), tiles)
    tile_indices = sort(map(parse_name, tiles))
    tile_indices = filter(x->x[1:2] == (wafer,sec), tile_indices)

    if length(tile_indices) == 0 continue; end

    fixed_indices = similar(tile_indices, 0)
    images = map(get_image, tile_indices)
    images = map(Images.restrict, images);
    images = map(Images.restrict, images);

    #fix the first one
    update_offset(tile_indices[1], [0, 0], [8000,8000]);
    push!(fixed_indices, tile_indices[1])
    oset = Points(length(tile_indices))
    oset[1] = [0,0];

    sigma = [1,1]

    for index in tile_indices[2:end]
	target_index = fixed_indices[findfirst(this -> is_adjacent(index, this), fixed_indices)]

	xc = normxcorr2(images[findfirst(ind -> ind == index, tile_indices)], images[findfirst(ind -> ind == target_index, tile_indices)]; shape="full")
	xcd = xc - Images.imfilter_gaussian(xc, sigma);

	rm, ind = findmax(xcd)
	i_max, j_max = ind2sub(xcd, ind);
	i_diff = 4 * (i_max - median(1:size(xcd, 1)))
	j_diff = 4 * (j_max - median(1:size(xcd, 2)))

#    	update_offset(index, [i_diff, j_diff] + get_offset(target_index), [size(images[findfirst(ind -> ind == index, tile_indices)])...]);
    	update_offset(index, [i_diff, j_diff] + oset[findfirst(ind -> ind == target_index, tile_indices)], [8000,8000]);
    	oset[findfirst(ind -> ind == index, tile_indices)] = [i_diff, j_diff] + oset[findfirst(ind -> ind == target_index, tile_indices)]
    	push!(fixed_indices, index)
    end
# end
 end
end
#=
function premontage(wafer::Int, section_range::UnitRange{Int64})
  for sec in section_range
    index = (wafer, sec, 0, 0)
    overview_path = get_path(get_overview_index(index))
    dir,name = splitdir(overview_path)
    println(dir)
    tiles = sort_dir(dir, "tif");
    tiles = filter(x->contains(x,"Tile"), tiles)

    save_fused_img_to = name[1:end-4]"_fused.jpg"
    save_xcorr_img_to = name[1:end-4]"_xcorr.png"
    if cur_dataset == "zebrafish"   ##################
      scale = 0.05
    else  # piriform
      #scale = 0.07
      scale = 0.3
    end
    offsets, = tiles_to_overview(tiles, overview_path, scale; tile_img_dir = dir,
        save_fused_img_to = joinpath(PREMONTAGED_DIR, save_fused_img_to),
        save_xcorr_img_to = joinpath(PREMONTAGED_DIR, save_xcorr_img_to),
        show_review_imgs = false)

    offset_file = joinpath(PREMONTAGED_DIR, "registry_premontaged.txt")
    f = open(offset_file, "a")
    for pair in offsets
      line = join((pair[1], pair[2]..., 8000, 8000), " ")
      write(f, line, "\n")
    end
    close(f)
  end
end
=#
