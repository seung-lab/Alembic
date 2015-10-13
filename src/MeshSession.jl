
function montage_section(wafer_num, n)
  @time Ms, images = load_section(PREMONTAGED_OFFSETS, wafer_num, n);
  @time add_all_matches!(Ms, images);
  @time solve_meshset!(Ms);
  save(Ms);
  images = 0;
end

function montage_sections(wafer_num, k::UnitRange{Int64})
  optimize_all_cores(PARAMS_MONTAGE);
  for n in k
    @time montage_section(wafer_num, n);
  end
end

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
  solve_meshset!(Ms);
  save(Ms);
  return Ms;
end

function prealign_section(src_wafer_num, src_sec_num)
  src_index = (src_wafer_num, src_sec_num, MONTAGED_INDEX, MONTAGED_INDEX);
  dst_index = find_preceding(src_index);
  println("Prealigning $src_index -> $dst_index:")
  @time images = affine_load_section_pair(src_index, dst_index)
  @time Ms = make_stack(dst_index, src_index);
  @time add_pair_matches_with_thumbnails!(Ms, src_index, dst_index, images);
  stats(Ms)
  save(Ms);
  return Ms;
end

function prealign_stack(first_wafer_num, first_sec_num, last_wafer_num, last_sec_num)
  println("Prealigning $first_wafer_num, $first_sec_num -> ... -> $last_wafer_num, $last_sec_num:")
end

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
      scale = 0.07
    end
    offsets, = tiles_to_overview(tiles, overview_path, scale; tile_img_dir = dir,
        save_fused_img_to = joinpath(PREMONTAGED_DIR, save_fused_img_to),
        save_xcorr_img_to = joinpath(PREMONTAGED_DIR, save_xcorr_img_to),
        show_review_imgs = false)

    offset_file = joinpath(PREMONTAGED_DIR, "premontaged_offsets_tilescale.txt")
    f = open(offset_file, "a")
    for pair in offsets
      line = join((pair[1], pair[2]...), " ")
      write(f, line, "\n")
    end
    close(f)
  end
end
