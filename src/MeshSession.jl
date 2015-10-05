
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

function align_stack(wafer_num, k::UnitRange{Int64})
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, k);
  for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
    @time add_pair_matches!(Ms, a, b);
    @time add_pair_matches!(Ms, b, a); 
  end
  save(Ms)
  solve_meshset!(Ms);
  save(Ms);
end

function align_stack(wafer_num, k::UnitRange{Int64}, fixed_interval)
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, k, fixed_interval);
  for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
    @time add_pair_matches!(Ms, a, b);
    @time add_pair_matches!(Ms, b, a); 
  end
  save(Ms)
  solve_meshset!(Ms);
  save(Ms);
end

function align_batch_to_fixed(wafer_num, aligned, batch::UnitRange{Int64})
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, aligned, batch);
  for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
    @time add_pair_matches!(Ms, a, b);
    @time add_pair_matches!(Ms, b, a); 
  end
  save(Ms)
  solve_meshset!(Ms);
  save(Ms);
end

function prealign(wafer_num, dst, src)# k::UnitRange{Int64})
  @time Ms = affine_make_stack(MONTAGED_OFFSETS, wafer_num, dst, src);
  @time add_pair_matches!(Ms, src, dst);
  @time affine_solve_meshset!(Ms);
  save(Ms);
  return Ms;
end

function prealign(wafer_num_a, sec_num_a, wafer_num_b, sec_num_b)# k::UnitRange{Int64})
  if wafer_num_a != wafer_num_b println("No support for different wafers yet.")
    return
  end
  optimize_all_cores(PARAMS_PREALIGNMENT);
  for src in sec_num_a:sec_num_b
    dst = src - 1;
    @time Ms = affine_make_stack(MONTAGED_OFFSETS, wafer_num_a, dst, src, false);
    @time add_pair_matches!(Ms, src, dst);
    @time affine_solve_meshset!(Ms);
    save(Ms);
  end
end

function align_to_fixed(wafer_num, aligned, prealigned)
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, aligned, prealigned);
  for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
    @time add_pair_matches!(Ms, a, b);
    @time add_pair_matches!(Ms, b, a); 
  end
  save(Ms)
  solve_meshset!(Ms);
  save(Ms);
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
