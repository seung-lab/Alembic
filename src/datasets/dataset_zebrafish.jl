global BUCKET = "/home/ubuntu"
global DATASET = "datasets/zebrafish"
global ROI_FIRST = (2,33,0,0);
global ROI_LAST = (9,164,0,0);
global DATASET_RESOLUTION = [5,5,45]

function get_name_legacy(index)
    if is_overview(index)
      if index[1] < 10
        return string("MontageOverviewImage_W00", index[1], "_sec", index[2])
      else
        return string("MontageOverviewImage_W0", index[1], "_sec", index[2])
      end
    elseif is_montaged(index)
        return string(index[1], ",", index[2], "_montaged")
    elseif is_prealigned(index)
      	if is_subsection(index)
        	return string(index[1], ",", index[2], "_prealigned_", index[4])
	end
        return string(index[1], ",", index[2], "_prealigned")
    elseif is_aligned(index)
        return string(index[1], ",", index[2], "_aligned")
    elseif is_finished(index)
        return string(index[1], ",", index[2], "_finished")
    else
      if index[1] < 10
    return string("Tile_r", index[3], "-c", index[4], "_W00", index[1], "_sec", index[2])
  	else
    return string("Tile_r", index[3], "-c", index[4], "_W0", index[1], "_sec", index[2])
  end
    end
end

# function get_path()
# methods: 
#     
# extensions:
# Mesh.jl: get_path(mesh::Mesh)
function get_path_legacy(index, ext = ".h5")
    name = get_name_legacy(index)
    if is_overview(index)
      if index[1] < 10
        section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
      else
        section_folder = string("W0", index[1], "_Sec", index[2], "_Montage")
      end
        #path = joinpath(BUCKET, WAFER_DIR_PATH_DICT[index[1]], section_folder, string(name, ext))
        path = joinpath(PREMONTAGED_DIR_PATH, string(name, ext))
    elseif is_montaged(index)
        path = joinpath(MONTAGED_DIR_PATH, string(name, ext))
    elseif is_prealigned(index)
        path = joinpath(PREALIGNED_DIR_PATH, string(name, ext))
    elseif is_aligned(index)
        path = joinpath(ALIGNED_DIR_PATH, string(name, ext))
    elseif is_finished(index)
        path = joinpath(FINISHED_DIR_PATH, string(name, ext))
    else
#        section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        #path = joinpath(BUCKET, WAFER_DIR_PATH_DICT[index[1]], section_folder, string(name, ext))
        path = joinpath(PREMONTAGED_DIR_PATH, string(name, ext))
    end
    return path
end

function migrate_legacy()
  registry_premontaged_legacy_path = joinpath(PREMONTAGED_DIR_PATH, "registry_premontaged.txt")
  registry = readdlm(registry_premontaged_legacy_path)
  for i in 1:size(registry, 1)
	ind = parse_name_legacy(registry[i,1])
	if isfile(get_path_legacy(ind))
	run(`mv $(get_path_legacy(ind)) $(get_path(ind))`)
      end
	registry[i,1] = get_name(parse_name_legacy(registry[i,1]))
  end
  writedlm(get_registry_path((1,1,1,1)), registry)
  run(`mv $(joinpath(MONTAGED_DIR_PATH, "registry_montaged.txt")) $(get_registry_path(montaged(1,1)))`)
  run(`mv $(joinpath(PREALIGNED_DIR_PATH, "registry_prealigned.txt")) $(get_registry_path(prealigned(1,1)))`)
  run(`mv $(joinpath(ALIGNED_DIR_PATH, "registry_aligned.txt")) $(get_registry_path(aligned(1,1)))`)
end

function parse_name_legacy(name::String)

    ret = (0, 0, 0, 0)
    # singleton tile
    m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2])
    end

    # overview image
    m = match(r"MontageOverviewImage_S2-W00(\d*)_sec(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), OVERVIEW_INDEX, OVERVIEW_INDEX   
    end

    # montaged section
    m = match(r"(\d*),(\d*)_montaged", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), MONTAGED_INDEX, MONTAGED_INDEX   
    end

    # prealigned section
    m = match(r"(\d*),(\d*)_prealigned", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), PREALIGNED_INDEX, PREALIGNED_INDEX 
    end

    # prealigned_subsection
    m = match(r"(\d*),(\d*)_prealigned_(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), PREALIGNED_INDEX, parse(Int, m[3])
    end

    # aligned-section
    m = match(r"(\d*),(\d*)_aligned", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), ALIGNED_INDEX, ALIGNED_INDEX 
    end

    return ret
    
end


