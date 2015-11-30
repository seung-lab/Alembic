

function waferpaths_to_dict(waferpath_filename)
    wdict = Dict()
    if isfile(waferpath_filename)
        warray = readdlm(waferpath_filename)
        for (idx, path) in zip(warray[:,1], warray[:,2])
            wdict[idx] = path
        end
    end
    return wdict
end

function parse_name(name::String)

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

    # aligned-section
    m = match(r"(\d*),(\d*)_aligned", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), ALIGNED_INDEX, ALIGNED_INDEX 
    end

    return ret
    
end

function parse_registry(path::String)
    registry = cell(0, 0)
    if isfile(path)
        file = readdlm(path)
        registry = cell(size(file, 1), size(file, 2) + 1) # name, index, dx, dy
        for i in 1:size(registry, 1)
            index = parse_name(file[i, 1])
            registry[i, 1] = get_name(index)
            registry[i, 2] = index
            for j in 3:size(registry, 2)
                registry[i, j] = file[i, j-1]
            end
        end
      end
    return registry
end

bucket_dir_path = "/mnt/bucket/labs/seung"
#=if isfile("bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("bucket_dir_path.txt"), '\n')
elseif isfile("../bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("../bucket_dir_path.txt"), '\n')
end
=#
datasets_dir_path = "research/Julimaps/datasets"
cur_dataset = "piriform"
#cur_dataset = "zebrafish"
affine_dir_path = "~"

premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"

wafer_filename = "wafer_paths.txt"
premontaged_registry_filename = "registry_premontaged.txt"
montaged_registry_filename = "registry_montaged.txt"
prealigned_registry_filename = "registry_prealigned.txt"
aligned_registry_filename = "registry_aligned.txt"

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, PREMONTAGED_OFFSETS, PREMONTAGE_DIR, ALIGNMENT_DIR

global BUCKET = bucket_dir_path
global AFFINE_DIR = affine_dir_path
global DATASET_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset)
global PREMONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path)
global MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path)
global PREALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path)
global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path)

waferpath_filename = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, wafer_filename)
global WAFER_DIR_DICT = waferpaths_to_dict(waferpath_filename)

premontaged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path, premontaged_registry_filename)
global REGISTRY_PREMONTAGED = parse_registry(premontaged_registry_path)

montaged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path, montaged_registry_filename)
global REGISTRY_MONTAGED = parse_registry(montaged_registry_path)

prealigned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path, prealigned_registry_filename)
global REGISTRY_PREALIGNED = parse_registry(prealigned_registry_path)

aligned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path, aligned_registry_filename)
global REGISTRY_ALIGNED = parse_registry(aligned_registry_path)

global GLOBAL_BB = BoundingBox(0,0,36000,36000)

show_plot = false
