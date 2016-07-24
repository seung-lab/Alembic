function get_name(index::Index)
    if is_overview(index)
        if cur_dataset == "zebrafish"
            return string("MontageOverviewImage_W00", index[1], "_sec", index[2])
        else
            return string("MontageOverviewImage_S2-W00", index[1], "_sec", index[2])
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
    return string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
    end
end

# function get_path()
# methods: 
#     
# extensions:
# Mesh.jl: get_path(mesh::Mesh)
function get_path(index, ext = ".h5")
    name = get_name(index)
    if is_overview(index)
        if cur_dataset == "zebrafish"
            section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
        else
            section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        end
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ext))
    elseif is_montaged(index)
        path = joinpath(MONTAGED_DIR, string(name, ext))
    elseif is_prealigned(index)
        path = joinpath(PREALIGNED_DIR, string(name, ext))
    elseif is_aligned(index)
        path = joinpath(ALIGNED_DIR, string(name, ext))
    elseif is_finished(index)
        path = joinpath(FINISHED_DIR, string(name, ext))
    else
        if cur_dataset == "zebrafish"
            section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
            path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ext))
        elseif cur_dataset == "piriform"
            section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
            path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ext))
        else
            path = joinpath(BUCKET, RAW_DIR, string(name, ext))
        end
    end
    # println(path)
    return path
end

function get_path(name::String)
    return get_path(parse_name(name))
end

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

function parse_index(s::String)
    m = Base.match(r"(\d*),(\d*),(\-\d+|\d+),(\-\d+|\d+)", s)
    return (parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4]))
end

function parse_name(name::String)

    ret = (0, 0, 0, 0)
    # singleton tile
    m = Base.match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2])
    end

    # overview image
    m = Base.match(r"MontageOverviewImage_S2-W00(\d*)_sec(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), OVERVIEW_INDEX, OVERVIEW_INDEX   
    end

    # montaged section
    m = Base.match(r"(\d*),(\d*)_montaged", name)
    if typeof(m) != Void
    ret = montaged(parse(Int, m[1]), parse(Int, m[2]))
    end

    # prealigned section
    m = Base.match(r"(\d*),(\d*)_prealigned", name)
    if typeof(m) != Void
    ret = prealigned(parse(Int, m[1]), parse(Int, m[2]))
    end

    # prealigned_subsection
    m = match(r"(\d*),(\d*)_prealigned_(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), PREALIGNED_INDEX, parse(Int, m[3])
    end

    # aligned-section
    m = Base.match(r"(\d*),(\d*)_aligned", name)
    if typeof(m) != Void
    ret = aligned(parse(Int, m[1]), parse(Int, m[2]))
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

if contains(gethostname(), "seunglab") || contains(gethostname(), "seungom") || contains(gethostname(), "spock") 
 bucket_dir_path = "/mnt/bucket/labs/seung/"
end

if contains(gethostname(), "seungworkstation")
 bucket_dir_path = joinpath(homedir(), "seungmount/")
end

tracerhostnames = ["seungworkstation04", "mycroft", "seungworkstation10", "hogwarts"]
if gethostname() in tracerhostnames
 bucket_dir_path = joinpath(homedir(), "seungmount/Omni/alignment/datasets/")
end

if isdefined(:omni)
    if omni
        bucket_dir_path = joinpath(homedir(), "seungmount/Omni/alignment/datasets/")
    end
end  

#=if isfile("bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("bucket_dir_path.txt"), '\n')
elseif isfile("../bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("../bucket_dir_path.txt"), '\n')
else
    bucket_dir_path = joinpath(homedir(), "seungmount/Omni/alignment/datasets")
    if isdefined(:training)
        if training
            println("TRAINING PATHS LOADED")
            bucket_dir_path = joinpath(homedir(), "seungmount/Omni/alignment/training")
        end
    end
    if isdefined(:seunglabs)
        if seunglabs
            println("SEUNGLABS PATHS LOADED")
            bucket_dir_path = "/mnt/bucket/labs/seung/Omni/alignment/datasets"
        end
    end
end
=#

datasets_dir_path = "research/Julimaps/datasets"
# cur_dataset = "piriform"
# cur_dataset = "AIBS_practice_234251S6R_01_01_aligned_01"
cur_dataset = "AIBS_actual_trial10_full"
# cur_dataset = "align_net"
# cur_dataset = "AIBS_pilot_v1"
# cur_dataset = "align_net"
# cur_dataset = "elastic_test_crack"
# cur_dataset = "elastic_real_crack_with_cropping"
# cur_dataset = "elastic_real_compression"
# cur_dataset = "elastic_no_crack"
# cur_dataset = "elastic_compression_automated"
in_alignment_test = false
# test_dataset = "AIBS_practice_spring_constants"
# test_dataset = "AIBS_practice_broken_springs"
# test_dataset = "AIBS_practice_broken_springs_no_bug"
# test_dataset = "AIBS_practice_broken_springs_spring_constants"
# test_dataset = "AIBS_practice_broken_springs_fixed_springs"
# test_dataset = "AIBS_practice_broken_section"
# test_dataset = "AIBS_practice_broken_section_no_bug"
# test_dataset = "elastic_test_crack_removed_matches"
# test_dataset = "elastic_test_crack_broken_springs"
# test_dataset = "elastic_compression_automated_debug"
# test_dataset = "elastic_test_crack_removed_matches_broken_springs"
# test_dataset = "elastic_test_crack_removed_matches_broken_springs_zeros"
# test_dataset = "elastic_test_crack_removed_matches_broken_springs_finer_mesh"
#cur_dataset = "zebrafish"
affine_dir_path = "~"

raw_dir_path = "0_raw"
premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"
finished_dir_path = "5_finished"

expunged_dir = "expunged"
stack_dir = "stacks"

wafer_filename = "wafer_paths.txt"
premontaged_registry_filename = "registry_premontaged.txt"
montaged_registry_filename = "registry_montaged.txt"
prealigned_registry_filename = "registry_prealigned.txt"
aligned_registry_filename = "registry_aligned.txt"
expunged_registry_filename = "expunged_aligned.txt"

function check_dataset_dir(dataset_name)
    
    function setup_dir(dir)
        if !isdir(dir)
            println("Creating $dir")
            mkdir(dir)
        end
    end

    dataset_dir = joinpath(bucket_dir_path, datasets_dir_path, dataset_name)
    dirs = [raw_dir_path, premontaged_dir_path, montaged_dir_path, 
                prealigned_dir_path, aligned_dir_path, finished_dir_path]
    setup_dir(dataset_dir)
    for d in dirs
        path = joinpath(dataset_dir, d)
        setup_dir(path)
        review_path = joinpath(path, "review")
        setup_dir(review_path)
        mask_path = joinpath(path, "mask")
        setup_dir(mask_path)
    end
end

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, PREMONTAGED_OFFSETS, PREMONTAGE_DIR, ALIGNMENT_DIR, INSPECTION_DIR

global BUCKET = bucket_dir_path
global AFFINE_DIR = affine_dir_path
global DATASET_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset)
global RAW_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, raw_dir_path)
global PREMONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path)
global MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path)
global PREALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path)
if in_alignment_test
    check_dataset_dir(test_dataset)
    global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, test_dataset, aligned_dir_path)
    global FINISHED_DIR = joinpath(bucket_dir_path, datasets_dir_path, test_dataset, finished_dir_path)
    global STACKS_DIR = joinpath(bucket_dir_path, datasets_dir_path, test_dataset, finished_dir_path, stack_dir)
    aligned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, test_dataset, aligned_dir_path, aligned_registry_filename)
    global REGISTRY_ALIGNED = parse_registry(aligned_registry_path)
else
    check_dataset_dir(cur_dataset)
    global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path)
    global FINISHED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, finished_dir_path)
    global STACKS_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, finished_dir_path, stack_dir)
    aligned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path, aligned_registry_filename)
    global REGISTRY_ALIGNED = parse_registry(aligned_registry_path)
end

global EXPUNGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, raw_dir_path, expunged_dir)

waferpath_filename = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, wafer_filename)
global WAFER_DIR_DICT = waferpaths_to_dict(waferpath_filename)

premontaged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path, premontaged_registry_filename)
global REGISTRY_PREMONTAGED = parse_registry(premontaged_registry_path)
#REGISTRY_PREMONTAGED = hcat(REGISTRY_PREMONTAGED, fill(8000, size(REGISTRY_PREMONTAGED)[1], 2));

montaged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path, montaged_registry_filename)
global REGISTRY_MONTAGED = parse_registry(montaged_registry_path)

prealigned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path, prealigned_registry_filename)
global REGISTRY_PREALIGNED = parse_registry(prealigned_registry_path)

expunged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, raw_dir_path, expunged_dir, expunged_registry_filename)
global REGISTRY_EXPUNGED = parse_registry(expunged_registry_path)

show_plot = false

if cur_dataset == "piriform"
    global ROI_FIRST = (1,2,0,0);
    global ROI_LAST = (8,173,0,0);
    global DATASET_RESOLUTION = [7,7,40]
else
    global ROI_FIRST = (1,1,0,0);
    global ROI_LAST = (1,102,0,0);
    global DATASET_RESOLUTION = [4,4,40]
end

if cur_dataset == "piriform" || cur_dataset == "zebrafish"
    global TILE_SIZE = [8000, 8000]
    global EXPECTED_TILE_OVERLAP = 0.145
else
    global TILE_SIZE = [3840, 3840]
    global EXPECTED_TILE_OVERLAP = 0.18
end