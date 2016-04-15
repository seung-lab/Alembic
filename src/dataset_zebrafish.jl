global ROI_FIRST = (2,31,0,0);
global ROI_LAST = (9,163,0,0);

function get_name(index)
    if is_overview(index)
        return string("MontageOverviewImage_W00", index[1], "_sec", index[2])
    elseif is_montaged(index)
        return string(index[1], ",", index[2], "_montaged")
    elseif is_prealigned(index)
        return string(index[1], ",", index[2], "_prealigned")
    elseif is_aligned(index)
        return string(index[1], ",", index[2], "_aligned")
    elseif is_finished(index)
        return string(index[1], ",", index[2], "_finished")
    else
    return string("Tile_r", index[3], "-c", index[4], "_W00", index[1], "_sec", index[2])
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
        section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
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
#        section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        #path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ext))
        path = joinpath(BUCKET, PREMONTAGED_DIR, string(name, ext))
    end
    return path
end

function get_path(name::String)
    return get_path(parse_name(name))
end


function get_image(index::Index)
    return get_image(get_path(index))
end

function get_thumbnail_path(index::Index)
  fn = string(join(index[1:2], ","), "_thumb.jpg")
  println(fn)
  return joinpath(MONTAGED_DIR, "review", fn) 
end

function get_thumbnail_path(indexA::Index, indexB::Index)
  fn = string(join(indexA[1:2],","), "-", join(indexB[1:2],","), "_thumb.png")
  println(fn)
  if is_prealigned(indexA)
    return joinpath(PREALIGNED_DIR, "review", fn) 
  else
    return joinpath(ALIGNED_DIR, "review", fn) 
  end
end

function parse_index(s::String)
    m = match(r"(\d*),(\d*),(\-\d+|\d+),(\-\d+|\d+)", s)
    return (parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4]))
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

bucket_dir_path = ""

if contains(gethostname(), "seunglab") || contains(gethostname(), "spock") 
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

datasets_dir_path = joinpath(homedir(), "datasets")
#datasets_dir_path = joinpath(homedir(), "seungmount", "research", "Julimaps", "datasets")
cur_dataset = "zebrafish"
affine_dir_path = "~"

premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"
finished_dir_path = "5_finished"

premontaged_registry_filename = "registry_premontaged.txt"
montaged_registry_filename = "registry_montaged.txt"
prealigned_registry_filename = "registry_prealigned.txt"
aligned_registry_filename = "registry_aligned.txt"

export BUCKET, DATASET_DIR, AFFINE_DIR, PREMONTAGED_OFFSETS, PREMONTAGE_DIR, ALIGNMENT_DIR, INSPECTION_DIR

global BUCKET = bucket_dir_path
global AFFINE_DIR = affine_dir_path
global DATASET_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset)
global PREMONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path)
global MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path)
global PREALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path)
global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path)
global FINISHED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, finished_dir_path)

premontaged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path, premontaged_registry_filename)
global REGISTRY_PREMONTAGED = parse_registry(premontaged_registry_path)
#REGISTRY_PREMONTAGED = hcat(REGISTRY_PREMONTAGED, fill(8000, size(REGISTRY_PREMONTAGED)[1], 2));

montaged_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path, montaged_registry_filename)
global REGISTRY_MONTAGED = parse_registry(montaged_registry_path)

prealigned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path, prealigned_registry_filename)
global REGISTRY_PREALIGNED = parse_registry(prealigned_registry_path)

aligned_registry_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path, aligned_registry_filename)
global REGISTRY_ALIGNED = parse_registry(aligned_registry_path)


