function is_overview(index::Index)
    if index[3:4] == (OVERVIEW_INDEX, OVERVIEW_INDEX)   return true else return false end
end

function is_premontaged(index::Index)
    if index[3] > 0 && index[4] > 0 return true else return false end
end

function is_montaged(index::Index)
    if index[3:4] == (MONTAGED_INDEX, MONTAGED_INDEX)   return true else return false end
end

function is_prealigned(index::Index)
    if index[3:4] == (PREALIGNED_INDEX, PREALIGNED_INDEX)   return true else return false end
end

function is_aligned(index::Index)
    if index[3:4] == (ALIGNED_INDEX, ALIGNED_INDEX) return true else return false end
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
    return string(index[1], ",", index[2], "_prealigned")
    elseif is_aligned(index)
    return string(index[1], ",", index[2], "_aligned")
    else
    return string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
    end
end

# function get_path()
# methods: 
#     
# extensions:
# Mesh.jl: get_path(mesh::Mesh)
function get_path(index::Index)
    name = get_name(index)
    if is_overview(index)
        if cur_dataset == "zebrafish"
            section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
        else
            section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        end
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"))
    elseif is_montaged(index)
        path = joinpath(MONTAGED_DIR, string(name, ".tif"))
    elseif is_prealigned(index)
        path = joinpath(PREALIGNED_DIR, string(name, ".tif"))
    elseif is_aligned(index)
        path = joinpath(ALIGNED_DIR, string(name, ".tif"))
    else
        if cur_dataset == "zebrafish"
            section_folder = string("W00", index[1], "_Sec", index[2], "_Montage")
        else
            section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        end
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"))
    end
    println(path)
    return path
end

function get_path(name::String)
    return get_path(parse_name(name))
end



# extensions:
# Mesh.jl get_image(mesh::Mesh)
function get_image(path::String)
    img = imread(path)
    img = img[:, :, 1]
    img.properties["timedim"] = 0
    return convert(Array{UInt8, 2}, round(convert(Array, img)*255))
end

function get_image(index::Index)
    return get_image(get_path(index))
end

function get_thumbnail_path(index::Index)
  fn = string(join(index[1:2], ","), "_thumb.png")
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

function waferpaths2dict(waferpath_filename)
    wdict = Dict()
    if isfile(waferpath_filename)
        warray = readdlm(waferpath_filename)
        for (idx, path) in zip(warray[:,1], warray[:,2])
            wdict[idx] = path
        end
    end
    return wdict
end

function parse_offsets(path::String)
    offsets = cell(0, 0)
    if isfile(path)
        file = readdlm(path)
        offsets = cell(size(file, 1), size(file, 2) + 1) # name, index, dx, dy
        for i in 1:size(offsets, 1)
            index = parse_name(file[i, 1])
            offsets[i, 1] = get_name(index)
            offsets[i, 2] = index
            for j in 3:size(offsets, 2)
                offsets[i, j] = file[i, j-1]
            end
        end
      end
    return offsets
end

bucket_dir_path = ""
if isfile("bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("bucket_dir_path.txt"), '\n')
elseif isfile("../bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("../bucket_dir_path.txt"), '\n')
end
datasets_dir_path = "research/Julimaps/datasets"
cur_dataset = "piriform"
#cur_dataset = "zebrafish"
affine_dir_path = "~"

premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"

wafer_filename = "wafer_paths.txt"
premontaged_offsets_filename = "premontaged_offsets.txt"
montaged_offsets_filename = "montaged_offsets.txt"
prealigned_offsets_filename = "prealigned_offsets.txt"
aligned_offsets_filename = "aligned_offsets.txt"

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, PREMONTAGED_OFFSETS, PREMONTAGE_DIR, ALIGNMENT_DIR

global BUCKET = bucket_dir_path
global AFFINE_DIR = affine_dir_path
global DATASET_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset)
global PREMONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path)
global MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path)
global PREALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path)
global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path)

waferpath_filename = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, wafer_filename)
global WAFER_DIR_DICT = waferpaths2dict(waferpath_filename)

premontaged_offsets_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path, premontaged_offsets_filename)
global PREMONTAGED_OFFSETS = parse_offsets(premontaged_offsets_path)

montaged_offsets_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path, montaged_offsets_filename)
global MONTAGED_OFFSETS = parse_offsets(montaged_offsets_path)

prealigned_offsets_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path, prealigned_offsets_filename)
global PREALIGNED_OFFSETS = parse_offsets(prealigned_offsets_path)

aligned_offsets_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path, aligned_offsets_filename)
global ALIGNED_OFFSETS = parse_offsets(aligned_offsets_path)

global GLOBAL_BB = BoundingBox(0,0,36000,36000)

show_plot = false
num_procs = nprocs()

