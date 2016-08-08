# EXAMPLE ENTRY FOR THINGS THAT MUST BE SPECIFIED IN dataset_*.jl
# global BUCKET = "/home/ubuntu/datasets" 	# if BUCKET is different for each computer, then use some if statements
#
# global DATASET = "zebrafish	"		# each DATASET lives under BUCKET, i.e. BUCKET/DATASET
#						# i.e. if there are two datasets named "zebrafish" and "piriform"
#						# zebrafish folders will exist as /home/ubuntu/datasets/zebrafish/0_overview, /home/ubuntu/datasets/zebrafish/1_premontaged...
#						# piriform folders will exist as /home/ubuntu/datasets/piriform/0_overview, /home/ubuntu/datasets/piriform/1_premontaged...
#
# global DATASET_RESOLUTION = [7,7,40] 	     	# dataset resolution in each dimension (i, j, k) or (y, x, z)
#
# global ROI_FIRST = (2,3,0,0)			# where the dataset starts - used for prealignment
# global ROI_LAST = (9,163,0,0)			# where the dataset ends

# gets the canonical name of the image associated with the index
function get_name(index::Index)
    if is_overview(index)	    	return string(index[1], ",", index[2], "_overview")
    elseif is_premontaged(index)	return string(index)
    elseif is_montaged(index)		return string(index[1], ",", index[2], "_montaged")
    elseif is_prealigned(index)
      if is_subsection(index)		return string(index[1], ",", index[2], "_prealigned_", index[4])
      else 				return string(index[1], ",", index[2], "_prealigned") end
    elseif is_aligned(index)		return string(index[1], ",", index[2], "_aligned")
    elseif is_finished(index) 		return string(index[1], ",", index[2], "_finished")
    end
end
# used for loading - i.e. getting the name of the Mesh by calling get_name((2,3,1,4), Mesh) returns "Mesh((2,3,1,4))"
function get_name(object_type::Union{DataType, String}, index)	
	if object_type == "MeshSet" || string(object_type) == "MeshSet"
	  # hack for old names
	  if typeof(index) == Tuple{Index, Index}
		firstindex = index[1]
		lastindex = index[2]
		if is_premontaged(firstindex) return "$(firstindex[1]),$(firstindex[2])_montaged";
		elseif is_montaged(firstindex) return "$(firstindex[1]),$(firstindex[2])_prealigned";
		else  return "$(firstindex[1]),$(firstindex[2])-$(lastindex[1]),$(lastindex[2])_aligned"; end
	      else
		if is_montaged(index) return "$(index[1]),$(index[2])_montaged";
		elseif is_prealigned(index) return "$(index[1]),$(index[2])_prealigned"; end
	      end
	      end
	      #end hack
  	# singleton case
	if typeof(index) == Index 			return string(object_type, "(", index, ")")	
        elseif typeof(index) == Tuple{Index, Index} 	return string(object_type, index)	end
end
# used for saving - gets the canonical name of the object
function get_name(object)		
	  # hack for old names
	#if object == "MeshSet" || string(typeof(object)) == "MeshSet"
	if string(typeof(object)) == "MeshSet"
		firstindex = get_index(object.meshes[1]);
		lastindex = get_index(object.meshes[end]);
		if is_premontaged(firstindex) return "$(firstindex[1]),$(firstindex[2])_montaged";
		elseif is_montaged(firstindex) return "$(firstindex[1]),$(firstindex[2])_prealigned";
		else  return "$(firstindex[1]),$(firstindex[2])-$(lastindex[1]),$(lastindex[2])_aligned"; end
	else
	      #end hack
	index = get_index(object)
	if typeof(index) == Index 			return string(typeof(object), "(", index, ")")	
        elseif typeof(index) == Tuple{Index, Index} 	return string(typeof(object), index)	end
      end
end

# gets the full directory of the index as if it were an image - e.g. .../1_premontaged
function get_dir_path(index::Union{Index, Tuple{Index, Index}})
  if typeof(index) != Index index = index[1] end
    if is_overview(index)		return OVERVIEW_DIR_PATH
    elseif is_premontaged(index) 	return PREMONTAGED_DIR_PATH
    elseif is_montaged(index) 		return MONTAGED_DIR_PATH
    elseif is_prealigned(index) 	return PREALIGNED_DIR_PATH
    elseif is_aligned(index) 		return ALIGNED_DIR_PATH
    elseif is_finished(index) 		return FINISHED_DIR_PATH	end
end

# gets the full directory of the object based on get_index
function get_dir_path(object)
  index = get_index(object); if typeof(index) != Index index = index[1]	end
  return get_dir_path(index)
end

# gets the full directory of the object - e.g. .../1_premontaged
function get_subdir(object)  			return get_subdir(typeof(object))	end
function get_subdir(object_type::DataType)	return get_subdir(string(object_type));	end

function get_subdir(string::String)
  if 	 string == "Mesh"		return MESH_DIR, ".jls"
  elseif string == "Match"     		return MATCH_DIR, ".jls"
  elseif string == "MeshSet"    	return MESHSET_DIR, ".jls"
  elseif string == "review"     	return REVIEW_DIR, ".h5"
  elseif string == "import"     	return IMPORT_DIR, ".txt"
  elseif string == "contrast_bias"      return CONTRAST_BIAS_DIR, ".h5"
  elseif string == "contrast_auto"      return CONTRAST_AUTO_DIR, ".txt"
  elseif string == "stats"     		return STATS_DIR, ".txt"
  elseif string == "mask"     		return MASK_DIR, ".png"
  elseif string == "outline"     	return OUTLINE_DIR, ".png"
  elseif string == "expunge"     	return EXPUNGED_DIR, ".h5"
  elseif string == "thumbnail"     	return THUMBNAIL_DIR, ".h5"
  end
end

# 

# function get_path()
# methods: 
#     
# extensions:
# Mesh.jl: get_path(mesh::Mesh)
function get_path(index::Index, ext = ".h5")
  return joinpath(get_dir_path(index), string(get_name(index), ext))
end
function get_path(object_type::Union{DataType, String}, index)
  # hack to support singleton load for meshsets
  if (object_type == "MeshSet" || string(object_type) == "MeshSet") && typeof(index) == Index
  return joinpath(get_dir_path(prevstage(index)), get_subdir(object_type)[1], string(get_name(object_type, index), get_subdir(object_type)[2]))
  end
  return joinpath(get_dir_path(index), get_subdir(object_type)[1], string(get_name(object_type, index), get_subdir(object_type)[2]))
end
function get_path(object)
  index = get_index(object); if typeof(index) != Index index = index[1]	end
  return joinpath(get_dir_path(index), get_subdir(object)[1], string(get_name(object), get_subdir(object)[2]))
end
function get_path(name::String)
    return get_path(parse_name(name))
end

function parse_index(s::String)
    m = Base.match(r"(\d+),(\d+),(\-\d+|\d+),(\-\d+|\d+)", s)
    return (parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4]))
end

function parse_name(name::String)

    ret = (0, 0, 0, 0)
    # singleton tile
    m = Base.match(r"(\d+),(\d+),(\d+),(\d+)", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4])
    end

    # overview image
    m = Base.match(r"(\d+),(\d+)_overview", name)
    if typeof(m) != Void
    ret = overview(parse(Int, m[1]), parse(Int, m[2]))
    end

    # montaged section
    m = Base.match(r"(\d+),(\d+)_montaged", name)
    if typeof(m) != Void
    ret = montaged(parse(Int, m[1]), parse(Int, m[2]))
    end

    # prealigned section
    m = Base.match(r"(\d+),(\d+)_prealigned", name)
    if typeof(m) != Void
    ret = prealigned(parse(Int, m[1]), parse(Int, m[2]))
    end

    # prealigned_subsection
    m = match(r"(\d+),(\d+)_prealigned_(\d+)", name)
    if typeof(m) != Void
    ret = subsection(parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]))
    end

    # aligned-section
    m = Base.match(r"(\d+),(\d+)_aligned", name)
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

global OVERVIEW_DIR = "0_overview"
global PREMONTAGED_DIR = "1_premontaged"
global MONTAGED_DIR = "2_montaged"
global PREALIGNED_DIR = "3_prealigned"
global ALIGNED_DIR = "4_aligned"
global FINISHED_DIR = "5_finished"

global MESH_DIR = "mesh"
global MATCH_DIR = "match"
global MESHSET_DIR = "meshset"
global EXPUNGED_DIR = "expunged"
global REVIEW_DIR = "review"
global MASK_DIR = "mask"
global IMPORT_DIR = "import"
global CONTRAST_BIAS_DIR = "bias"
global CONTRAST_AUTO_DIR = "auto"
global STATS_DIR = "stats"
global OUTLINE_DIR = "outline"
global THUMBNAIL_DIR = "thumbnail"

global OVERVIEW_DIR_PATH = joinpath(BUCKET, DATASET, OVERVIEW_DIR)
global PREMONTAGED_DIR_PATH = joinpath(BUCKET, DATASET, PREMONTAGED_DIR)
global MONTAGED_DIR_PATH = joinpath(BUCKET, DATASET, MONTAGED_DIR)
global PREALIGNED_DIR_PATH = joinpath(BUCKET, DATASET, PREALIGNED_DIR)
global ALIGNED_DIR_PATH = joinpath(BUCKET, DATASET, ALIGNED_DIR)
global FINISHED_DIR_PATH = joinpath(BUCKET, DATASET, FINISHED_DIR)

global REGISTRY_FILENAME = "registry.txt"
global REGISTRY_PREMONTAGED = parse_registry(joinpath(PREMONTAGED_DIR_PATH, REGISTRY_FILENAME))
global REGISTRY_MONTAGED = parse_registry(joinpath(MONTAGED_DIR_PATH, REGISTRY_FILENAME))
global REGISTRY_PREALIGNED = parse_registry(joinpath(PREALIGNED_DIR_PATH, REGISTRY_FILENAME))
global REGISTRY_ALIGNED = parse_registry(joinpath(ALIGNED_DIR_PATH, REGISTRY_FILENAME))

function get_registry_path(index)
    if is_premontaged(index) 		return joinpath(PREMONTAGED_DIR_PATH, REGISTRY_FILENAME)
    elseif is_montaged(index) 		return joinpath(MONTAGED_DIR_PATH, REGISTRY_FILENAME)
    elseif is_prealigned(index) 	return joinpath(PREALIGNED_DIR_PATH, REGISTRY_FILENAME)
    elseif is_aligned(index) 		return joinpath(ALIGNED_DIR_PATH, REGISTRY_FILENAME)
    end
end

function check_dirs(dataset_name::String = DATASET)
    function setup_dir(dir)
        if !isdir(dir)
            println("Creating $dir")
            mkdir(dir)
        end
    end

    dataset_dir = joinpath(BUCKET, dataset_name)
    dirs = [ OVERVIEW_DIR, PREMONTAGED_DIR, MONTAGED_DIR, PREALIGNED_DIR, ALIGNED_DIR, FINISHED_DIR ]
    subdirs = [ MESH_DIR, MATCH_DIR, MESHSET_DIR, EXPUNGED_DIR, REVIEW_DIR, MASK_DIR, IMPORT_DIR, CONTRAST_BIAS_DIR, CONTRAST_AUTO_DIR, STATS_DIR, OUTLINE_DIR, THUMBNAIL_DIR ]

    setup_dir(dataset_dir)
    for d in dirs
        path = joinpath(dataset_dir, d)
        setup_dir(path)
	for sd in subdirs
	  subpath = joinpath(path, sd)
       	  setup_dir(subpath)
	end
    end
end

check_dirs()
