
function is_overview(index::Index)
	if index[3:4] == (OVERVIEW_INDEX, OVERVIEW_INDEX)	return true; else return false; end
end

function is_montaged(index::Index)
	if index[3:4] == (MONTAGED_INDEX, MONTAGED_INDEX)	return true; else return false; end
end

function is_pre_aligned(index::Index)
	if index[3:4] == (PRE_ALIGNED_INDEX, PRE_ALIGNED_INDEX)	return true; else return false; end
end

function is_aligned(index::Index)
	if index[3:4] == (ALIGNED_INDEX, ALIGNED_INDEX)	return true; else return false; 
end
end


function parseName(name::String)

	ret = (0, 0, 0, 0);
	# singleton tile
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", name)
	if typeof(m) != Void
	ret = parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2])
	end

	# overview image
	m = match(r"MontageOverviewImage_S2-W00(\d*)_sec(\d*)", name);
	if typeof(m) != Void
	ret = parse(Int, m[1]), parse(Int, m[2]), OVERVIEW_INDEX, OVERVIEW_INDEX; 	
	end

	# montaged section
	m = match(r"(\d*),(\d*)_montaged", name)
	if typeof(m) != Void
	ret = parse(Int, m[1]), parse(Int, m[2]), MONTAGED_INDEX, MONTAGED_INDEX; 	
	end

	# pre-aligned section
	m = match(r"(\d*),(\d*)_pre-aligned", name)
	if typeof(m) != Void
	ret = parse(Int, m[1]), parse(Int, m[2]), PRE_ALIGNED_INDEX, PRE_ALIGNED_INDEX; 
	end

	# aligned-section
	m = match(r"(\d*),(\d*)_aligned", name)
	if typeof(m) != Void
	ret = parse(Int, m[1]), parse(Int, m[2]), ALIGNED_INDEX, ALIGNED_INDEX; 
	end

	return ret;
	
end

function getName(index::Index)
	if is_overview(index)
	return string("MontageOverviewImage_S2-W00", index[1], "_sec", index[2]);
	elseif is_montaged(index)
	return string(index[1], ",", index[2], "_montaged");
	elseif is_pre_aligned(index)
	return string(index[1], ",", index[2], "_pre-aligned");
	elseif is_aligned(index)
	return string(index[1], ",", index[2], "_aligned")
	else
	return string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2]);
	end
end

# function getPath()
# methods: 
#	  
# extensions:
# Mesh.jl: getPath(mesh::Mesh)
function getPath(index::Index)
	name = getName(index);

	if is_overview(index)
	section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage");
	return joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"));

	elseif is_montaged(index)
	return joinpath(MONTAGED_DIR, string(name, ".tif"));
	
	elseif is_pre_aligned(index)
	return joinpath(PRE_ALIGNED_DIR, string(name, ".tif"));

	elseif is_aligned(index)
	return joinpath(ALIGNED_DIR, string(name, ".tif"));

	else
	section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage");
	return joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"));
	end
end

function getPath(name::String)
	return getPath(parseName(name));
end



# extensions:
# Mesh.jl getImage(mesh::Mesh)
function getImage(path::String)
	img = imread(path);
	img = img[:, :, 1];
	img.properties["timedim"] = 0;
	return convert(Array{UInt8, 2}, round(convert(Array, img)*255));
end

function getImage(index::Index)
	return getImage(getPath(index));
end


function waferpaths2dict(wafer_path_file)
	warray = readdlm(wafer_path_file)
	wdict = Dict()
	for (idx, path) in zip(warray[:,1], warray[:,2])
		wdict[idx] = path
	end
	return wdict
end

function parse_offsets(info_path::String)
	file = readdlm(info_path);
	session = cell(size(file, 1), 4); # name, index, dx, dy
	for i in 1:size(file, 1)
	index = parseName(file[i, 1]);
	session[i, 1] = getName(index);
	session[i, 2] = index;
	session[i, 3] = file[i, 2];
	session[i, 4] = file[i, 3];
	end
	return session;
end

bucket_dir_path = "/usr/people/dih/seungmount";
datasets_dir_path = "research/Julimaps/datasets";
cur_dataset = "piriform";
affine_dir_path = "~";

pre_montaged_dir_path = "1_pre-montaged";
montaged_dir_path = "2_montaged";
pre_aligned_dir_path = "3_pre-aligned";
aligned_dir_path = "4_aligned";
wafer_filename = "wafer_paths.txt";
pre_montaged_offsets_filename = "pre-montaged_offsets.txt";
pre_aligned_offsets_filename = "pre-aligned_offsets.txt";

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, PRE_MONTAGED_OFFSETS, PRE_MONTAGE_DIR, ALIGNMENT_DIR

global BUCKET = bucket_dir_path;
global AFFINE_DIR = affine_dir_path;
global DATASET_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset);
global PRE_MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, pre_montaged_dir_path);
global MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path);
global PRE_ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, pre_aligned_dir_path);
global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path);
global WAFER_DIR_DICT = waferpaths2dict(joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, wafer_filename));
global PRE_MONTAGED_OFFSETS = parse_offsets(joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, pre_montaged_dir_path, pre_montaged_offsets_filename));
global PRE_ALIGNED_OFFSETS = parse_offsets(joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, pre_aligned_dir_path, pre_aligned_offsets_filename));



export tile_size, block_size, search_r, min_r, mesh_length, mesh_coeff, match_coeff, eta_grad, eta_newton, show_plot, num_procs, ftol_grad, ftol_newton, num_tiles, num_rows, num_cols, mesh_length_alignment, min_r_alignment, search_r_alignment, block_size_alignment;

tile_size = 8000;
block_size = 40;
search_r = 80;
min_r = 0.75;
mesh_length = 200;
block_size_alignment = 150;
search_r_alignment = 500;
min_r_alignment = 0.15;
mesh_length_alignment = 1500;
mesh_coeff = 1;
match_coeff = 20;
eta_grad = 0.01;
eta_newton = .5; 
show_plot = false;
num_procs = nprocs();
ftol_grad = 1/500;
ftol_newton = 1/1000000;
num_tiles = 16;
num_rows = 4;
num_cols = 4;




