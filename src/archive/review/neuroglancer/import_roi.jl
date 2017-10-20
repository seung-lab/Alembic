"""
Convert & save all CSV files in dir to mask polygon points files
"""
function import_neuroglancer_rois(dir)
    filenames = [fn for fn in readdir(dir) if fn[end-2:end] == "csv"]
    for fn in filenames
        println("Converting ", fn)
        convert_neuroglancer_csv(dir, fn)
    end
end

"""
Take neuroglancer roi_select csv & convert to image coords for mask

Inputs:
    dir: directory with the neuroglancer CSV files
    fn: filename of the particular section
    start_file_z: starting z index of neuroglancer (sometimes not 1)

Outpus:
    None; writes converted polygon points file to masks dir
"""
function convert_neuroglancer_csv(dir, fn, start_file_z=16384)
    i = Int(parse(fn[1:end-3])) - start_file_z + 1
    index = REGISTRY_ALIGNED[i,2]
    pts = readdlm(joinpath(dir,fn), ',', Int64)[:,1:2]
    new_pts = copy(pts)
    offset = round(Int64, get_offset(index))
    img_size = get_image_size(index)
    for i in 1:size(pts,1)
        new_pts[i,1] = min(max(pts[i,1]-16384-offset[1], 1), img_size[1]-1)
        new_pts[i,2] = min(max(pts[i,2]-16384-offset[2], 1), img_size[2]-1)
    end
    new_fn = get_path("mask", index)
    writedlm(new_fn, new_pts)
end

