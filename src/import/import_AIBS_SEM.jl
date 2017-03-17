using LightXML

function init_sem_load_dict()
    sem_load_fn = joinpath(BUCKET, DATASET, "SEM_z2sess_again.csv")
    sem_list = readdlm(sem_load_fn,',')
    sem = Dict()
    for i in 1:size(sem_list,1)
        sem[sem_list[i,1]] = sem_list[i,2]
    end
    return sem
end

function load_sem_meta_files(z_index)
    sem = init_sem_load_dict()
    src_fn = joinpath(GCLOUD_RAW_DIR, sem[z_index], "MosaicInfo*.ve-updates")
    dst_fn = get_sem_meta_filepath(z_index)
    cmd = `sudo gsutil -m cp $src_fn $dst_fn`
    # println(cmd)
    Base.run(cmd)
end

function import_sem_meta_files()
    z_indices = [4019, 4337, 3624, 4485, 3533, 3814, 3907]
    for z_index in z_indices
        load_sem_meta_files(z_index)
    end
end

function get_sem_meta_filepath(z_index)
    return joinpath(BUCKET, DATASET, "$(z_index).xml")
end

function load_sem_meta_filepath(z_index)
    sem_filepath = get_sem_meta_filepath(z_index)
    xml = parse_file(sem_filepath)
    xroot = LightXML.root(xml)



end
