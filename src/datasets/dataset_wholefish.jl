global BUCKET = "/home/ubuntu"
#global DATASET = "datasets/cremi"
global DATASET = "datasets/wholefish"
#global ROI_FIRST = (1,1,0,0);
global ROI_FIRST = (1,1,0,0);
global ROI_LAST = (6,200,0,0);
global DATASET_RESOLUTION = [4,4,40]

global TASKS_LOCALE = "gcs"
#global TASKS_LOCALE = "aws"

if TASKS_LOCALE == "aws"
global TASKS_TASK_QUEUE_NAME = "task-queue-TEST";
global TASKS_ERROR_QUEUE_NAME = "error-queue-TEST";
global TASKS_DONE_QUEUE_NAME = "done-queue-TEST";
global TASKS_REGISTRY_QUEUE_NAME = "registry-queue-TEST";
global TASKS_BUCKET_NAME = "seunglab";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;
end

if TASKS_LOCALE == "gcs"
global TASKS_TASK_QUEUE_NAME = "task-queue-GCS";
global TASKS_ERROR_QUEUE_NAME = "error-queue-GCS";
global TASKS_DONE_QUEUE_NAME = "done-queue-GCS";
global TASKS_REGISTRY_QUEUE_NAME = "registry-queue-GCS";
global TASKS_BUCKET_NAME = "image_assembly";
global TASKS_CACHE_DIRECTORY = BUCKET;
global TASKS_BASE_DIRECTORY = DATASET;
global TASKS_POLL_FREQUENCY = 10;
end

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
		if index[2] < 10
    return string("Tile_r", index[3], "-c", index[4], "_W00", index[1], "_sec00", index[2])
    		elseif index[2] < 100
    return string("Tile_r", index[3], "-c", index[4], "_W00", index[1], "_sec0", index[2])
    		else
    return string("Tile_r", index[3], "-c", index[4], "_W00", index[1], "_sec", index[2])
  end
  	else
		if index[2] < 10
    return string("Tile_r", index[3], "-c", index[4], "_W0", index[1], "_sec00", index[2])
    		elseif index[2] < 100
    return string("Tile_r", index[3], "-c", index[4], "_W0", index[1], "_sec0", index[2])
    		else
    return string("Tile_r", index[3], "-c", index[4], "_W0", index[1], "_sec", index[2])
  end
  end
    end
end


function import_flat(index)
          path = joinpath(PREMONTAGED_DIR_PATH, string(get_name_legacy(index), ".tif"))
	          if isfile(path)
		                    img = Images.load(path)
				                    save(get_path(index), raw(img[:,:,1]'))
						                    update_registry(index, image_size = [8000,8000], offset = [7200 * index[3], 7200 * index[4]])
								            end
									  end
