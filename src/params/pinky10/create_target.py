from cloudvolume import CloudVolume
from cloudvolume.lib import Vec, Bbox
import numpy as np
from copy import deepcopy

padding = Vec(128, 128, 0)
start = Vec(20002, 12009, 65) - padding
stop = Vec(22050, 14057, 165) + padding
bbox = Bbox(start, stop)

src = CloudVolume('gs://neuroglancer/pinky10/image_rechunked', mip=0, fill_missing=False)
info = deepcopy(src.info)
scale = None
for s in info['scales']:
	if s['key'] == '4_4_40':
		scale = s
		break

scale['chunk_sizes'] = [list(map(int, [bbox.size3()[0], bbox.size3()[1], 1]))]
scale['size'] = list(map(int, bbox.size3()))
scale['voxel_offset'] = list(map(int, start))

info['scales'] = [scale]
dst = CloudVolume('gs://neuroglancer/pinky10/ground_truth_registration/target', mip=0,
									info=info, cdn_cache=False)													
dst.commit_info()

img = src[bbox.to_slices()]
dst[:,:,:,:] = img


# create destination for warped image & segmentations
dirs = {'img': 'image',
		'mit_new': 'mit_new',
		'seg_original': 'seg_original',
		'seg_invagination_new': 'seg_invagination_new',
		'psd': 'psd'}
for old_dir, new_dir in dirs.items():
	src = CloudVolume('gs://neuroglancer/kisuk/pinky/ground_truth/stitched_vol19-vol34/{0}'.format(old_dir), mip=0, fill_missing=True, bounded=False)
	info = deepcopy(src.info)
	padding = Vec(256, 256, 0)
	start = Vec(20002, 12009, 65) - padding
	stop = Vec(22050, 14057, 165) + padding
	bbox = Bbox(start, stop)
	scale['chunk_sizes'] = [list(map(int, [bbox.size3()[0], bbox.size3()[1], 1]))]
	scale['size'] = list(map(int, bbox.size3()))
	scale['voxel_offset'] = list(map(int, start))
	info['scales'] = [scale]
	new_cv = CloudVolume('gs://neuroglancer/pinky10/ground_truth_registration/dst/{0}'.format(new_dir), mip=0,
										info=info, cdn_cache=False)	
	new_cv.commit_info()