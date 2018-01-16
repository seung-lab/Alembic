from cloudvolume import CloudVolume
import numpy as np
from scipy import ndimage

src = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/cfsplit_manual', mip=5, cdn_cache=False);
roi = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/roicc', mip=6);
dst = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/cfsplit_manual_cc', mip=5, cdn_cache=False);
src_val = 0
roi_val = 1
# cc_filter = np.ones((3,3))

def get_whole_slice(cv):
	offset = cv.voxel_offset
	img_size = cv.shape
	x_slice = slice(offset[0], offset[0]+img_size[0])
	y_slice = slice(offset[1], offset[1]+img_size[1])
	return x_slice, y_slice

def get_image(cv, z):
	s = get_whole_slice(cv) + (z, )
	return cv[s][:,:,0,0]

for z in l: #range(1,2176):
	print(z)
	src_img = get_image(src, z)
	roi_img = get_image(roi, z).repeat(2, axis=0).repeat(2, axis=1)
	mask = np.logical_and((src_img == src_val), (roi_img == roi_val))
	labeled, nr_objects = ndimage.label(mask) #, cc_filter)
	for k in range(nr_objects):
		cc = labeled == k+1
		# print((k, np.sum(cc)))
		if np.sum(cc) < 625:
			labeled[cc] = 0
	labeled, nr_objects = ndimage.label(labeled, cc_filter)
	labeled = labeled.astype(np.uint8)
	dst[get_whole_slice(dst) + (z,)] = np.reshape(labeled, labeled.shape+(1,))

