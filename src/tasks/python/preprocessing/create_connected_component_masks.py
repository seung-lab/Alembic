from cloudvolume import CloudVolume
import numpy as np
from scipy import ndimage

mip = 6
src = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices', mip=mip)
dst = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/roicc', 
													cdn_cache=True, mip=mip)

# Get bounding box of a total slice
offset = src.voxel_offset
size = tuple(src.shape[:3])
x_slice = slice(offset[0], offset[0]+size[0])
y_slice = slice(offset[1], offset[1]+size[1])
z_slice = slice(offset[2], offset[2]+size[2])
cc_filter = np.ones((3,3))

# Find connected component
for z in range(z_slice.start, z_slice.stop):
	src_image = src[x_slice, y_slice, z]
	src_image = np.reshape(src_image, src_image.shape[:2])
	labeled, nr_objects = ndimage.label(src_image <= 1, structure=cc_filter)
	labeled = (labeled != labeled[0,0]).astype(np.uint8)
	dst[x_slice, y_slice, z] = np.reshape(labeled, labeled.shape+(1,))
