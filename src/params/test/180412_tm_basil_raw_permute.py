from cloudvolume import CloudVolume

src_path = 'gs://neuroglancer/basil_v0/raw_image_cropped'
dst_path = 'gs://neuroglancer/basil_v0/raw_image_cropped/permuted'

def get_xy_slice(cv):
	o = cv.voxel_offset
	s = cv.shape
	return slice(o[0], o[0]+s[0]), slice(o[1], o[1]+s[1])

for mip in range(9,-1,-1):
	src_cv = CloudVolume(src_path, mip=mip, fill_missing=True, parallel=64)
	dst_cv = CloudVolume(dst_path, mip=mip, cdn_cache=False)
	print(mip)
	for src_z in range(133,134):
		dst_z = src_z - 4
		src_slice = get_xy_slice(src_cv) + (src_z,)
		dst_slice = get_xy_slice(dst_cv) + (dst_z,)
		print((src_slice, dst_slice))
		dst_cv[dst_slice] = src_cv[src_slice]
	for src_z in range(129,133):
		dst_z = src_z + 1
		src_slice = get_xy_slice(src_cv) + (src_z,)
		dst_slice = get_xy_slice(dst_cv) + (dst_z,)
		print((src_slice, dst_slice))
		dst_cv[dst_slice] = src_cv[src_slice]
	for src_z in [127,128,135,136]:
		dst_z = src_z
		src_slice = get_xy_slice(src_cv) + (src_z,)
		dst_slice = get_xy_slice(dst_cv) + (dst_z,)
		print((src_slice, dst_slice))
		dst_cv[dst_slice] = src_cv[src_slice]