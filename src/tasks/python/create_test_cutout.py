from cloudvolume import CloudVolume

image_in = 'gs://neuroglancer/pinky100_v0/image_single_slices'
image_out = 'gs://neuroglancer/pinky100_v0/aligned_test_v12'
image_mip = 0
roi_in = 'gs://neuroglancer/pinky100_v0/image_single_slices/roicc'
roi_out = 'gs://neuroglancer/pinky100_v0/aligned_test_v11/roicc'
roi_mip = 6
cfsplit_in = 'gs://neuroglancer/pinky100_v0/image_single_slices/cfsplit'
cfsplit_out = 'gs://neuroglancer/pinky100_v0/aligned_test_v11/cfsplit'
cfsplit_mip = 2
cfmanual_in = 'gs://neuroglancer/pinky100_v0/image_single_slices/cfmanual'
cfmanual_out = 'gs://neuroglancer/pinky100_v0/aligned_test_v11/cfmanual'
cfmanual_mip = 5
match_in = 'gs://neuroglancer/pinky100_v0/image_single_slices/nccnet'
match_out = 'gs://neuroglancer/pinky100_v0/aligned_test_v11/nccnet'
match_mip = 2
dst_in = 'gs://neuroglancer/pinky100_v0/aligned_test_v5'
dst_mip = 0

src_dst = [(cfmanual_in, cfmanual_out, cfmanual_mip)]

z_slice = slice(199, 208)
src_mip = 0

def scale_slice(s, src_mip, dst_mip):
	scale = 1/2**(dst_mip - src_mip)
	return slice(int(s.start*scale), int(s.stop*scale))

def scale_slices(x_slice, y_slice, z_slice, src_mip, dst_mip):
	return (scale_slice(x_slice, src_mip, dst_mip), 
			scale_slice(y_slice, src_mip, dst_mip), 
			scale_slice(z_slice, src_mip, dst_mip))

def get_cloudvolume(path, mip): 
	return CloudVolume(path, mip=mip)

def update_info_mips(cv, no_of_mips=6):
	print("updating info mips")
	for mip in range(1,no_of_mips+1):
		factor = (2**mip, 2**mip, 1)
		cv.add_scale(factor)
		cv.commit_info()

def get_xy_slice(cv):
	o = cv.voxel_offset
	s = cv.shape
	return slice(o[0], o[0]+s[0]), slice(o[1], o[1]+s[1])

for (src_path, dst_path, mip) in src_dst:
	print(src_path)
	print(dst_path)
	print(mip)
	cv = get_cloudvolume(dst_path, 0)
	update_info_mips(cv, 6)
	dst_cv = get_cloudvolume(dst_path, mip)
	src_cv = get_cloudvolume(src_path, mip)
	sl = get_xy_slice(dst_cv) + (z_slice,)
	print(sl)
	dst_cv[sl] = src_cv[sl]