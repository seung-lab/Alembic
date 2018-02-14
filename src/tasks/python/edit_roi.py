from cloudvolume import CloudVolume
import numpy as np

roi5 = CloudVolume("gs://neuroglancer/pinky100_v0/father_of_alignment_v3/roi", mip=5)
roi6 = CloudVolume("gs://neuroglancer/pinky100_v0/father_of_alignment_v3/roi", mip=6)

dst5 = CloudVolume("gs://neuroglancer/pinky100_v0/aligned_test_v28", mip=5, cdn_cache=False)
dst6 = CloudVolume("gs://neuroglancer/pinky100_v0/aligned_test_v28", mip=6, cdn_cache=False)

def get_xy_slice(cv):
	o = cv.voxel_offset
	s = cv.shape
	return slice(o[0], o[0]+s[0]), slice(o[1], o[1]+s[1])

z_sections  = [207]
z_sections += range(1193, 1196)
z_sections += range(305, 315)
z_sections += range(544, 552)

# copy roi to dst
roi5_xy_slice = get_xy_slice(roi5)
roi6_xy_slice = get_xy_slice(roi6)
for z in z_sections:
	print(z)
	sl5 = roi5_xy_slice + (z,)
	sl6 = roi6_xy_slice + (z,)
	dst5[sl5] = roi5[sl5]
	dst6[sl6] = roi6[sl6]


# all slices defined in mip 0 coordinates
fix207 = (slice(63784, 65900), slice(71312, 76414), 207)

def adjust_slice(sl, dst_mip=5, src_mip=0):
	scale = 2**(dst_mip - src_mip)
	dst_sl = ()
	for s in sl[0:2]:
		dst_sl += (slice(s.start/scale, s.stop/scale), )
	return dst_sl + (sl[2], )

def fill_zeros(sl):
	x = sl[0].stop - sl[0].start
	y = sl[1].stop - sl[1].start
	return np.zeros((x, y, 1, 1)).astype(np.uint8)

def round_slice_to_chunk(cv, sl):
	chunks = cv.underlying
	offset = cv.voxel_offset
	x_start = int((np.floor((sl[0].start - offset[0])*1.0 / chunks[0])) * chunks[0] + offset[0])
	x_stop  = int((np.ceil((sl[0].stop - offset[0])*1.0 / chunks[0])) * chunks[0] + offset[0])
	y_start = int((np.floor((sl[1].start - offset[1])*1.0 / chunks[1])) * chunks[1] + offset[1])
	y_stop  = int((np.floor((sl[1].stop - offset[1])*1.0 / chunks[1])) * chunks[1] + offset[1])
	return slice(x_start, x_stop), slice(y_start, y_stop), sl[2]

def subtract_slice()

def fill_cv_zeros(sl):
	for (mip, cv) in [(5, dst5), (6, dst6)]:
		mip_sl = adjust_slice(sl, dst_mip=mip, src_mip=0)
		sl_rounded = round_slice_to_chunk(cv, mip_sl)
		data = cv[sl_rounded]
		adj_sl = ()
		data