from cloudvolume import CloudVolume
import numpy as np

roi5 = CloudVolume("gs://neuroglancer/pinky100_v0/father_of_alignment_v3/roi", mip=5)
roi6 = CloudVolume("gs://neuroglancer/pinky100_v0/father_of_alignment_v3/roi", mip=6)

# dst5 = CloudVolume("gs://neuroglancer/pinky100_v0/aligned_test_v28", mip=5, cdn_cache=False)
# dst6 = CloudVolume("gs://neuroglancer/pinky100_v0/aligned_test_v28", mip=6, cdn_cache=False)
dst5 = roi5
dst6 = roi6

def get_xy_slice(cv):
	o = cv.voxel_offset
	s = cv.shape
	return slice(o[0], o[0]+s[0]), slice(o[1], o[1]+s[1])

z_sections  = [207]
z_sections += range(1193, 1196)
z_sections += range(305, 315)
z_sections += range(544, 552)

# copy roi to dst
def copy_src_to_dst(z_sections):
	roi5_xy_slice = get_xy_slice(roi5)
	roi6_xy_slice = get_xy_slice(roi6)
	for z in z_sections:
		print(z)
		sl5 = roi5_xy_slice + (z,)
		sl6 = roi6_xy_slice + (z,)
		dst5[sl5] = roi5[sl5]
		dst6[sl6] = roi6[sl6]

def scale_slice(sl, dst_mip=5, src_mip=0):
	scale = 2**(dst_mip - src_mip)
	dst_sl = ()
	for s in sl[0:2]:
		dst_sl += (slice(s.start/scale, s.stop/scale), )
	return dst_sl + (sl[2], )

def create_zero_array(sl):
	"""Create 2D array of zeros
	"""
	x = sl[0].stop - sl[0].start
	y = sl[1].stop - sl[1].start
	return np.zeros((x, y, 1)).astype(np.uint8)

def round_slice_to_chunk(sl, chunks, offset):
	"""Find smallest chunk-aligned slice containing slice

	Args:
		sl: tuple of two slice objects and an index (x_slice, y_slice, z)
		chunks: collection of chunk dimensions
		offset: minimum coordinate of all chunks

	Output:
		tuple of two slice objects and an index (x_slice, y_slice, z)
	"""	
	x_start = int((np.floor((sl[0].start - offset[0])*1.0 / chunks[0])) * chunks[0] + offset[0])
	x_stop  = int((np.ceil((sl[0].stop - offset[0])*1.0 / chunks[0])) * chunks[0] + offset[0])
	y_start = int((np.floor((sl[1].start - offset[1])*1.0 / chunks[1])) * chunks[1] + offset[1])
	y_stop  = int((np.ceil((sl[1].stop - offset[1])*1.0 / chunks[1])) * chunks[1] + offset[1])
	return slice(x_start, x_stop), slice(y_start, y_stop), sl[2]

def round_slice_to_cv_chunk(cv, sl):
	"""Find smallest chunk-aligned slice containing slice, using cv properties
	"""
	chunks = cv.underlying
	offset = cv.voxel_offset
	return round_slice_to_chunk(sl, chunks, offset)

def subtract_slice(sl1, sl2):
	"""Subtract sl2 from sl1 (sl1 - sl2), assume sl2 contained in sl1 & zero z
	"""
	new_sl = ()
	for (s1, s2) in zip(sl1[:2], sl2[:2]):	
		new_sl += (slice(s1.start - s2.start, s1.stop - s2.start),)
	return new_sl + (0,)

def fill_cv_zeros(sl, src_mip=0):
	"""Fill slice range with zeros at appropriate mip levels
	"""
	for (mip, cv) in [(5, dst5), (6, dst6)]:
		scaled_slice = scale_slice(sl, dst_mip=mip, src_mip=src_mip)
		rounded_slice = round_slice_to_cv_chunk(cv, scaled_slice)
		data = cv[rounded_slice]
		local_slice = subtract_slice(scaled_slice, rounded_slice)
		data[local_slice] = create_zero_array(local_slice)
		cv[rounded_slice] = data


# all slices defined in mip 0 coordinates
for z in range(305, 314):
	print(z)
	fix305 = (slice(114600, 118000), slice(47045, 51563), z)
	fill_cv_zeros(fix305)

for z in range(544, 552):
	print(z)
	fix544 = (slice(42000, 59201), slice(47647, 78000), z)
	fill_cv_zeros(fix544)

fix207_01 = (slice(63784, 65900), slice(71312, 77414), 207)
fix207_02 = (slice(61248, 63786), slice(72651, 77414), 207)
fix207_03 = (slice(57892, 61250), slice(74422, 77414), 207)
fix207_04 = (slice(56397, 57894), slice(75466, 77414), 207)
fix1193 = (slice(114118, 117940), slice(39294, 43411), 1193)
fix1194 = (slice(111662, 118000), slice(38433, 79676), 1194)
fix1195_01 = (slice(111000, 118600), slice(37793, 56456), 1195)
fix1195_02 = (slice(109393, 115000), slice(37793, 44273), 1195)
fix1195_03 = (slice(108237, 109400), slice(37793, 40847), 1195)
fixes = [fix207_01,
	fix207_02,
	fix207_03,
	fix207_04,
	fix1193,
	fix1194,	
	fix1195_01,
	fix1195_02,
	fix1195_03]
for f in fixes:
	print(f)
	fill_cv_zeros(f)