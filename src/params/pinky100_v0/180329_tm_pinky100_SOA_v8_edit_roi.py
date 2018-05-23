from cloudvolume import CloudVolume
import numpy as np

roi6 = CloudVolume("gs://neuroglancer/pinky100_v0/son_of_alignment_v8/roi_eroded", mip=6, cdn_cache=False)
dst6 = roi6

def get_xy_slice(cv):
	o = cv.voxel_offset
	s = cv.shape
	return slice(o[0], o[0]+s[0]), slice(o[1], o[1]+s[1])

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

def create_one_array(sl):
	"""Create 2D array of zeros
	"""
	x = sl[0].stop - sl[0].start
	y = sl[1].stop - sl[1].start
	return np.ones((x, y, 1)).astype(np.uint8)

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
	for (mip, cv) in [(6, dst6)]:
		scaled_slice = scale_slice(sl, dst_mip=mip, src_mip=src_mip)
		rounded_slice = round_slice_to_cv_chunk(cv, scaled_slice)
		data = cv[rounded_slice]
		local_slice = subtract_slice(scaled_slice, rounded_slice)
		data[local_slice] = create_zero_array(local_slice)
		cv[rounded_slice] = data

def fill_cv_ones(sl, src_mip=0):
	"""Fill slice range with zeros at appropriate mip levels
	"""
	for (mip, cv) in [(6, dst6)]:
		scaled_slice = scale_slice(sl, dst_mip=mip, src_mip=src_mip)
		rounded_slice = round_slice_to_cv_chunk(cv, scaled_slice)
		data = cv[rounded_slice]
		local_slice = subtract_slice(scaled_slice, rounded_slice)
		data[local_slice] = create_one_array(local_slice)
		cv[rounded_slice] = data


# all slices defined in mip 0 coordinates
for z in range(1, 10):
	f = (slice(54986, 70501), slice(36858, 41293), z)
	print(f)
	fill_cv_ones(f)

fix1185 = (slice(41831, 44190), slice(59456, 62398), 1185)
fix1013 = (slice(96279, 99986), slice(38010, 40144), 1013)
fix851 = (slice(93110, 96378), slice(38346, 39967), 851)
fix603 = (slice(50967, 52257), slice(51639, 54243), 603)
fix533 = (slice(70962, 81476), slice(37306, 40791), 533)
fix533_0 = (slice(70860, 77578), slice(40098, 41775), 533)

for z in range(436, 534):
	f = (slice(114000, 118120), slice(41883, 79698), z)
	print(f)
	fill_cv_ones(f)

for z in range(135, 137):
	f = (slice(56243, 67231), slice(37097, 39747), z)
	print(f)
	fill_cv_ones(f)

for z in range(136, 435):
	f = (slice(52972, 67664), slice(35941, 39242), z)
	print(f)
	fill_cv_ones(f)

for z in range(1, 106):
	f = (slice(52899, 63708), slice(34321, 38007), z)
	print(f)
	fill_cv_ones(f)

fix95 = (slice(52899, 116176), slice(34321, 38450), 95)