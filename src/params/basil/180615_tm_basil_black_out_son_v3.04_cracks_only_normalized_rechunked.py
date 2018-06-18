from cloudvolume import CloudVolume
import numpy as np

cv_path = "gs://neuroglancer/basil_v0/son_of_alignment/v3.04_cracks_only_normalized_rechunked"

def scale_slice(sl, dst_mip=1, src_mip=0):
	scale = 2**(dst_mip - src_mip)
	dst_sl = ()
	for s in sl[0:2]:
		dst_sl += (slice(int(s.start/scale), int(s.stop/scale)), )
	return dst_sl + (sl[2], )

def create_zero_array(sl):
	"""Create 2D array of zeros
	"""
	x = sl[0].stop - sl[0].start
	y = sl[1].stop - sl[1].start
	return np.zeros((x, y, 1)).astype(np.uint8)

def fill_cv_zeros(sl, src_mip=0):
	"""Fill slice range with zeros at appropriate mip levels
	"""
	mips = range(1,7)
	for mip in mips:
		print(mip)
		cv = CloudVolume(cv_path, mip=mip, cdn_cache=False, non_aligned_writes=True)
		scaled_slice = scale_slice(sl, dst_mip=mip, src_mip=src_mip)
		print(scaled_slice)
		cv[scaled_slice] = create_zero_array(scaled_slice)

# all slices defined in mip 0 coordinates
fills = [(slice(220000, 220200), slice(231400, 232000), 101),
		# (slice(151848, 180000), slice(237283, 267446), 101),
		# (slice(180000, 235773), slice(237283, 267446), 101),
		# (slice(157658, 233612), slice(206646, 237903), 101),
		# (slice(199542, 227829), slice(211022, 233612), 101),
		# (slice(176099, 201417), slice(186642, 214773), 101),
		# (slice(200818, 227242), slice(175060, 208089), 101),
		(slice(183323, 186907), slice(122283, 127407), 104),
		(slice(57837, 81391), slice(17484, 29396), 186),
		(slice(39319, 59624), slice(14019, 61559), 186),
		(slice(35999, 39923), slice(45358, 47298), 186),
		(slice(121859, 134318), slice(78120, 126053), 196),
		(slice(114098, 130260), slice(95204, 109530), 196),
		(slice(140584, 173482), slice(37046, 76586), 198) ,
		(slice(136755, 174732), slice(75570, 105264), 198),
		(slice(134271, 154527), slice(115146, 149375), 201),
		(slice(40451, 52240), slice(71818, 103644), 228),
		(slice(184826, 203190), slice(6114, 77067), 294),
		(slice(168770, 199570), slice(68913, 262740), 294),
		(slice(146998, 170761), slice(85773, 262740), 294),
		(slice(124694, 146998), slice(99979, 262873), 294),
		(slice(104117, 124827), slice(124008, 263005), 294),
		(slice(78229, 104117), slice(154277, 263138), 294),
		(slice(67188, 72100), slice(12379, 77032), 354),
		(slice(17977, 55311), slice(55971, 99093), 874)]

for f in fills:
	print(f)
	fill_cv_zeros(f)