from cloudvolume import CloudVolume
import numpy as np
import pandas as pd

# key slices @ mip=0
film01 = (slice(218500, 219524), slice(54000, 55024))
film02 = (slice(213400, 214424), slice(78200, 79224))
film03 = (slice(216400, 217424), slice(105600, 106624))
film04 = (slice(217700, 218724), slice(59900, 60924))

resin01 = (slice(212300, 213324), slice(54000, 55024))
resin02 = (slice(211600, 212624), slice(78200, 79224))
resin03 = (slice(213800, 214824), slice(105600, 106624))
resin04 = (slice(209200, 210224), slice(60900, 61924))

z_range = range(301, 460)

src_mip = 0
dst_mip = 5
img = CloudVolume("gs://neuroglancer/basil_v0/father_of_alignment/v3", 
												mip=dst_mip, cdn_cache=False)

def scale_slice(sl, dst_mip=dst_mip, src_mip=src_mip):
	scale = 2**(dst_mip - src_mip)
	dst_sl = ()
	for s in sl[0:2]:
		dst_sl += (slice(s.start/scale, s.stop/scale), )
	return dst_sl + (sl[2], )

def cv_max(sl, ):
	return np.max(img[sl])

def cv_min(sl):
	return find_nonzero_min(img[sl])

def find_nonzero_min(arr):
	b = arr[arr != 0]
	if len(b) == 0:
		return 0
	else:
		return np.min(b)

# output dataframe
columns = ['film_min', 'resin_max']
df = pd.DataFrame(index=z_range, columns=columns)

for z in df.index:
	film = [film01 + (z,), film02 + (z,), film03 + (z,), film04 + (z,)]
	resin = [resin01 + (z,), resin02 + (z,), resin03 + (z,), resin04 + (z,)]
	film = map(scale_slice, film)
	resin = map(scale_slice, resin)
	df.loc[z,'film_min'] = np.array(map(cv_min, film))
	df.loc[z,'resin_max'] = np.array(map(cv_max, resin))