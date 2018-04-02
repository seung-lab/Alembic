from cloudvolume import CloudVolume
import numpy as np
import pandas as pd
from random import randint
from os.path import expanduser, join
import tiffile

def get_xy_slice(cv):
	o = cv.voxel_offset
	s = cv.shape
	return slice(o[0], o[0]+s[0]), slice(o[1], o[1]+s[1])

label_track_cv = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/cfmanual', mip=5)
label_cv = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/cfmanual', mip=4)
image_cv = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices', mip=5)

def defect_unique_count(img, z_range):
	columns = [0, 1, 2, 3]
	section_df = pd.DataFrame(index=z_range, columns=columns);
	section_df.fillna(0, inplace=True)

	for z in z_range:
		print(z)
		z_slice = img[get_xy_slice(img) + (z,)]
		vals, counts = np.unique(z_slice, return_counts=True)
		for (v, c) in zip(vals, counts):
			section_df.loc[z, v] = c
	return section_df

# section_df = defect_unique_count(label_track_cv, range(1,2176))
# section_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_counts.csv'))
section_df = pd.read_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_counts.csv'), index_col=0)

def get_random_chunk(cv, z, chunk_dim=(256,256)):
	o = cv.voxel_offset
	s = cv.shape
	x_start = randint(o[0], o[0]+s[0]-chunk_dim[0])
	y_start = randint(o[1], o[1]+s[1]-chunk_dim[1])
	x_stop = x_start+chunk_dim[0]
	y_stop = y_start+chunk_dim[1]
	return slice(x_start, x_stop), slice(y_start, y_stop), z

def defect_chunk_counts(cv, section_df; chunks_per_section=50, chunk_dim=(256, 256)):
	crack_df = section_df[section_df['1'] > 0]
	no_of_chunks = chunks_per_section*len(crack_df.index)
	columns = ['x_start', 'y_start', 'x_stop', 'y_stop', 'z', '0', '1', '2', '3']
	chunk_df = pd.DataFrame(index=range(no_of_chunks), columns=columns)
	chunk_df.fillna(0, inplace=True)

	for (k, z) in enumerate(crack_df.index):
		for n in range(chunks_per_section):
			i = k*chunks_per_section + n
			print((k,z,n,i))
			sl = get_random_chunk(cv, z, chunk_dim=chunk_dim)
			chunk = cv[sl]
			vals, counts = np.unique(chunk, return_counts=True)
			chunk_df.loc[i, 'x_start'] = sl[0].start
			chunk_df.loc[i, 'x_stop'] = sl[0].stop
			chunk_df.loc[i, 'y_start'] = sl[1].start
			chunk_df.loc[i, 'y_stop'] = sl[1].stop
			chunk_df.loc[i, 'z'] = sl[2]
			for (v, c) in zip(vals, counts):
				chunk_df.loc[i, str(v)] = c
	return chunk_df

chunks_per_section = 50
chunk_dim = (256, 256)
# chunk_df = defect_chunk_counts(label_track_cv, section_df, chunks_per_section=chunks_per_section, chunk_dim=chunk_dim)
# chunk_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_random_chunk_counts01.csv'))
# chunk_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_random_chunk_counts02.csv'))
# chunk_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_random_chunk_counts03.csv'))
# chunk_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_random_chunk_counts04.csv'))

chunk_filenames = [join(expanduser('~'), 'seunglab/defect_mask_net/data/pinky100_cfmanual_random_chunk_counts{:02d}.csv'.format(i)) for i in range(1,5)]
chunk_dfs = [pd.read_csv(fn, index_col=0) for fn in chunk_filenames]
chunk_df = pd.concat(chunk_dfs, ignore_index=True)
chunk_df['positive_ground_truth'] = chunk_df['1'] > 0
chunk_df['negative_ground_truth'] = False
chunk_df.loc[chunk_df[(chunk_df['2'] > 0) & (chunk_df['1'] == 0)].sample(n=365).index, 'negative_ground_truth'] = True
chunk_df.loc[chunk_df[(chunk_df['2'] == 0) & (chunk_df['1'] == 0)].sample(n=366).index, 'negative_ground_truth'] = True
# chunk_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/180401_pinky100_cfmanual_ground_truth_log.csv'))
chunk_df = pd.read_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/180401_pinky100_cfmanual_ground_truth_log.csv'), index_col=0)

ground_truth_df = chunk_df[chunk_df['negative_ground_truth'] | (chunk_df['positive_ground_truth'])]
ground_truth_df.reset_index(inplace=True, drop=True)
ground_truth_df.to_csv(join(expanduser('~'), 'seunglab/defect_mask_net/data/180402_pinky100_ground_truth_only.csv'))
positive_gt = ground_truth_df[ground_truth_df['positive_ground_truth']].reset_index(inplace=True, drop=True)
negative_gt = ground_truth_df[ground_truth_df['positive_ground_truth']].reset_index(inplace=True, drop=True)

def compile_ground_truth(cv, df, src_mip=5, dst_mip=5):
	s = 2**(src_mip - dst_mip)
	gt_shape = (chunk_dim[0]*s, chunk_dim[1]*s)
	gt = np.ndarray(gt_shape + (len(df.index),), np.uint8)
	for (k, row) in df.iterrows():
		sl = slice(row['x_start']*s, row['x_stop']*s), slice(row['y_start']*s, row['y_stop']*s), row['z']
		print((k, sl))
		gt[:,:,k] = cv[sl].reshape(gt_shape)


tiffile.imsave(join(expanduser('~'), 'seunglab/defect_mask_net/data/180401_pinky100_positive_gt_image.tif'), gt)
