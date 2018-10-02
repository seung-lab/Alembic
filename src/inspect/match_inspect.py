import os
from os.path import join
import numpy as np
import pandas as pd
from correspondences import Controller


from scipy import spatial

class MatchInspect():
	"""Object to control display of correspondence point pairs with filters so
	they can be manually inspected and editted.

	df: pandas dataframe with pre & post points along with current filters
	points_controller: points.Controller object for state server interaction
	"""
	def __init__(self, path, ng, mip, scale):
		self.path = path 
		self.mip = mip
		self.scale = scale
		self.ng = ng
		self.load()
		self.display()

	def set_scale(self, s):
		self.scale = s
		self.display()

	def reset_scale(self):
		self.set_scale(1)

	def load(self):
		df = pd.read_csv(self.path)
		df['filter'] = df['filter'].astype('bool')
		self.df = df

	def scale_list(self, arr):
		return list(zip(map(int,arr[:,0]*(2**self.mip)), 
	                    map(int,arr[:,1]*(2**self.mip)), 
	                    map(int,arr[:,2])))

	def get_correspondences(self):
		"""Return point pairs for synapses, scaled accordingly, along w/ indices
		"""
		sub_df = self.df.copy()
		sub_df = sub_df[sub_df['filter']]
		post = sub_df[['post_x','post_y','post_z']].values
		pre = sub_df[['pre_x','pre_y','pre_z']].values
		post[:,:2] = pre[:,:2] + self.scale*(post[:,:2]-pre[:,:2])
		pts = list(zip(self.scale_list(pre), self.scale_list(post)))
		return pts, sub_df.index
	
	def set_filter(self, bool_arr):
		self.df['filter'] = bool_arr
		self.display()

	def invert(self):
		print('Inverting filter')
		self.set_filter(~self.df['filter'])

	def display(self):
                self.ng.set_correspondences(self.get_correspondences()[0])

	def pull(self):
		"""Determine which synapses were removed from neuroglancer
			**depends on neuroglancer maintaining order
		"""
		new_syn = self.ng.get_correspondences()
		all_syn, indices = self.get_correspondences()
		# all_syn = self.get_synapses(use_filter=False)
		# new_filter = [syn in new_syn for syn in all_syn]
		removed_indices = []
		i, j = 0, 0
		while i < len(all_syn):
			if i-j >= len(new_syn):
				# print((indices[i], all_syn[i], "NA"))
				removed_indices.append(indices[i])
			else:
				# print((indices[i], all_syn[i], new_syn[i-j]))
				if np.array_equal(all_syn[i], new_syn[i-j]):
					removed_indices.append(indices[i])
					j += 1
			i += 1
		self.df.loc[removed_indices, 'filter'] = False

	def get_nearest(self):
		"""Return the index of the correspondence that's nearest to the position
			of neuroglancer.
		"""
		ctr = np.array(self.ng.get_position()[:2])
		print('Finding match nearest to {0}'.format(ctr))
		pt =  ctr / (2**self.mip)
		fpts = self.df[self.df['filter']]
		dist, i = spatial.KDTree(fpts[['pre_x', 'pre_y']]).query(pt)
		print('{0} @ ({1},{2})'.format(fpts.iloc[i].name, 
							*fpts.iloc[i][['pre_x', 'pre_y']]*(2**self.mip)))
		return fpts.iloc[i]

	def save(self):
		self.pull()
		df = self.df.copy()
		df['filter'] = df['filter'].astype('float')
		df.to_csv(self.path, index=False)



class MeshSetInspect():
	"""Inspect correspondences of the match objects in a meshset.
	The MeshSet must have been serialized to local disk (use `to_csv(ms)`).
	"""
	def __init__(self, data_dir, img_path, mip=0, scale=1, start_id=0):
		self.ng = Controller()
		self.mip = mip
		self.scale = scale
		self.id = start_id
		self.img_path = img_path
		self.load(data_dir)
		self.set_image_path(img_path)

	def set_image_path(self, path):
		self.img_path = path
		layer = {'name': 'image',
			'url': join('precomputed://', self.img_path)}
		self.ng.set_layer(layer)
		print(self.ng.viewer)

	def load(self, dir):
		self.dir = dir
		self.files = os.listdir(dir)
		self.files.sort()

	def go_to(self, i):
		if (i < 0) or (i >= len(self.files)):
			raise ValueError('index outside of files')
		else:
			self.id = i
		fn = self.files[self.id]
		filepath = join(self.dir, fn)
		print('Inspecting {0}: {1} ({2}/{3})'.format(self.id, fn, 
													self.id+1, len(self.files)))
		self.ng.set_z(int(fn.split()[0][1:-1]))
		return MatchInspect(filepath, self.ng, self.mip, self.scale)

	def get_next(self):
		return self.go_to(self.id+1)

	def get_prev(self):
		return self.go_to(self.id-1)

	def save(self):
		self.match_obj.save()

	def save_next(self):
		self.save()
		return self.get_next()
