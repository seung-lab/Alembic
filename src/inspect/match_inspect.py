import os
import numpy as np
import pandas as pd
import correspondences

from scipy import spatial

class MatchInspect():
	"""Object to control display of correspondence point pairs with filters so
	they can be manually inspected and editted.

	df: pandas dataframe with pre & post points along with current filters
	points_controller: points.Controller object for state server interaction
	"""
	def __init__(self, name, path, controller, mip, scale):
		self.name = name
		self.path = path
		self.mip = mip
		self.scale = scale
		self.c = controller
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

	def np_to_ng_list(self, arr):
		return list(zip(map(int,arr[:,0]*(2**self.mip)), 
	                    map(int,arr[:,1]*(2**self.mip)), 
	                    map(int,arr[:,2])))

	def get_synapses(self, use_filter=True):
		"""Return point pairs for synapses, scaled accordingly, along w/ indices
		"""
		sub_df = self.df.copy()
		if use_filter:
			sub_df = sub_df[sub_df['filter']]
		post = sub_df[['post_x','post_y','post_z']].as_matrix()
		pre = sub_df[['pre_x','pre_y','pre_z']].as_matrix()
		post[:,:2] = pre[:,:2] + self.scale*(post[:,:2]-pre[:,:2])
		pts = zip(self.np_to_ng_list(pre), self.np_to_ng_list(post))
		return pts, sub_df.index
	
	def syn_to_ng(self, use_filter=True):
		syn, _ = self.get_synapses(use_filter=use_filter)
		# list of interleaved pre & post tuples
		return [p for s in syn for p in s]

	def ng_to_syn(self):
		syn = self.c.get()
		return [(syn[i], syn[i+1]) for i in range(0,len(syn),2)]

	def set_filter(self, bool_arr):
		self.df['filter'] = bool_arr
		self.display()

	def invert(self):
		print('Inverting filter')
		self.set_filter(~self.df['filter'])

	def display(self):
		self.c.set(self.syn_to_ng(use_filter=True))
		self.c.set_z(self.df['pre_z'][0])

	def update_df_from_ng(self):
		"""Determine which synapses were removed from neuroglancer
			**depends on neuroglancer maintaining order
		"""
		new_syn = self.ng_to_syn()
		all_syn, indices = self.get_synapses(use_filter=True)
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
				if all_syn[i] != new_syn[i-j]:
					removed_indices.append(indices[i])
					j += 1
			i += 1
		self.df.loc[removed_indices, 'filter'] = False

	def get_nearest(self):
		"""Return the index of the correspondence that's nearest to the position
			of neuroglancer.
		"""
		ctr = np.array(self.c.get_position()[:2])
		print('Finding match nearest to {0}'.format(ctr))
		pt =  ctr / (2**self.mip)
		fpts = self.df[self.df['filter']]
		dist, i = spatial.KDTree(fpts[['pre_x', 'pre_y']]).query(pt)
		print('{0}: {1} @ ({2},{3})'.format(self.name, fpts.iloc[i].name, 
							*fpts.iloc[i][['pre_x', 'pre_y']]*(2**self.mip)))
		return fpts.iloc[i]

	def save(self):
		self.update_df_from_ng()
		df = self.df.copy()
		df['filter'] = df['filter'].astype('float')
		df.to_csv(self.path, index=False)



class MeshSetInspect():
	"""Inspect correspondences of the match objects in a meshset.
	The MeshSet must have been serialized to local disk (use `to_csv(ms)`).
	"""
	def __init__(self, dir, port=9999, mip=0, scale=1, start_id=0):
		self.controller = correspondences.controller(port)
		self.mip = mip
		self.scale = scale
		self.id = start_id
		self.load(dir)

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
		filepath = os.path.join(self.dir, fn)
		print('Inspecting {0}: {1} ({2}/{3})'.format(self.id, fn, 
													self.id+1, len(self.files)))
		return MatchInspect(fn, filepath, self.controller, self.mip, self.scale)

	def get_next(self):
		return self.go_to(self.id+1)

	def get_prev(self):
		return self.go_to(self.id-1)

	def save(self):
		self.match_obj.save()

	def save_next(self):
		self.save()
		return self.get_next()
