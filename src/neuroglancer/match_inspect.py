import os
import numpy as np
import pandas as pd
import correspondences

class MatchInspect():
	"""Object to control display of correspondence point pairs with filters so
	they can be manually inspected and editted.

	df: pandas dataframe with pre & post points along with current filters
	points_controller: points.Controller object for state server interaction
	"""
	def __init__(self, name, df, controller):
		print('Inspecting {0}'.format(name))
		self.name = name
		self.df = df
		self.c = controller
		self.display()

	def get_df(self):
		self.update()
		return self.df

	def get_pre(self, sub_df):
		return list(zip(map(int,sub_df.pre_x), 
	                    map(int,sub_df.pre_y), 
	                    map(int,sub_df.pre_z)))

	def get_post(self, sub_df):
		return list(zip(map(int,sub_df.post_x), 
	                    map(int,sub_df.post_y), 
	                    map(int,sub_df.post_z)))

	def get_synapses(self, use_filter=True):
		sub_df = self.df
		if use_filter:
			sub_df = sub_df[sub_df['filter'] == True]
		return zip(self.get_pre(sub_df), self.get_post(sub_df))
	
	def syn_to_ng(self, use_filter=True):
		syn = self.get_synapses(use_filter=use_filter)
		# list of interleaved pre & post tuples
		return [p for s in syn for p in s]

	def ng_to_syn(self):
		syn = self.c.get()
		return [(syn[i], syn[i+1]) for i in range(0,len(syn),2)]

	def flip(self):
		print('Flipping filter')
		accepted = self.df['filter']
		self.df['filter'] = [not x for x in accepted]
		self.display()

	def display(self):
		self.c.set(self.syn_to_ng(use_filter=True))

	def update(self):
		new_syn = self.ng_to_syn()
		all_syn = self.get_synapses(use_filter=False)
		self.df['filter'] = [syn in new_syn for syn in all_syn]
		self.display()



class MeshSetInspect():
	"""Inspect correspondences of the match objects in a meshset.
	The MeshSet must have been serialized to local disk (use `to_csv(ms)`).
	"""
	def __init__(self, dir, port=9999):
		self.dir = dir
		self.files = os.listdir(dir)
		self.files.sort()
		self.controller = correspondences.controller(port)
		self.id = 0
		self.match_obj = None

	def set_id(self, i):
		if (i < 0) or (i >= len(self.files)):
			raise ValueError('index outside of files')
		else:
			self.id = i

	def get_next(self):
		self.set_id(self.id+1)
		return self.inspect()

	def get_fn(self):
		return self.files[self.id]

	def get_path(self):
		return os.path.join(self.dir, self.get_fn())

	def inspect(self):
		fn = self.get_fn()
		filepath = self.get_path()
		df = pd.read_csv(filepath)
		df['filter'] = df['filter'].astype('bool')
		self.match_obj = MatchInspect(fn, df, self.controller)
		return self.match_obj

	def save(self):
		df = self.match_obj.get_df()
		df['filter'] = df['filter'].astype('float')
		df.to_csv(self.get_path())

	def save_next(self):
		self.save()
		return self.get_next()







