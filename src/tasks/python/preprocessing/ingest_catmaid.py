from __future__ import division
from cloudvolume import CloudVolume, Storage
from PIL import Image
import numpy as np
import math
import os

class IngestCATMAID():

	def __init__(self, catmaid_dir, inlayer, outlayer, 
									tile_x=1024, tile_y=1024, z_offset=2266):
		self.catmaid_dir = catmaid_dir
		self.tile_x = tile_x
		self.tile_y = tile_y
		self.z_offset = z_offset
		self.vol = CloudVolume(inlayer)
		self.out = CloudVolume(outlayer)

	def floor_chunk(self, x, origin_x, chunk_x):
		"""Round coordinate up to the nearest chunk
		"""
		chunks = np.floor((x - origin_x) / chunk_x).astype(np.int)
		return origin_x + chunks*chunk_x

	def floor_point(self, x, y, z):
		"""Round point up to the nearest chunk
		"""
		origin = self.vol.voxel_offset
		chunk_dims = self.vol.underlying
		x_new = self.floor_chunk(x, origin[0], chunk_dims[0])
		y_new = self.floor_chunk(y, origin[1], chunk_dims[1])
		z_new = self.floor_chunk(z, origin[2], chunk_dims[2])
		return x_new, y_new, z_new

	def ceil_chunk(self, x, origin_x, chunk_x):
		"""Round coordinate up to the nearest chunk
		"""
		chunks = np.ceil((x - origin_x) / chunk_x).astype(np.int)
		return origin_x + chunks*chunk_x

	def ceil_point(self, x, y, z):
		"""Round point up to the nearest chunk
		"""
		origin = self.vol.voxel_offset
		chunk_dims = self.vol.underlying
		x_new = self.ceil_chunk(x, origin[0], chunk_dims[0])
		y_new = self.ceil_chunk(y, origin[1], chunk_dims[1])
		z_new = self.ceil_chunk(z, origin[2], chunk_dims[2])
		return x_new, y_new, z_new

	def ng_to_catmaid_range(self, x, y, z_range):
		"""Take ng x,y,z_range and return col,row,slice_range for CATMAID
		"""
		return (int(math.floor(y/self.tile_y)), 
			int(math.floor(x/self.tile_x)), np.array(z_range)+self.z_offset)

	def ng_to_catmaid(self, x, y, z):
		"""Take ng x,y,z and return col,row,slice for CATMAID
		"""
		return (int(math.floor(y/self.tile_y)), 
					int(math.floor(x/self.tile_x)), z+self.z_offset)

	def catmaid_to_ng(self, c, r, s):
		"""Take CATMAID col,row,slice to ng x,y,z
		"""
		return self.tile_y*r, self.tile_x*c, s-self.z_offset

	def catmaid_to_ng_range(self, c, r, s_range):
		"""Take CATMAID col,row,slice_range to ng x,y,z chunk-aligned slices
		"""
		start = self.catmaid_to_ng(c, r, s_range[0])
		stop = self.catmaid_to_ng(c+1, r+1, s_range[-1])
		x_start, y_start, z_start = self.floor_point(*start)
		x_stop, y_stop, z_stop = self.ceil_point(*stop)
		x_slice = slice(x_start, x_stop)
		y_slice = slice(y_start, y_stop)
		z_slice = slice(z_start, z_stop)
		return x_slice, y_slice, z_slice

	def get_catmaid_url(self, c, r, s):
		return '{0}/0/{1}/{2}/{3}.png'.format(self.catmaid_dir, s, c, r)

	def get_local_dir(self):
		home = os.path.expanduser('~')
		return '{0}/catmaid/'.format(home)
	
	def get_local_url(self, c, r, s):
		local_dir = self.get_local_dir()
		return '{0}/{3}_{1}_{2}.png'.format(local_dir, c, r, s)

	def make_local_dir(self):
		path = self.get_local_dir()
		if not os.path.exists(path):
			os.makedirs(path)

	def get_catmaid_files(self, x, y, z_range):
		self.make_local_dir()
		for z in z_range:
			c, r, s = self.ng_to_catmaid(x,y,z)
			catmaid_url = self.get_catmaid_url(c,r,s)
			local_url = self.get_local_url(c,r,s)
			os.system('gsutil -m cp {0} {1}'.format(catmaid_url, local_url))

	def make_catmaid_stack(self, x, y, z_range):
		cstack = np.zeros((self.tile_x, self.tile_y, len(z_range)), np.uint8)
		for k, z in enumerate(z_range):
			c, r, s = self.ng_to_catmaid(x,y,z)
			local_url = self.get_local_url(c,r,s)
			cstack[:,:,k] = np.transpose(np.array(Image.open(local_url)))
		return cstack

	def upload_catmaid_stack(self, x, y, z_range):
		c, r, s_range = self.ng_to_catmaid_range(x, y, z_range)
		xs_start, ys_start, zs_start = self.catmaid_to_ng(c, r, s_range[0])
		xs_stop, ys_stop, zs_stop = self.catmaid_to_ng(c+1, r+1, s_range[-1]+1)
		x_slice, y_slice, z_slice = self.catmaid_to_ng_range(c, r, s_range)
		ngvol = self.vol[x_slice, y_slice, z_slice]
		ngvol = ngvol.reshape(ngvol.shape[0],ngvol.shape[1],ngvol.shape[2])
		cstack = self.make_catmaid_stack(x, y, z_range)
		x_adj = slice(xs_start-x_slice.start, xs_stop-x_slice.start)
		y_adj = slice(ys_start-y_slice.start, ys_stop-y_slice.start)
		z_adj = slice(zs_start-z_slice.start, zs_stop-z_slice.start)
		ngvol[x_adj,y_adj,z_adj] = cstack
		self.out[x_slice, y_slice, z_slice] = ngvol


def main():
	# x = 107965
	# y = 57923
	# z_range = range(256,272)

	x = 63744
	y = 37056
	z_range = range(385,450)

	# x = 71229
	# y = 62974
	# z_range = range(572,573)

	catmaid_dir = 'gs://aibs_alignment/20170927_PINKY100_BOSSPRECURSOR'
	inlayer = 'gs://neuroglancer/pinky100_v0/image_corrected'
	outlayer = 'gs://neuroglancer/pinky100_v0/image_test'
	ic = IngestCATMAID(catmaid_dir, inlayer, inlayer)
	ic.get_catmaid_files(x, y, z_range)
	ic.upload_catmaid_stack(x, y, z_range)

if __name__ == "__main__":
    main()