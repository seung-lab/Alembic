from cloudvolume import CloudVolume, Storage
import numpy as np
from tqdm import tqdm
from scipy.stats.mstats import mquantiles
import io
import csv
import sys

class ContrastQuantiles():
	"""Build quantiles for each section for contrast correction task to follow
	"""
	def __init__(self, inlayer, outlayer, x_slice, y_slice, z_slice, n_samples):
		self.vol = CloudVolume(inlayer)
		self.out = Storage(outlayer)
		self.min_x = x_slice.start
		self.max_x = x_slice.stop
		self.min_y = y_slice.start
		self.max_y = y_slice.stop
		self.min_z = z_slice.start
		self.max_z = z_slice.stop
		self.n_samples = n_samples

	def round_to_chunk(self, x, origin_x, chunk_x):
		"""Round coordinate up to the nearest chunk
		"""
		chunks = np.ceil((x - origin_x) / chunk_x).astype(np.int)
		return origin_x + chunks*chunk_x

	def round_point(self, x, y):
		"""Round point up to the nearest chunk
		"""
		origin = self.vol.voxel_offset
		chunk_dims = self.vol.underlying
		x_new = self.round_to_chunk(x, origin[0], chunk_dims[0])
		y_new = self.round_to_chunk(y, origin[1], chunk_dims[1])
		return x_new, y_new

	def get_random_point(self):
		"""Get randome point within the x_slice, y_slice ranges of the class
		"""
		chunk_dims = self.vol.underlying
		x = np.random.randint(self.min_x, self.max_x-chunk_dims[0]-1)
		y = np.random.randint(self.min_y, self.max_y-chunk_dims[1]-1)
		return self.round_point(x, y)

	def accept_sample(self, smp, zero_threshold=0.1):
		"""Accept samples that are not (zero_threshold)% black pixels
		"""
		chunk_dims = self.vol.underlying
		n_pixels = chunk_dims[0]*chunk_dims[1]
		return np.sum(smp == 0) < zero_threshold

	def compile_samples(self, z):
		"""Collect n_samples of chunk intensities from section z
		"""
		chunk_dims = self.vol.underlying
		samples = []
		k = 0
		with tqdm(total=100) as pbar:
			while k < self.n_samples:
				try:
					x, y = self.get_random_point()
					x_slice = slice(x, x+chunk_dims[0])
					y_slice = slice(y, y+chunk_dims[1])
					smp = self.vol[x_slice, y_slice, z].flatten()
					if self.accept_sample(smp):
						k += 1
						pbar.update(100/self.n_samples)
						samples.append(smp)
					else:
						print('chunk not accepted')
				except:
					print(sys.exc_info()[0])
		return np.concatenate(samples)

	def get_quantiles(self, samples, bins=20):
		"""Create quantiles of collected samples with n bins
		"""
		n = 1.0/bins
		prob = np.arange(0.0,1.0,n)
		return mquantiles(samples, prob=prob)

	def get_quantiles_and_upload(self, z):
		samples = self.compile_samples(z)
		quantiles = self.get_quantiles(samples)
		fn = '{0}'.format(z)
		output = io.BytesIO()
		np.savetxt(output, quantiles, delimiter=',', newline='\n')
		self.out.put_file(fn, output.getvalue())

	def loop_through_z(self):
		for z in tqdm(range(self.min_z, self.max_z)):
			print(z)
			self.get_quantiles_and_upload(z)


class ContrastQuantilesChunked():
	"""Build quantiles for each section for contrast correction task to follow
	"""
	def __init__(self, inlayer, outlayer, x_slice, y_slice, z_slice, n_samples):
		self.vol = CloudVolume(inlayer, cache=True)
		self.out = Storage(outlayer)
		self.min_x = x_slice.start
		self.max_x = x_slice.stop
		self.min_y = y_slice.start
		self.max_y = y_slice.stop
		self.min_z = z_slice.start
		self.max_z = z_slice.stop
		self.n_samples = n_samples
		self.chunk_locations = []
		self.chunks = {}

	def round_to_chunk(self, x, origin_x, chunk_x):
		"""Round coordinate up to the nearest chunk
		"""
		chunks = np.ceil((x - origin_x) / chunk_x).astype(np.int)
		return origin_x + chunks*chunk_x

	def round_point(self, x, y):
		"""Round point up to the nearest chunk
		"""
		origin = self.vol.voxel_offset
		chunk_dims = self.vol.underlying
		x_new = self.round_to_chunk(x, origin[0], chunk_dims[0])
		y_new = self.round_to_chunk(y, origin[1], chunk_dims[1])
		return x_new, y_new

	def get_random_point(self):
		"""Get randome point within the x_slice, y_slice ranges of the class
		"""
		chunk_dims = self.vol.underlying
		x = np.random.randint(self.min_x, self.max_x-chunk_dims[0]-1)
		y = np.random.randint(self.min_y, self.max_y-chunk_dims[1]-1)
		return self.round_point(x, y)

	def accept_sample(self, smp, zero_threshold=0.1):
		"""Accept samples that are not (zero_threshold)% black pixels
		"""
		chunk_dims = self.vol.underlying
		n_pixels = chunk_dims[0]*chunk_dims[1]
		return np.sum(smp == 0) < zero_threshold

	def select_chunk_locations(self):
		self.chunk_locations = []
		for k in range(self.n_samples):
			self.chunk_locations.append(self.get_random_point())

	def select_chunks(self, z):
		chunk_dims = self.vol.underlying
		self.chunks = {}
		for x,y in tqdm(self.chunk_locations):
			x_slice = slice(x, x+chunk_dims[0])
			y_slice = slice(y, y+chunk_dims[1])
			z_slice = slice(z, z+chunk_dims[2])
			self.chunks[x,y] = self.vol[x_slice, y_slice, z_slice]

	def compile_samples(self, chunk_z, z):
		"""Collect n_samples of chunk intensities from section z
		"""
		samples = []
		for x,y in self.chunk_locations:
			chunk = self.chunks[x,y]
			smp = chunk[:,:,z-chunk_z].flatten()
			if self.accept_sample(smp):
				samples.append(smp)
			else:
				print('sample not accepted')
		return samples

	def get_quantiles(self, samples, bins=20):
		"""Create quantiles of collected samples with n bins
		"""
		n = 1.0/bins
		prob = np.arange(0.0,1.0,n)
		return mquantiles(samples, prob=prob)

	def get_quantiles_and_upload(self, chunk_z, z):
		samples = self.compile_samples(chunk_z, z)
		if len(samples) > 0:
			quantiles = self.get_quantiles(np.concatenate(samples))
			fn = '{0}'.format(z)
			output = io.BytesIO()
			np.savetxt(output, quantiles, delimiter=',', newline='\n')
			self.out.put_file(fn, output.getvalue())

	def loop_through_z(self):
		chunk_dims = self.vol.underlying
		self.select_chunk_locations()
		z_range = range(self.min_z, self.max_z, chunk_dims[2])
		for z_start, z_stop in zip(z_range[:-1], z_range[1:]):
			self.select_chunks(z_start)
			for z in range(z_start, z_stop):
				print(z)
				self.get_quantiles_and_upload(z_start, z)


class CorrectContrastTask():
	"""Use precomputed quantiles for each section to adjust the contrast of each
	section within a chunk.
	"""

	def __init__(self, inlayer, outlayer, contrastlayer, x_slice, y_slice, z_slice):
		self.vol = CloudVolume(inlayer)
		self.out = CloudVolume(outlayer)
		self.contrast = Storage(contrastlayer)
		self.x_slice = x_slice
		self.y_slice = y_slice
		self.z_slice = z_slice
		self.chunk = self.vol[x_slice, y_slice, z_slice]

	def get_quantiles(self):
		"""Return dict of quantile arrays, indexed by z
		Files are in contrast_layer named by their z index
		"""
		z_range = list(range(self.z_slice.start, self.z_slice.stop))
		z_str = list(map(str, z_range))
		files = self.contrast.get_files(z_str)
		quantiles = {}
		for f in files:
			if f['content'] is not None:
				arr = np.fromstring(f['content'], sep='\n', dtype=np.int) #.astype(np.int)
				quantiles[int(f['filename'])] = arr
		return quantiles

	def stretch_contrast_section(self, im, v_min, v_max, 
													minval=1/255., maxval=1.):
		"""Linearly transform im pixel values to fit within v_min and v_max

		Args:
			im: np.uint8 image
			v_min: lower intensity threshold (anything below set to minval)
			v_max: upper intensity threshold (anything above set to maxval)
			minval: lower bound of output intensities (0 reserved for masks)
			maxval: upper bound of output intensities

		Retruns:
			np.uint8 image
		"""
		img = im / 255.
		vmin = v_min / 255.
		vmax = v_max / 255.
		img = np.minimum(maxval, np.maximum(minval, (img-vmin)/(vmax-vmin)))
		return (img*255).astype(np.uint8)

	def correct_and_upload(self):
		quantiles = self.get_quantiles()
		for z in range(self.z_slice.start, self.z_slice.stop):
			chunk_z = z-self.z_slice.start
			im = self.chunk[:,:,chunk_z]
			v_min = quantiles[z][1]
			v_max = quantiles[z][-2]
			self.chunk[:,:,chunk_z] = self.stretch_contrast_section(im, 
																v_min, v_max)
		self.out[self.x_slice, self.y_slice, self.z_slice] = self.chunk



def calculate_contrast_quantiles():
	inlayer = 'gs://neuroglancer/pinky100_v0/image'
	outlayer = 'gs://neuroglancer/pinky100_v0/contrast'
	x_slice = slice(60000,100000)
	y_slice = slice(50000, 70000)
	z_slice = slice(1153,   1218)
	c = ContrastQuantilesChunked(inlayer, outlayer, x_slice, y_slice, z_slice, 200)
	c.loop_through_z()

def correct_contrast():
	inlayer = 'gs://neuroglancer/pinky100_v0/image'
	outlayer = 'gs://neuroglancer/pinky100_v0/image_corrected'
	contrastlayer = 'gs://neuroglancer/pinky100_v0/contrast'
	x_slice = slice(65920+64*4,65920+64*8)
	y_slice = slice(55488,55488+64*4)
	z_slice = slice(1,65)
	c = CorrectContrastTask(inlayer, outlayer, contrastlayer, 
											x_slice, y_slice, z_slice)
	c.correct_and_upload()

def loop_over_xy():
	inlayer = 'gs://neuroglancer/pinky100_v0/image'
	outlayer = 'gs://neuroglancer/pinky100_v0/image_corrected'
	contrastlayer = 'gs://neuroglancer/pinky100_v0/contrast'
	z_slice = slice(1,65)

	x_range = range(65920,68000,512)
	y_range = range(55488,55488+64*10,512)
	for x_start, x_stop in zip(x_range[:-1], x_range[1:]):
		for y_start, y_stop in zip(y_range[:-1], y_range[1:]):
			x_slice = slice(x_start, x_stop)
			y_slice = slice(y_start, y_stop)
			c = CorrectContrastTask(inlayer, outlayer, contrastlayer, 
												x_slice, y_slice, z_slice)
			c.correct_and_upload()
