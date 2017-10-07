from cloudvolume import CloudVolume, Storage
from shapely.geometry import Polygon, box
from shapely.geometry.collection import GeometryCollection
import numpy as np
from io import BytesIO
from PIL import Image, ImageDraw
import json
import sys

vol = None
mask = None
contrast = None
order = None
out = None

class PreprocessChunkTask():
	"""Run the following image processing tasks on a per chunk basis:
		* Stretch the contrast (requires precomputed quantiles: see contrast.py)
		* Apply edge masks (requires point clouds)
		* Remove & transpose slices (requires dst_to_src mapping)

	Those processes are hard-coded in that order.
	"""
	def __init__(self, inlayer, outlayer, contrastlayer, masklayer, orderlayer,
						x_slice, y_slice, z_slice):
		global vol
		global mask
		global contrast
		global order
		global out
		if vol is None:
			vol = CloudVolume(inlayer)
		if mask is None:
			mask = Storage(masklayer)
		if contrast is None:
			contrast = Storage(contrastlayer)
		if order is None:
			order = Storage(orderlayer)
		if out is None:
			out = CloudVolume(outlayer)
		self.vol = vol
		self.mask = mask
		self.contrast = contrast
		self.order = order
		self.out = out
		self.x_slice = slice(*x_slice)
		self.y_slice = slice(*y_slice)
		self.dst_z_slice = slice(*z_slice)
		self.name  = str(self.x_slice.start)+'-'+str(self.x_slice.stop)+'_'
		self.name += str(self.y_slice.start)+'-'+str(self.y_slice.stop)+'_'
		self.name += str(self.dst_z_slice.start)+'-'+str(self.dst_z_slice.stop)

		self.bbox = box(self.x_slice.start, self.y_slice.start, 
										self.x_slice.stop, self.y_slice.stop)

		f = self.order.get_file('dst_to_src.csv')
		order_arr = np.genfromtxt(BytesIO(f), dtype=np.int, 
											delimiter=',', skip_header=1)
		self.dst_to_src = {order_arr[i,0]:order_arr[i,1] 
									for i in range(order_arr.shape[0])}
		src_z_range = []
		for dst_z in range(self.dst_z_slice.start, self.dst_z_slice.stop):
			if dst_z in self.dst_to_src:
				src_z_range.append(self.dst_to_src[dst_z])
		self.src_z_slice = slice(min(src_z_range), max(src_z_range)+1)

		self.src_chunk = self.vol[self.x_slice, self.y_slice, self.src_z_slice]
		self.src_shape = self.src_chunk.shape
		self.src_chunk = np.reshape(self.src_chunk, self.src_shape[:3])
		self.dst_chunk = np.zeros((self.x_slice.stop-self.x_slice.start,
						self.y_slice.stop-self.y_slice.start,
						self.dst_z_slice.stop-self.dst_z_slice.start), np.uint8)


		z_range = list(range(self.src_z_slice.start, self.src_z_slice.stop))
		z_str = list(map(str, z_range))
		contrast_files = self.contrast.get_files(z_str)
		self.quantiles = {}
		for f in contrast_files:
			if f['content'] is not None:
				arr = np.fromstring(f['content'], sep='\n').astype(np.int)
				self.quantiles[int(f['filename'])] = arr
		self.mask_polygons = {}
		contrast_files = self.mask.get_files(z_str)
		for f in contrast_files:
			if f['content'] is not None:
				pts = np.genfromtxt(BytesIO(f['content']), 
										dtype=np.int, delimiter=',')
				poly = Polygon(map(tuple, pts[:,:2]))
				self.mask_polygons[int(f['filename'])] = poly

	def upload(self):
		"""Upload dst_chunk to gcloud
		"""
		x, y, z = self.x_slice, self.y_slice, self.dst_z_slice
		self.out[x,y,z] = np.reshape(self.dst_chunk, self.dst_chunk.shape+(1,))

	def stretch_contrast(self, im, v_min, v_max, minval=1/255., maxval=1.):
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

	def correct_contrast_src_chunk(self):
		"""Autostretch contrast between 5 and 95% quantile values
		Only correct if contrast file is available
		"""
		for src_z in range(self.src_z_slice.start, self.src_z_slice.stop):
			if src_z in self.quantiles:
				z = src_z - self.src_z_slice.start
				im = self.src_chunk[:,:,z]
				v_min = self.quantiles[src_z][1]  #  5% -> 0
				v_max = self.quantiles[src_z][-2] # 95% -> 1
				self.src_chunk[:,:,z] = self.stretch_contrast(im, v_min, v_max)

	def get_section_polygon_mask(self, z):
		"""Retrieve points list from gcloud for given slice & create polygon
		"""
		if z in self.mask_polygons:
			return self.mask_polygons[z]
		else:
			return Polygon([])

	def get_src_chunk_mask(self, z):
		"""For slice, get binary mask of pixels contained in section polygon
		"""
		sm = self.get_section_polygon_mask(z)
		cm = sm.intersection(self.bbox)
		if type(cm) is Polygon:
			pts = [(int(a[0])-self.x_slice.start, 
				int(a[1])-self.y_slice.start) for a in list(cm.exterior.coords)]
			mask = Image.new('1', self.src_shape[:2], 0)
			ImageDraw.Draw(mask).polygon(pts, outline=1, fill=1)
			return np.transpose(np.array(mask))
		else:
			return np.zeros(self.src_shape[:2], dtype=bool)

	def mask_src_section(self, z):
		"""Apply make to src_chunk for slice z
		"""
		mask = self.get_src_chunk_mask(z)
		z_off = z-self.src_z_slice.start
		self.src_chunk[:,:,z_off] = np.multiply(self.src_chunk[:,:,z_off], mask)

	def mask_src_chunk(self):
		"""Apply slice-wise mask to src_chunk
		"""
		for z in range(self.src_z_slice.start, self.src_z_slice.stop):
			self.mask_src_section(z)

	def copy_src_to_dst(self):
		"""Copy over the slices from src to dst as specified by dst_to_src
		"""
		for z in range(self.dst_z_slice.start, self.dst_z_slice.stop):
			if z in self.dst_to_src:
				src_z = self.dst_to_src[z] - self.src_z_slice.start
				dst_z = z - self.dst_z_slice.start
				self.dst_chunk[:,:,dst_z] = self.src_chunk[:,:,src_z]

	def execute(self):
		print('preprocessing {0}'.format(self.name))
		self.correct_contrast_src_chunk()
		self.mask_src_chunk()
		self.copy_src_to_dst()
		self.upload()

def test():
	inlayer = 'gs://neuroglancer/pinky100_v0/image'
	outlayer = 'gs://neuroglancer/pinky100_v0/image_corrected'
	contrastlayer = 'gs://neuroglancer/pinky100_v0/contrast'
	masklayer = 'gs://neuroglancer/pinky100_v0/edge_mask'
	orderlayer = 'gs://neuroglancer/pinky100_v0/z_order_corrected'
	x_slice = (78272-24*64, 78272-0*64)
	y_slice = (76800-12*64, 76800-0*64)
	z_slice = (129, 193)
	c = PreprocessChunkTask(inlayer, outlayer, contrastlayer, 
						masklayer, orderlayer,
						x_slice, y_slice, z_slice)
	c.execute()

def main(arg):
	params = json.loads(arg)
	task = PreprocessChunkTask(**params)
	task.execute()

if __name__ == "__main__":
	main(sys.argv[1])