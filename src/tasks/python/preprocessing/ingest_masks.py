from cloudvolume import CloudVolume, Storage
from shapely.geometry import Polygon, box
import numpy as np
from io import BytesIO
from PIL import Image, ImageDraw

mip = 6
mask = Storage('gs://neuroglancer/pinky100_v0/edge_mask')
order = Storage('gs://neuroglancer/pinky100_v0/z_order_corrected')
out = CloudVolume('gs://neuroglancer/pinky100_v0/image_single_slices/roi', 
													cdn_cache=True, mip=mip)

# Get bounding box of a total slice
offset = out.voxel_offset
size = tuple(out.shape[:2])
x_slice = slice(offset[0], offset[0]+size[0])
y_slice = slice(offset[1], offset[1]+size[1])

# Build z remap dict
f = order.get_file('dst_to_src.csv')
order_arr = np.genfromtxt(BytesIO(f), dtype=np.int, 
											delimiter=',', skip_header=1)
dst_to_src = {order_arr[i,0]:order_arr[i,1] 
									for i in range(order_arr.shape[0])}

# Compile all ROI mask polygons (indexed by src_z)
# ROI is translated by offset
src_z_range = dst_to_src.values()
src_z_filenames = list(map(str, src_z_range))
mask_files = mask.get_files(src_z_filenames)
mask_polygons = {}
for f in mask_files:
	if f['content'] is not None:
		pts = np.genfromtxt(BytesIO(f['content']), 
								dtype=np.float, delimiter=',')
		poly = Polygon(map(tuple, pts[:,:2] / 2**mip))
		pts = [(int(round(a[0]))-offset[0], int(round(a[1]))-offset[1]) for a 
												in list(poly.exterior.coords)]
		mask_polygons[int(f['filename'])] = pts

# Create mask image for each dst_z
dst_z_range = dst_to_src.keys()
for dst_z in dst_z_range:
	src_z = dst_to_src[dst_z]
	pts = mask_polygons[src_z]
	img_mask = Image.new('1', size, 0)
	ImageDraw.Draw(img_mask).polygon(pts, outline=1, fill=1)
	img_mask = np.transpose(np.array(img_mask)).astype(np.uint8)
	out[x_slice, y_slice, dst_z] = np.reshape(img_mask, img_mask.shape+(1,))
