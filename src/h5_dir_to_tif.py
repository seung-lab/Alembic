#!/usr/bin/python
"""
T Macrina
160314

Make TIF images of all H5 files in directory
**All H5 files must have "img" group

Args:
	sys.argv[1]: full path to the H5 image directory

Returns:
	TIF files (extension changed to .tif) saved in the same directory
"""

from PIL import Image
import numpy as np
import h5py
import os
import sys

def h5_to_array(fn):
	"""Open H5 file with "img" group & convert to numpy ndarray of dtype

	Args:
		fn: filename (full path) of the image

	Returns:
		An ndarray of dtype
	"""
	f = h5py.File(fn, "r")
	return np.array(f["/img"]).T

def write_to_tif(fn, arr):
	"""Write ndarray to tif file
	"""
	img = Image.fromarray(arr)
	img.save(fn)

def main():
	"""Make TIF images of all H5 matrices in directory
	"""
	dir = os.getcwd()
	if len(sys.argv) > 1:
		dir = sys.argv[1]
	files = os.listdir(dir)
	files.sort()
	for file in files:
		if file.endswith(".h5") or file.endswith(".hdf5"):
			fn = os.path.join(dir, file)
			array = h5_to_array(fn)
			new_fn = os.path.splitext(fn)[0] + ".tif"
			write_to_tif(new_fn, array)

if __name__ == '__main__':
	main()