from PIL import Image
import numpy as np
import h5py
import os
import sys
import cv2

# Maybe consider implemeting more involved auto-balancing
# http://wiki.cmci.info/documents/120206pyip_cooking/python_imagej_cookbook#automatic_brightnesscontrast_button

def open(fn):
	print "Opening " + fn
	ext = os.path.splitext(fn)
	if ext == '.h5':
		f = h5py.File(fn, 'r')
		img = f['/img']
		return np.array(img)
	else:
		img = Image.open(fn)
		return np.array(img)

def write_to_h5(fn, arr):
	"""Write ndarray to H5 file under group "main"
	"""
	print "Writing to " + fn
	sz = np.asarray(arr.shape)
	f = h5py.File(fn, "w")
	f.create_dataset("/img", data=arr.T, dtype=arr.dtype)
	f.create_dataset("/size", data=sz, dtype=sz.dtype)
	f.close()

def main(src_dir, dst_dir):
	"""Make H5 files of all TIF images in directory
	"""
	# src_dir = "/media/tmacrina/667FB0797A5072D7/3D_align"
	# dst_dir = "/media/tmacrina/4BED39E032CF5004/datasets/AIBS_import/3_prealigned"
	# src_dir = "/usr/people/tmacrina/seungmount/research/Julimaps/datasets/AIBS_aligned"	
	# dst_dir = "/usr/people/tmacrina/seungmount/research/Julimaps/datasets/AIBS_aligned/3_prealigned"	
	# src_dir = "/media/tmacrina/4BED39E032CF5004/datasets/AIBS_import/0_overview"
	# dst_dir = "/media/tmacrina/4BED39E032CF5004/datasets/stage_stitch/3_prealigned"

	files = os.listdir(src_dir)
	for file in files:
		if file.endswith(".tif"):
			print "Importing " + file
			src_fn = os.path.join(src_dir, file)
			dst_fn = os.path.join(dst_dir, "1," + str(file[2:6]) + "_montaged.h5")
			# dst_fn = os.path.join(dst_dir, "1," + str(file[1:5]) + "_prealigned.h5")
			if os.path.isfile(src_fn) and not os.path.isfile(dst_fn):
				# print(src_fn, dst_fn)
				arr = open(src_fn)
				write_to_h5(dst_fn, arr)

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])