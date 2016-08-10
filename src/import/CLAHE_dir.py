from PIL import Image
import numpy as np
import h5py
import os
import sys
import cv2

# Maybe consider implemeting more involved auto-balancing
# http://wiki.cmci.info/documents/120206pyip_cooking/python_imagej_cookbook#automatic_brightnesscontrast_button

def apply_clahe_to_H5(fn, clahe):
	f = h5py.File(fn, "r+")
	img = f["/img"]
	# apply clahe
	arr = clahe.apply(np.array(img))
	# stretch distribution across 0-255 range
	max_a = np.max(arr)
	min_a = np.min(arr)
	alpha = 255.0/(max_a - min_a)
	beta = -alpha*min_a
	arr = (alpha*arr + beta).astype(np.uint8)
	# resave image
	img[...] = arr
	f.close()

def get_H5_array(fn):
	f = h5py.File(fn, "r")
	return np.array(f["/img"])

def main():
	"""Make TIF images of all H5 matrices in directory
	"""
	dir = os.getcwd()
	# file = sys.argv[1]
	files = os.listdir(dir)
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(127,127))
	for file in files:
		if file.endswith("1,1_prealigned.h5"):
			print "Applying CLAHE to " + file
		# if file == 'Tile_r1-c7_S2-W001_sec15.h5':
			fn = os.path.join(dir, file)
			apply_clahe_to_H5(fn, clahe)

# if __name__ == '__main__':
# 	main()