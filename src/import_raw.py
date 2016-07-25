from PIL import Image
import numpy as np
import h5py
import os
import sys
import cv2

# Maybe consider implemeting more involved auto-balancing
# http://wiki.cmci.info/documents/120206pyip_cooking/python_imagej_cookbook#automatic_brightnesscontrast_button

def open(fn):
	ext = os.path.splitext(fn)
	if ext == '.h5'
		f = h5py.File(fn, 'r')
		img = f['/img']
		return np.array(img)
	else
		img = Image.open(fn)
		return np.array(img)

def auto_adjust(img, minpercent=5, maxpercent=95):
	minval = np.percentile(img[:], minpercent)
	maxval = np.percentile(img[:], maxpercent)
		distr = nquantile(img[:], 16)
		minval = distr[2]; maxval = distr[16];
		return min(1, max(0, (img-minval) / (maxval-minval)))


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

def main():
	"""Make TIF images of all H5 matrices in directory
	"""
	dir = os.getcwd()
	# file = sys.argv[1]
	files = os.listdir(dir)
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(63,63))
	for file in files:
		if file.endswith("1,1_prealigned.h5"):
			print "Applying CLAHE to " + file
		# if file == 'Tile_r1-c7_S2-W001_sec15.h5':
			fn = os.path.join(dir, file)
			apply_clahe_to_H5(fn, clahe)

if __name__ == '__main__':
	main()