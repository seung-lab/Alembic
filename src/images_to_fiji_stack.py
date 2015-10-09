# images_to_fiji_stack.py
# Load image directory into FIJI, apply translations, and merge into stack

import os
import csv
from ij import IJ, ImagePlus, ImageStack

resize_factor = 0.4
sec_start = 11
sec_end = 40

bucket_dir_path = '/usr/people/tmacrina/seungmount'
datasets_dir_path = "research/Julimaps/datasets"
cur_dataset = "piriform"

premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"

premontaged_offsets_filename = 'premontaged_offsets.txt'
montaged_offsets_filename = 'montaged_offsets.txt'
prealigned_offsets_filename = 'prealigned_offsets.txt'
aligned_offsets_filename = 'aligned_offsets.txt'

dir_path = os.path.join(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path, "1,1-1,40_aligned")
offsets_csv = open(os.path.join(dir_path, aligned_offsets_filename))
image_reader = csv.reader(offsets_csv, delimiter='\t')

image_stack = None

for idx, row in enumerate(image_reader):
  if idx+1 in range(sec_start,sec_end+1):
    filename = row[0]
    # i_offset = int(row[1])
    # j_offset = int(row[2])
    path = os.path.join(dir_path, filename)
    imp = IJ.openImage(path)
    name = imp.getTitle()
    ip = imp.getProcessor()
    ip.setInterpolationMethod(ImageProcessor.BILINEAR)
    # ip.translate(j_offset, i_offset)
    ip = ip.resize(int(imp.getWidth()*resize_factor)) # will scale & crop
    if image_stack is None:
        image_stack = ImageStack(ip.width, ip.height)
    image_stack.addSlice(name, ip)
    # new_imp = ImagePlus(name, ip)
    # new_imp.show()

imp = ImagePlus("stack", image_stack)
imp.show()

img_width = imp.width
img_height = imp.height

win_width = 400
win_height = 400

width = img_width/win_width
height = img_height/win_height

# index for tile to inspect
n = 0

def n2xy(n):
  return (n/height)*win_width, (n%height)*win_height

def start():
  go_to(0, 0)
  return 0

def go_next(n=n):
  n += 1
  min(n, width*height)
  go_to(*n2xy(n))
  return n

def go_back(n=n):
  n -= 1
  max(n, 0)
  go_to(*n2xy(n))
  return n

def go_to(x, y, w=win_width, h=win_height):
  x = max(x, 0)
  y = max(y, 0)
  x = min(x, img_width - w)
  y = min(y, img_height - h)

  IJ.makeRectangle(x, y, w, h)
  IJ.log(str((x, y, w, h)))
  IJ.run("To Selection")

# def inspect_movie(imp):
#   win = imp.getWindow()
#   c = win.getCanvas()