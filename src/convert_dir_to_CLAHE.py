#!/usr/bin/env jython

from ij import IJ
import os

from mpicbg.ij.clahe import Flat
from ij.process import ImageConverter

# http://fiji.sc/wiki/index.php/Enhance_Local_Contrast_(CLAHE)
# http://fiji.sc/cgi-bin/gitweb.cgi?p=mpicbg.git;a=blob;f=mpicbg/ij/clahe/PlugIn.java;h=663153764493547de560c08ee11f2e6b1e7e1a32;hb=HEAD

dir = "/usr/people/tmacrina/seungmount/research/Julimaps/datasets/AIBS_pilot_v1/0_raw/"

blocksize = 63
histogram_bins = 255
maximum_slope = 3
mask = "*None*"
composite = False
mask = None

# files = os.listdir(dir)
# files.sort()
# for file in files:
#      if file.endswith(".tif")
fn = os.path.join(dir, 'original.tif')
imp = IJ.openImage(fn)
output_fn = os.path.splitext(fn)[0] + "_CLAHE_8bit.tif"
imp = IJ.openImage(fn)
  
Flat.getFastInstance().run( imp, 
                        blocksize,
                        histogram_bins,
                        maximum_slope,
                        mask,
                        composite )
ImageConverter(imp).convertToGray8()

IJ.save(imp, output_fn)
