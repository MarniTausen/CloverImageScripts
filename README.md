Clover genotype analysis supplementary
================
11/06/2018 - 15:54:04

-   [Introduction](#introduction)
-   [Image Analysis](#image-analysis)

Introduction
============

This repository contains scripts for the clover and rhizobia image analysis pipeline. It also includes scripts of the analysis part of the pipeline.

Image Analysis
==============

The main image filtering script containing the generation of the masks.

``` python
from sys import argv
from skimage import io
from skimage.color import hsv2rgb, rgb2hsv, lab2rgb, rgb2lab, gray2rgb
from skimage.measure import label
from skimage import exposure
import numpy as np
from time import clock
from optparse import OptionParser

# C settings hue=[0.085, 0.36], sat=[0.15, 0.9], val=[0, 1]
def LABMaskImage(rgb, L=[10, 90], a=[-128, -4], b=[4, 128]):
    LABimg = rgb2lab(rgb)
    r = np.zeros(LABimg.shape[:2], dtype=int)
    minv = [L[0], a[0], b[0]]
    maxv = [L[1], a[1], b[1]]
    tmask = LABimg>=minv
    tmask = (LABimg<=maxv)==tmask
    mask = np.zeros(LABimg.shape[:2], dtype=bool)
    mask = tmask.sum(2)==3
    r[:] = 255
    r[mask] = 0
    return r

# C settings hue=[0.085, 0.36], sat=[0.15, 0.9], val=[0, 1]
def HSVMaskImage(rgb, hue=[0.085, 0.42], sat=[0, 1], val=[0, 1]):
    HSVimg = rgb2hsv(rgb)
    r = np.zeros(HSVimg.shape[:2], dtype=int)
    minv = [hue[0], sat[0], val[0]]
    maxv = [hue[1], sat[1], val[1]]
    tmask = HSVimg>=minv
    tmask = (HSVimg<=maxv)==tmask
    mask = np.zeros(HSVimg.shape[:2], dtype=bool)
    mask = tmask.sum(2)==3
    r[:] = 255
    r[mask] = 0
    return r

def combinemasks(mask1, mask2):
    combined = mask1+mask2
    combined[combined>255] = 255
    return combined

def savesummary(mask, img, outfile):
    n = mask[:]
    n[n==0] = 1
    n[n==255] = 0
    TotalArea = n.sum()
    w, h = mask.shape
    PArea = TotalArea/float(w*h)
    MeanColor = img[mask==0].mean()

    f = open(outfile, "w")
    f.write("TotalArea\t%Area\tMean\n")
    f.write("%i\t%f\t%f" % (TotalArea, PArea, MeanColor))
    f.close()

if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-i', type="string", nargs=1, dest="input", help="input image")
    parser.add_option('-o', type="string", nargs=1, dest="output", help="output file")
    options, args = parser.parse_args()

    if options.input==None:
        raise "No input image given, please give input image: -i image"
    else:
        file_name = options.input

    if options.output==None:
        out_file = "overlay_"+file_name.split("/")[-1]
    else:
        out_file = options.output
        
    inimg = io.imread(file_name)
    labmask = LABMaskImage(inimg)
    hsvmask = HSVMaskImage(inimg)
    mask = combinemasks(labmask, hsvmask)
    #io.imsave(out_file, mask)
    overlay = inimg
    overlay[mask==0] = (55, 255, 55)
    io.imsave(out_file, overlay)
    #savesummary(mask, inimg, argv[1]+".xls")
```
