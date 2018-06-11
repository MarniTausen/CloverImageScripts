Clover genotype analysis supplementary
================
11/06/2018 - 16:14:19

-   [Introduction](#introduction)
-   [Image Analysis](#image-analysis)
    -   [Cropping](#cropping)
    -   [Finding Pots](#finding-pots)

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

Cropping
--------

Crop script for specifically doing the crop

``` python
from sys import argv
from skimage import io
from skimage.color import hsv2rgb, rgb2hsv, rgb2grey
from skimage.measure import label
from ImageFiltering import *
from scipy.cluster.vq import kmeans
from skimage.draw import circle, line
from Signature import *
import numpy as np

def drawline(img, r0, c0, r1, c1, t, d, color=(255, 25, 25)):

    thickness = range(0-(t-1), 0+t)

    for t in thickness:
        if d==1:
            rr, cc = line(r0+t, c0, r1+t, c1)
            img[rr,cc] = color
        if d==2:
            rr, cc = line(r0, c0+t, r1, c1+t)
            img[rr,cc] = color

def rectangle(center, width, height, img, t, color=(255, 25, 25)):

    im_height, im_width, _ = img.shape
    
    left = (center[1]-(width/2))
    right = (center[1]+(width/2))
    top = (center[0]-(height/2))
    bottom = (center[0]+(height/2))
    if left<0:
        diff = 0-left
        left += diff+t+1
        right += diff+t+1
    if right>im_width:
        diff = im_width-right
        left += diff-t-1
        right += diff-t-1
    if top<0:
        diff = 0-top
        top += diff+t+1
        bottom += diff+t+1
    if bottom>im_height:
        diff = im_height-bottom
        top += diff-t-1
        bottom += diff-t-1

    # Top line
    drawline(img, top, left, top, right, t=t, d=1, color=color)
    # Left line
    drawline(img, top, left, bottom, left, t=t, d=2, color=color)
    # Bottom line
    drawline(img, bottom, left, bottom, right, t=t, d=1, color=color)
    # Right line
    drawline(img, top, right, bottom, right, t=t, d=2, color=color)
            
def square(center, dim, img, name, t, color=(255, 25, 25)):

    rectangle(center, dim, dim, img, t, color)

    left = (center[1]-(dim/2))
    right = (center[1]+(dim/2))
    top = (center[0]-(dim/2))
    bottom = (center[0]+(dim/2))

    try:
        for off in range(1, name+1):
            drawline(img, top, right-(off*20), top-30, right-(off*20), t=2, d=2, color=color)
    except:
        pass

def crop_square(center, dim, img, i, name):

    im_height, im_width, _ = img.shape
    
    left = (center[1]-(dim/2))
    right = (center[1]+(dim/2))
    top = (center[0]-(dim/2))
    bottom = (center[0]+(dim/2))
    if left<0:
        diff = 0-left
        left += diff+t
        right += diff+t
    if right>im_width:
        diff = im_width-right
        left += diff-t
        right += diff-t
    if top<0:
        diff = 0-top
        top += diff+t
        bottom += diff+t
    if bottom>im_height:
        diff = im_height-bottom
        top += diff-t
        bottom += diff-t
    
    io.imsave("Pot"+str(i)+"_"+name.split("/")[-1],
              img[top:bottom,
                  left:right])
    
def loadcenters(filename):
    f = open(filename)
    centers = []

    #round = lambda x: sorted([int(x), int(x)+1], key=lambda y: abs(y-x))[0]
    
    for line in f:
        info = line.split(" ")[2:]
        center = info[2]
        center = center.split(",")
        center[0] = int(round(float(center[0])))
        center[1] = int(round(float(center[1])))

        if center[0]<100 or center[1]<100:
            continue
        if center[0]>(width-100) or center[1]>(height-100):
            continue
        
        size = [int(i) for i in info[1].split("+")[0].split("x")]

        if size[0]>350:
            centers.append(center[::-1])
            # WARNING?
        else:
            centers.append(center[::-1])

    return centers
    
if __name__=="__main__":
    dim = 400
    inimg = io.imread(argv[1])
    height, width, _ = inimg.shape
    print height, width
    
    #centers = kmeans(blackpixels.T.astype(float), 10, iter=100, thresh=1e-08)[0].astype(int)
    centers = loadcenters(argv[2])
    centers = sorted(centers, key=lambda x: x[0])
    abovecenters = sorted(centers[0:5], key=lambda x: x[1])
    belowcenters = sorted(centers[5:], key=lambda x: x[1])
    print abovecenters
    print belowcenters
    i = 10
    points = {}
    for center in abovecenters:

        square(center, dim, inimg, i, 3)
        points[str(i)] = center
        #crop_square(center, dim, inimg, i, argv[1])
        
        i -= 1
    for center in belowcenters:

        square(center, dim, inimg, i, 3)
        points[str(i)] = center
        #crop_square(center, dim, inimg, i, argv[1])

        i -= 1

    points['0'] = findCentre(points)
    square(points['0'], dim/9, inimg, i, 3)
    writeSignature(argv[4], points)
    io.imsave(argv[3], inimg)
```

Supporting crop functions

``` python
def groupCentres(centres):
    groups = []

    for i in range(len(centres)):
        temp_group = [centres[i]]
        for j in range(i, len(centres)):
            if i==j:
                continue
            temp_group.append(centres[j])
            n = len(temp_group)
            if n<3: continue
            if n>10: continue
            temp_group = orderGroup(temp_group)
            gheight, gwidth = dimension_of_group(temp_group, "dimension")
            if len(temp_group)!=n: continue
            if gheight>1300: continue
            if gwidth>2300: continue
            groups.append(orderGroup(temp_group))

    return groups

def closestGroupSignature(groups, signature, radius):

    number_of_matches = 0
    bestgroup = []

    i = 1
    j = 0
    for group in groups:
        psignature = overlaySignature(signature, group)
        match = []
        for point in psignature.values():
            for centre in group:
                if centre[0]<(point[0]+radius) and centre[0]>(point[0]-radius) and centre[1]<(point[1]+radius) and centre[1]>(point[1]-radius):
                    match.append(centre)
                    continue
        matches = len(match)
        if number_of_matches<matches:
            bestgroup = group
            j = i
            number_of_matches = matches
        i += 1

    return bestgroup


def getPointdiff(point1, point2):
    return point1[0]-point2[0], point1[1]-point2[1]

def dimension_of_group(group, return_value="centre"):
    max_y = max(group, key=lambda x: x[0])[0]
    min_y = min(group, key=lambda x: x[0])[0]
    max_x = max(group, key=lambda x: x[1])[1]
    min_x = min(group, key=lambda x: x[1])[1]
    if return_value=="centre" or return_value=="center":
        return [(max_y-min_y)/2+min_y, (max_x-min_x)/2+min_x]
    if return_value=="rectangle":
        return [[min_y, min_x], [min_y, max_x], [max_y, min_x], [max_y, max_x]]
    if return_value=="dimension":
        width = distance([min_y, min_x], [min_y, max_x])
        height = distance([min_y, min_x], [max_y, min_x])
        return height, width
    if return_value=="area":
        width = distance([min_y, min_x], [min_y, max_x])
        height = distance([min_y, min_x], [max_y, min_x])
        return width*height
    return max_y, min_y, max_x, min_x

def overlaySignature(signature, group):

    centre = dimension_of_group(group, "centre")

    sign_height, sign_width = dimension_of_group(signature.values(), "dimension")
    sign_height = int(sign_height/2)
    sign_width = int(sign_width/2)

    if (centre[0]-sign_height)<0:
        print "Above"
        centre[0] += abs(centre[0]-sign_height)+dim/2
    if (centre[0]+sign_height)>height:
        print "Below"
        centre[0] -= abs(height-(centre[0]+sign_height))+dim/2
    if (centre[1]-sign_width)<0:
        print "Left"
        centre[1] += abs(centre[1]-sign_width)+dim/2
    if (centre[1]+sign_width)>width:
        print "Right"
        centre[1] -= abs(width-(centre[1]+sign_width))+dim/2

    diff_y, diff_x = getPointdiff(signature['0'], centre)

    for point in signature:
        if point=='0':
            signature['0'] = centre
            continue
        signature[point] = [signature[point][0]-diff_y, signature[point][1]-diff_x]

    return signature

def splitGroup(group):
    max_y = max(group, key=lambda x: x[0])[0]
    min_y = min(group, key=lambda x: x[0])[0]
    max_x = max(group, key=lambda x: x[1])[1]
    min_x = min(group, key=lambda x: x[1])[1]
    abovecentres = []
    belowcentres = []
    group = sorted(group, key=lambda x: x[1])
    for centre in group:
        if abs(centre[0]-max_y)<abs(centre[0]-min_y):
            belowcentres.append(centre)
        if abs(centre[0]-max_y)>abs(centre[0]-min_y):
            abovecentres.append(centre)
    b_median = belowcentres[len(belowcentres)/2][0]
    a_median = abovecentres[len(abovecentres)/2][0]
    belowcentres = filter(lambda x: x[0]<(b_median+400) and x[0]>(b_median-400), belowcentres)
    abovecentres = filter(lambda x: x[0]<(a_median+400) and x[0]>(a_median-400), abovecentres)
    return abovecentres, belowcentres

def orderGroup(group):
    abovecentres, belowcentres = splitGroup(group)
    return abovecentres+belowcentres

def control_edge(point, cheight, cwidth):
    left = (point[1]-(cwidth/2))
    right = (point[1]+(cwidth/2))
    top = (point[0]-(cheight/2))
    bottom = (point[0]+(cheight/2))
    if left<0:
        diff = 0-left
        left += diff+1
        right += diff+1
    if right>width:
        diff = width-right
        left += diff-1
        right += diff-1
    if top<0:
        diff = 0-top
        top += diff+1
        bottom += diff+1
    if bottom>height:
        diff = height-bottom
        top += diff-1
        bottom += diff-1
    return top, bottom, left, right
```

Finding Pots
------------

Function for finding the main points on the screen

``` python
from sys import argv
from skimage import io
from skimage.color import hsv2rgb, rgb2hsv, rgb2grey
from skimage.measure import label
from ImageFiltering import *
from Signature import *
from CropPoints import *
import numpy as np
from optparse import OptionParser
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage import label

def loadcenters(filename, mask, potsize=500):
    f = open(filename)
    centers = []

    for line in f:
        info = line.split(" ")[2:]
        center = info[2]
        center = center.split(",")
        center[0] = int(round(float(center[0])))
        center[1] = int(round(float(center[1])))
        center = center[::-1]

        if center[0]<150 or center[1]<150:
            print "Centre on the edge, and will therefore be removed"
            continue
        if center[1]>(width-150) or center[0]>(height-150):
            print "Centre on the edge, and will therefore be removed"
            continue

        size = [int(i) for i in info[1].split("+")[0].split("x")]

        if size[0]>potsize:
            #centers.append(center)
            print "WARNING: Cluster size too big, potential merge of pots detected"
            print "Employing K-means on the merged blob to find the pots"

            c = (size[0]//potsize)+1

            top, bottom, left, right = control_edge(center, potsize, size[0]+100)
            npoint = local_Kmeans(mask[top:bottom,left:right], c)
            temp_centers = []
            for point in npoint:
                diff_x, diff_y = point[0]-int(potsize/2), point[1]-int((size[0]+100)/2)
                point =  [center[0]+diff_x, center[1]+diff_y]
                if point[0]<150:
                    continue
                if point[1]<150:
                    continue
                if point[0]>(height-150):
                    continue
                if point[1]>(width-150):
                    continue
                temp_centers.append(point)

            min_distance = 300
            cgroup = []

            ## Average centers that are to close together.
            for tcenter in temp_centers:
                for tcenter2 in temp_centers:
                    if tcenter==tcenter2: continue
                    dist = distance(tcenter, tcenter2)
                    if dist<min_distance:
                        cgroup = [tcenter, tcenter2]
                        min_distance = dist

            if cgroup==[]:
                for tcenter in temp_centers:
                    centers.append(tcenter)
            else:
                ncenter = findCentre(cgroup)
                centers.append(ncenter)
                for tcenter in temp_centers:
                    if tcenter not in cgroup:
                        centers.append(tcenter)

        else:
            centers.append(center)

    return centers

def groupCentres(centres):
    groups = []

    for i in range(len(centres)):
        temp_group = [centres[i]]
        for j in range(i, len(centres)):
            if i==j:
                continue
            temp_group.append(centres[j])
            n = len(temp_group)
            if n<3: continue
            if n>10: continue
            temp_group = orderGroup(temp_group)
            gheight, gwidth = dimension_of_group(temp_group, "dimension")
            if len(temp_group)!=n: continue
            if gheight>1300: continue
            if gwidth>2300: continue
            groups.append(orderGroup(temp_group))

    return groups

def closestGroupSignature(groups, signature, radius):

    number_of_matches = 0
    bestgroup = []

    i = 1
    j = 0
    for group in groups:
        psignature = overlaySignature(signature, group)
        match = []
        for point in psignature.values():
            for centre in group:
                if centre[0]<(point[0]+radius) and centre[0]>(point[0]-radius) and centre[1]<(point[1]+radius) and centre[1]>(point[1]-radius):
                    match.append(centre)
                    continue
        matches = len(match)
        if number_of_matches<matches:
            bestgroup = group
            j = i
            number_of_matches = matches
        i += 1

    return bestgroup


def getPointdiff(point1, point2):
    return point1[0]-point2[0], point1[1]-point2[1]

def dimension_of_group(group, return_value="centre"):
    max_y = max(group, key=lambda x: x[0])[0]
    min_y = min(group, key=lambda x: x[0])[0]
    max_x = max(group, key=lambda x: x[1])[1]
    min_x = min(group, key=lambda x: x[1])[1]
    if return_value=="centre" or return_value=="center":
        return [(max_y-min_y)/2+min_y, (max_x-min_x)/2+min_x]
    if return_value=="rectangle":
        return [[min_y, min_x], [min_y, max_x], [max_y, min_x], [max_y, max_x]]
    if return_value=="dimension":
        width = distance([min_y, min_x], [min_y, max_x])
        height = distance([min_y, min_x], [max_y, min_x])
        return height, width
    if return_value=="area":
        width = distance([min_y, min_x], [min_y, max_x])
        height = distance([min_y, min_x], [max_y, min_x])
        return width*height
    return max_y, min_y, max_x, min_x

def overlaySignature(signature, group):

    centre = dimension_of_group(group, "centre")

    sign_height, sign_width = dimension_of_group(signature.values(), "dimension")
    sign_height = int(sign_height/2)
    sign_width = int(sign_width/2)

    if (centre[0]-sign_height)<0:
        print "Above"
        centre[0] += abs(centre[0]-sign_height)+dim/2
    if (centre[0]+sign_height)>height:
        print "Below"
        centre[0] -= abs(height-(centre[0]+sign_height))+dim/2
    if (centre[1]-sign_width)<0:
        print "Left"
        centre[1] += abs(centre[1]-sign_width)+dim/2
    if (centre[1]+sign_width)>width:
        print "Right"
        centre[1] -= abs(width-(centre[1]+sign_width))+dim/2

    diff_y, diff_x = getPointdiff(signature['0'], centre)

    for point in signature:
        if point=='0':
            signature['0'] = centre
            continue
        signature[point] = [signature[point][0]-diff_y, signature[point][1]-diff_x]

    return signature

def splitGroup(group):
    max_y = max(group, key=lambda x: x[0])[0]
    min_y = min(group, key=lambda x: x[0])[0]
    max_x = max(group, key=lambda x: x[1])[1]
    min_x = min(group, key=lambda x: x[1])[1]
    abovecentres = []
    belowcentres = []
    group = sorted(group, key=lambda x: x[1])
    for centre in group:
        if abs(centre[0]-max_y)<abs(centre[0]-min_y):
            belowcentres.append(centre)
        if abs(centre[0]-max_y)>abs(centre[0]-min_y):
            abovecentres.append(centre)
    b_median = belowcentres[len(belowcentres)/2][0]
    a_median = abovecentres[len(abovecentres)/2][0]
    belowcentres = filter(lambda x: x[0]<(b_median+400) and x[0]>(b_median-400), belowcentres)
    abovecentres = filter(lambda x: x[0]<(a_median+400) and x[0]>(a_median-400), abovecentres)
    return abovecentres, belowcentres

def orderGroup(group):
    abovecentres, belowcentres = splitGroup(group)
    return abovecentres+belowcentres

def control_edge(point, cheight, cwidth):
    left = (point[1]-(cwidth/2))
    right = (point[1]+(cwidth/2))
    top = (point[0]-(cheight/2))
    bottom = (point[0]+(cheight/2))
    if left<0:
        diff = 0-left
        left += diff+1
        right += diff+1
    if right>width:
        diff = width-right
        left += diff-1
        right += diff-1
    if top<0:
        diff = 0-top
        top += diff+1
        bottom += diff+1
    if bottom>height:
        diff = height-bottom
        top += diff-1
        bottom += diff-1
    return top, bottom, left, right

def local_Kmeans(mask, k=1):
    blackpixels = np.where(mask==0)
    blackpixels = np.matrix(blackpixels)
    if len(blackpixels)==0: return [200, 200]
    centers = kmeans(blackpixels.T.astype(float), k, iter=20, thresh=1e-08)[0].astype(int)
    return centers

def local_center_of_mass(mask):
    blackpixels = np.where(mask==0)
    blackpixels = np.matrix(blackpixels)
    center = np.squeeze(np.asarray(blackpixels.mean(axis=1).astype(int).flatten('C')))
    return center

def claim_region(mask, top, bottom, left, right):
    mask[top:bottom, left:right] = 255

def testRadius(signature, group, radius, mask, dim):
    new_group = []

    del signature['0']

    for point in signature.values():
        match = []
        for centre in group:
            if centre[0]<(point[0]+radius) and centre[0]>(point[0]-radius) and centre[1]<(point[1]+radius) and centre[1]>(point[1]-radius):
                if match==[]:
                    match = centre
                if distance(match, point)>distance(centre, point):
                    match = centre
        if match==[]:
            #new_group.append(point)
            top, bottom, left, right = control_edge(point, dim, dim)
            try:
                npoint = local_Kmeans(mask[top:bottom,left:right])[0]
                npoint = local_center_of_mass(mask[top:bottom,left:right])
            except:
                new_group.append(point)
                continue
            diff_x, diff_y = npoint[0]-(dim/2), npoint[1]-(dim/2)
            claim_region(mask, top, bottom, left, right)
            new_group.append([point[0]+diff_x, point[1]+diff_y])
        else:
            top, bottom, left, right = control_edge(point, dim, dim)
            claim_region(mask, top, bottom, left, right)
            new_group.append(match)

    return new_group

def test_content(img, point, dim):
    #print img.shape
    top, bottom, left, right = control_edge(point, dim, dim)
    img = img[top:bottom,left:right]
    labmask = LABMaskImage(img)
    hsvmask = HSVMaskImage(img)
    mask = combinemasks(labmask, hsvmask)
    blackpixels = np.where(mask==0)
    blackpixels = np.matrix(blackpixels)
    return blackpixels.shape[1]

if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-i', type="string", nargs=1, dest="input", help="input image")
    parser.add_option('-o', type="string", nargs=1, dest="output", help="output file")
    parser.add_option('-c', type="string", nargs=1, dest="centers", help="centers")
    parser.add_option('-s', type="string", nargs=1, dest="signature", help="signature")
    parser.add_option('-d', type="int", nargs=1, dest="dimension", help="dimension")
    options, args = parser.parse_args()

    if options.input==None:
        raise "No input image given, please give input image: -i image"
    else:
        in_file = options.input

    if options.output==None:
        raise "No output name given, please give output name: -o filename"
    else:
        out_file = options.output

    if options.centers==None:
        raise "No centers given, unable to match center points: -c filename.log"
    else:
        centers_filename = options.centers

    if options.signature==None:
        raise "No signature given, unable to find correct center points: -s filename.sign"
    else:
        sign_filename = options.signature

    if options.dimension==None:
        print "Notice: No dimension given, using defualt size 400: -d int"
        dim = 400
    else:
        dim = options.dimension

    inimg = io.imread(in_file)
    height, width, _ = inimg.shape

    print height, width, dim

    ## Overall mask
    labmask = LABMaskImage(inimg)
    hsvmask = HSVMaskImage(inimg)
    mask = combinemasks(labmask, hsvmask)

    #centers = kmeans(blackpixels.T.astype(float), 10, iter=100, thresh=1e-08)[0].astype(int)
    print "Loading Centres"
    centres = loadcenters(centers_filename, mask)
    print "Grouping centres"
    groups = groupCentres(centres)
    print "Loading signature"
    signature, signature_value = readSignature(sign_filename)
    #centres = orderGroup(centres)
    print centres
    print "Find closest group to signature"
    #centres = closestGroup(groups, signature_value)
    centres =  closestGroupSignature(groups, signature, 250)
    print(centres)
    print "Overlay signature on the group"
    signature = overlaySignature(signature, centres)
    copy_centres = centres
    print "Test and fix fit of signature centers to group"
    centres = testRadius(signature, centres, 200, mask, 500)
    centres = orderGroup(centres)
    i = 10
    c = 0
    print "Cropping images"
    for centre in centres:

        count = test_content(inimg, centre, dim)
        if count<1000:
            i -= 1
            continue

        #crop_square(centre, dim, inimg, i, in_file, out_file)

        try:
            square(copy_centres[c], dim, inimg, 0, 3, color=(55, 255, 55))
            square(copy_centres[c], 40, inimg, 0, 5, color=(55, 255, 55))
        except:
            pass
        square(centre, 10, inimg, 0, 5, color=(55, 55, 255))
        square(centre, dim, inimg, i, 3, color=(55, 55, 255))
        square(signature[str(i)], 200, inimg, 0, 3, color=(255, 55, 55))
        i -= 1
        c += 1

    #square(signature['0'], 50, inimg, 0, 3, color=(255, 55, 55))
    io.imsave(out_file, inimg)
```

Support functions for dealing with the signature

``` python
from sys import argv

def readSignature(filename):
    f = open(filename).read().split(";")
    f, sign_value = f[:-1], float(f[-1].split(":")[-1])
    points = {}
    for point in f:
        if point=="": continue
        name, coord = point.split(":")
        points[name] = [int(i) for i in coord.split(",")]
    return points, sign_value

def writeSignature(filename, points):
    f = open(filename, "w")
    print points.keys()
    order = sorted(points.keys(), key=lambda x: int(x))
    for point in order:
        f.write(point+":"+",".join([str(i) for i in points[point]])+";")
    f.write("Value:"+str((relativeDistanceMeanSign(points)+centreMean(points['0'], points))/2))
    f.close()

def distance(point1, point2):
    return ((point1[0]-point2[0])**2+(point1[1]-point2[1])**2)**0.5

def centreSum(center, points):
    csum = 0
    if type(points)==dict:
        for point in points.values():
            csum += distance(center, point)
    if type(points)==list:
        for point in points:
            csum += distance(center, point)
    return csum

def centreMean(center, points):
    csum = 0
    if type(points)==dict:
        for point in points.values():
            csum += distance(center, point)
    if type(points)==list:
        for point in points:
            csum += distance(center, point)
    return csum/len(points)

def findCentre(points):
    if type(points)==dict:
        x_values = [xy[0] for xy in points.values()]
        y_values = [xy[1] for xy in points.values()]
    if type(points)==list:
        y_values = [xy[0] for xy in points]
        x_values = [xy[1] for xy in points]
    n = len(points)
    return [sum(x_values)/n, sum(y_values)/n]

def relativeDistanceSumSign(points):
    csum = 0
    for i in range(1, 11):
        if (i-1)>=1: csum += distance(points[str(i)], points[str(i-1)])
        if (i+1)<=5: csum += distance(points[str(i)], points[str(i+1)])
        if (i+5)<=10: csum += distance(points[str(i)], points[str(i+5)])
    return csum

def relativeDistanceSum(abovepoints, belowpoints):
    csum = 0
    n = len(abovepoints)
    m = len(belowpoints)
    for i, point in enumerate(abovepoints):
        if (i-1)>=1: csum += distance(point, abovepoints[i-1])
        if (i+1)<n: csum += distance(point, abovepoints[i+1])
        if i<m: csum += distance(point, belowpoints[i])
        else: csum += distance(point, belowpoints[-1])
    for i, point in enumerate(belowpoints):
        if (i-1)>=1: csum += distance(point, belowpoints[i-1])
        if (i+1)<m: csum += distance(point, belowpoints[i+1])
    return csum

def relativeDistanceMeanSign(points):
    csum = 0
    for i in range(1, 11):
        if (i-1)>=1: csum += distance(points[str(i)], points[str(i-1)])
        if (i+1)<=5: csum += distance(points[str(i)], points[str(i+1)])
        if (i+5)<=10: csum += distance(points[str(i)], points[str(i+5)])
    return csum/(len(points)-1)

def relativeDistanceMean(abovepoints, belowpoints):
    csum = 0
    n = len(abovepoints)
    m = len(belowpoints)
    for i, point in enumerate(abovepoints):
        if (i-1)>=1: csum += distance(point, abovepoints[i-1])
        if (i+1)<n: csum += distance(point, abovepoints[i+1])
        if i<m: csum += distance(point, belowpoints[i])
        else: csum += distance(point, belowpoints[-1])
    for i, point in enumerate(belowpoints):
        if (i-1)>=1: csum += distance(point, belowpoints[i-1])
        if (i+1)<m: csum += distance(point, belowpoints[i+1])
    return csum/(n+m)

if __name__=="__main__":
    #writeSignature("signature.sign", p)
    points, sign_value = readSignature(argv[1])
    bestcentre = points['0']
    del points['0']
    print "Signature Value:", sign_value
    print "Centre Sum:", centreSum(bestcentre, points)
    print "Centre Sum:", centreSum([bestcentre[0]+1, bestcentre[0]-1], points)
    print "Centre Sum:", centreSum([bestcentre[0], bestcentre[0]-1], points)
    print "Centre Sum:", centreSum([bestcentre[0]+45, bestcentre[0]-102], points)
```
