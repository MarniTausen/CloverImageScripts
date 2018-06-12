Clover image analysis supplementary + Clover grass analysis
================
12/06/2018 - 14:30:44

-   [Introduction](#introduction)
-   [Image Analysis](#image-analysis)
    -   [Cropping](#cropping)
    -   [Finding Pots](#finding-pots)
-   [Growth Curves](#growth-curves)
    -   [Time frame correction](#time-frame-correction)
    -   [Heritability](#heritability)

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

Growth Curves
=============

The main first script which loads and processes all of the data. Some scripts are not printed on here, but are located in the folder GrowthCurveAnalysis.

``` r
library(ggplot2)
library(reshape2)
library(dplyr)
library(drc)

## Contains functinos for converting strings into time.
# source("TimeConversion.R")

## Contains functions for converting the Combinations into placements, and vice versa.
## convertToPlacement(x) where x is the name of the combination
## convertToPot(x) ...
## convertToCamera(x) ...
source("CameraConversion.R")

## Supporting functions that make everything easier.
## mergeVectors(x, y) Merge 2 string vectors by a seperator "_" and returning z
## divide_by_harverst(df, combination_table) Takes a data.frame (df) and uses a combination_table which contains the harvest information
## pickone(x) Takes the first element from vector x, assuming all values in this vector are identical.
## predict_points(df) Takes the TAf table and predicts points linearly between points if there is more than 1 day between them.
source("supportfunctions.R")

##############################################
## Compiling Harvest information
##############################################

harvested <- read.csv("nchain_harvested.tab.csv", sep=";")

full_table <- read.table("clover.tab", header=TRUE, sep="\t")
full_table <- full_table[!is.na(full_table$ID),]
full_table <- full_table[!is.na(full_table$Table),]
rownames(full_table) <- full_table$ID
#full_table$ID <- NULL
full_table$Cutting <- NULL
full_table$Harvest <- NULL
full_table$Potting <- NULL
full_table$Status <- NULL
full_table %>% arrange(ID) -> full_table
full_table$Combination <- mergeVectors(full_table$Clover, full_table$Rhizobium)

HarvestData <- read.csv("harvest1.csv", sep=";")
names(HarvestData)[1] <- "barcode"
HarvestData$Combination <- full_table[HarvestData$barcode,"Combination"]
rownames(HarvestData) <- HarvestData$barcode
HarvestData$harvest_date <- as.Date(HarvestData$harvest_date, format="%d/%m/%y")
HarvestData$cutting_date <- as.Date(HarvestData$cutting_date, format="%d/%m/%y")
HarvestData$inoculation_date <- as.Date(HarvestData$inoculation_date, format="%d/%m/%y")
HarvestData$potting_date <- as.Date(HarvestData$potting_date, format="%d/%m/%y")
HarvestData <- HarvestData[as.character(sort(HarvestData$barcode)),]
rownames(HarvestData) <- HarvestData$Combination

#full_table <- full_table[sort(harvested$Pot), ]
full_table <- full_table[sort(HarvestData$barcode), ]

#row.names(harvested) <- harvested$Pot
#harvested$X.2 <- NULL
#harvested <- harvested[as.character(sort(as.numeric(rownames(harvested)))),]

#full_table$Harvest <- as.Date(harvested$date, format="%d/%m/%Y")
full_table$Harvest <- as.Date(HarvestData$harvest_date, format="%d/%m/%y")

combination_table <- full_table
row.names(combination_table) <- combination_table$Combination

##############################################
## Loading the growth data
##############################################

TA <- read.csv("CloverImageDataTA.csv")

all_names <- colnames(TA)[2:ncol(TA)]

## Filter by harvested
TA <- TA[,c("Time", full_table$Combination)]

## Convert into Date.
TA$Time <- as.Date(TA$Time, format="%Y/%m/%d - %H:%M")


## melt from matrix into continuous data.frame
TA <- melt(TA, id.vars = "Time", variable.name="Pot", value.name = "TotalArea")

## Convert pixels into cm^2
max_pixel_count <- 400*400
crop_in_plate_size <- 40*40

TA$TotalArea <- TA$TotalArea/max_pixel_count
TA$TotalArea <- TA$TotalArea*crop_in_plate_size

## Filter early days.
TA %>% filter(Time!="2017-05-23") %>% filter(Time!="2017-05-22") -> TA

## Add Camera ID and Pot ID to the data.frame
TA %>% mutate(Camera = unlist(convertToCamera(Pot)),
              Placement = unlist(convertToPot(Pot))) -> TA

## Remove all rows with missing data.
TA <- TA[!is.na(TA$TotalArea),]

## Create per day summary data
TA %>% group_by(Pot, Time) %>%
    summarise(Std=sd(TotalArea), TotalArea=median(TotalArea), Camera=pickone(Camera),
              Placement=pickone(Placement)) -> TAf

## Convert into data.frame
#TAf <- as.data.frame(TAf)

## Remove missing data if any.
TAf <- TAf[!is.na(TAf$TotalArea),]
TAf <- TAf[!is.na(TAf$Std),]

## Divided by harvest
TAf$Period <- divide_by_harvest(TAf, combination_table)

show <- function(Period){
    length(unique(Period))>1
}

TAf %>% group_by(Pot) %>% filter(show(Period)) -> TAf

##############################################
## Filtering days (Done manually)
##############################################

## Contains the manual filtering list and black list.
## manual_filter(pot, time) Uses the manual_list to filter if that day is in the list.
## blacklist_filter(pot) If pot is in black_list then filter it.
source('manual_list.R')

TAf %>% group_by(Pot, Time) %>% filter(blacklist_filter(Pot)) -> TAblacklist

in_filter <- Vectorize(function(Pot, Manual_list) {
    exists(as.character(Pot), where=Manual_list)
}, vectorize.args = "Pot")

TAf %>% group_by(Pot, Time) %>% filter(!blacklist_filter(Pot)) %>%
    filter(manual_filter(Pot, Time)) -> TAf

## Add predicted results
TAfpred <- predict_points(TAf)

TAfpred %>% group_by(Pot, Time) -> TAf

############################
# Getting the growth rates #
#####################################
## Get growth rates using the medians
#####################################

TAf %>% group_by(Pot, Period) %>%
    summarise(b=drm(TotalArea ~ Time, fct=L.4())$coefficients[1],
              c=drm(TotalArea ~ Time, fct=L.4())$coefficients[2],
              d=drm(TotalArea ~ Time, fct=L.4())$coefficients[3],
              e=drm(TotalArea ~ Time, fct=L.4())$coefficients[4]) -> GrowthRates

GrowthRates %>% arrange(b) -> GrowthRates

GrowthRates <- as.data.frame(GrowthRates)
#print(GrowthRates)
```

Collection of support functions to make the processing easier/faster

``` r
contains <- function(x, m){
    if(m==x) return(TRUE)
    FALSE
}
contains <- Vectorize(contains)

mergeVectors <- function(x, y){
    if(length(x)!=length(y)) stop("Vectors are of unequal length")
    z <- vector(length=length(x))
    for(i in 1:length(x)){
        tmp <- paste(x[i], y[i], sep=".")
        j <- 0
        while(any(z==tmp)){
            j <- j+1
            tmp <- paste(x[i], y[i], as.character(j), sep=".")
        }
        z[i] <- tmp
    }
    z
}

divide_by_harvest <- function(df, combination_table){
    v <- vector(length=nrow(df))
    for(i in 1:nrow(df)){
        if(df$Time[i]<combination_table[df$Pot[i], 'Harvest']){
            v[i] <- "growth1"
        } else {
            v[i] <- "growth2"
        }
    }
    v
}

pickone <- function(x) x[1]


predict_points <- function(df){

    new_df = data.frame()
    last_day = NULL
    new_pot = TRUE
    pot_name = NULL
    period = NULL
    phantom = vector()

    for(row in 1:nrow(df)){

        if(new_pot==TRUE){
            last_day = df$Time[row]
            pot_name = as.character(df$Pot[row])
            period = df$Period[row]
            new_pot = FALSE
            if(nrow(new_df)==0) new_df = df[row,]
            else new_df = rbind(new_df, df[row,])
            phantom <- c(phantom, "true data")
        } else {

            if(as.numeric(df$Time[row]-last_day) > 1){
                Days <- last_day+1:(as.numeric(df$Time[row]-last_day)-1)

                ## Find the results inbetween the Days.
                temp <- data.frame(Time=df$Time[(row-1):row],
                                   TotalArea=df$TotalArea[(row-1):row])
                ### Do linear regression inbetween points
                fit <- lm(TotalArea ~ Time, data=temp)
                ### Predict points between.
                new_points <- predict(fit, newdata=data.frame(Time=Days))
                new_days <- data.frame(Pot=rep(df$Pot[row], times=length(Days)),
                                       Time=Days,
                                       Std=rep(0, times=length(Days)),
                                       TotalArea=new_points,
                                       Camera=rep(df$Camera[row], times=length(Days)),
                                       Placement=rep(df$Placement[row], times=length(Days)),
                                       Period=rep(df$Period[row], times=length(Days)))

                new_days %>% group_by(Pot, Time) -> new_days
                new_df = rbind(new_df, new_days)
                phantom <- c(phantom, rep("predicted data", times=length(Days)))

            }

            last_day = df$Time[row]
            new_df = rbind(new_df, df[row,])
            phantom <- c(phantom, "true data")

        }

        if((row+1)<nrow(df) && pot_name!=as.character(df$Pot[row+1])){
            new_pot = TRUE
        }
        if((row+1)<nrow(df) && period!=df$Period[row+1]){
            new_pot = TRUE
        }
    }

    new_df$Predicted = phantom
    new_df
}
```

Script which generated all of the growth figures.

``` r
#TA %>% group_by(Pot, Time) %>%
#    summarise(Std=sd(TotalArea), TotalArea=median(TotalArea)) -> TAf

#Cams <- unique((TAf %>% group_by(Camera, Placement, Pot))$Pot)
Cams <- unique(GrowthRates$Pot)

library(drc)

for(Cam in Cams){
    #TAf %>% filter(Camera==Cam) -> TAplot
    TAf %>% filter(Pot==Cam) -> TAplot

    print(Cam)

    if(nrow(TAplot)==0) next

    TAplot %>% filter(Period=="growth1") -> TAplot1
    TAplot %>% filter(Period=="growth2") -> TAplot2

    figure <- ggplot(TAplot, aes(x=Time, y=TotalArea, color=Period)) + geom_line(size=1.05) +
        geom_errorbar(data=TAplot, aes(ymin=TotalArea-Std, ymax=TotalArea+Std),
                      alpha=0.6, width=0.1) + geom_point(aes(color=Predicted)) +
        labs(x="", title=paste(pickone(TAplot$Pot), pickone(TAplot$Camera),
                               pickone(TAplot$Placement), collapse=" "),
             y="Total Coverage (cm^2)") +
        scale_x_date(date_labels = "%b %d", date_breaks = "2 days") + theme_bw() +
        scale_y_continuous(limits=c(0, 1600)) +
        scale_color_manual(values=c("red", "blue", "gray", "black"))

    if(nrow(TAplot1)!=0){
        logfit1 <- drm(TotalArea ~ Time, data=TAplot1, fct=L.4())
        figure <- figure + geom_line(data=TAplot1, aes(x=Time, y=predict(logfit1)))
    }
    if(nrow(TAplot2)!=0){
        logfit2 <- drm(TotalArea ~ Time, data=TAplot2, fct=L.4())
        figure <- figure + geom_line(data=TAplot2, aes(x=Time, y=predict(logfit2)))
    }

    print(figure)

}
```

Time frame correction
---------------------

Main scripts containing and doing the time frame correction

``` r
calculate_area <- function(x, y, a, b, baseline=0){
    square <- (abs(min(y, b)-baseline))*(abs(a-x))
    triangle <- ((abs(b-y))*(abs(a-x)))/2
    square+triangle
}

area_under_the_curve <- function(TotalArea, Time){
    base <- min(TotalArea)
    n <- length(TotalArea)
    total_area <- 0
    for(i in 1:(n-1)){
        total_area <- total_area + calculate_area(1, TotalArea[i], 2, TotalArea[i+1], base)
    }
    total_area
}

area_ratio <- function(TotalArea, Time){
    total_area <- area_under_the_curve(TotalArea, Time)
    max_area <- (max(TotalArea)-min(TotalArea))*length(TotalArea)
    total_area/max_area
}

select_time_frame <- function(Pot, TotalArea, Time, func){

    #print(as.character(Pot[1]))

### Figure out whether this is the replicate .1, or has one.
    sequence <- unlist(strsplit(as.character(Pot[1]), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) replicate <- TRUE
    else replicate <- FALSE

### Get name of replicate
    if(replicate) other <- paste(head(sequence, n=n-2), collapse = "")
    else other <- paste(c(sequence, ".", "1"), collapse = "")

    H1 <- max(Time)
    H2 <- max((TAf %>% filter(Pot==other, Period=="growth1"))$Time)
    H1_t <- HarvestData[as.character(Pot)[1], "inoculation_date"]
    H2_t <- HarvestData[other, "inoculation_date"]

    if(H1<=H2) { H <- H1; Ht <- H1_t }
    else { H <- H2; Ht <- H2_t }

    in1 <- HarvestData[as.character(Pot)[1], "inoculation_date"]
    in2 <- HarvestData[other, "inoculation_date"]

    start_date <- as.Date("2017-05-24", format="%Y-%m-%d")

    ### Find common regions, which are comparable.
    ## Find biggest common region, with min(HarvestDate)-max(Inoculationdate)

    if(in1<in2){
        y <- in2 - in1
        p_start_date <- start_date

        #cat(as.character(Pot[1]), Ht-H, y, "\n")

        if((Ht-H)>y){
            #cat("This is true \n")
            x <- H - p_start_date
        } else {
            x <- H - p_start_date - y
        }

        p_end_date <- p_start_date + x
    } else {
        y <- in1 - in2
        p_start_date <- start_date+y

        if((Ht-H)>y){
            #cat("This is also true \n")
            x <- H - start_date
            #x <- x + y
            #print(x)
        } else {
            #cat("Why is this one true? \n")
            x <- H - p_start_date
        }

        p_end_date <- p_start_date + x
    }

    data.frame(Time=Time, TotalArea) %>%
        filter(Time>=p_start_date, Time<=p_end_date) -> df

    Time <- df$Time
    TotalArea <- df$TotalArea

    func(TotalArea, Time)
}

is_or_has_replicate <- function(Name, names){
    sequence <- unlist(strsplit(as.character(Name), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) is_replicate <- TRUE
    else is_replicate <- FALSE
    if(is_replicate){
        non_replicate <- paste(sequence[1:(n-2)], collapse="")
        if(any(names==non_replicate)) return(TRUE)
        else return(FALSE)
    } else {
        non_replicate <- paste(c(sequence, ".1"), collapse="")
        if(any(names==non_replicate)) return(TRUE)
        else return(FALSE)
    }
}

all_names <- as.character(unique(TAf$Pot))
TAf %>% filter(is_or_has_replicate(Pot, all_names)) %>% group_by(Pot, Period) %>%
    summarise(AreaRatio=area_ratio(TotalArea, Time),
              Area=area_under_the_curve(TotalArea, Time),
              incoluation_date=HarvestData[Pot, "inoculation_date"][1],
              harvest_date=HarvestData[Pot, "harvest_date"][1]) %>%
    arrange(Pot) -> TAarea

TAf %>% filter(is_or_has_replicate(Pot, all_names), Period=="growth1") %>%
    group_by(Pot, Period) %>%
    summarise(AreaRatio=select_time_frame(Pot, TotalArea, Time, area_ratio),
              Area=select_time_frame(Pot, TotalArea, Time, area_under_the_curve)) -> TAarea

TAarea$Pot <- as.character(TAarea$Pot)

print(as.data.frame(TAarea %>% filter(Period=="growth1") %>% arrange(Pot)))

#write.table(as.data.frame(TAarea %>% arrange(Pot)),
#            file="AreaMeasurements2.csv", sep=",", row.names=FALSE)





select_time_frame_filter <- function(Pot, TotalArea, Time){

    #print(as.character(Pot[1]))

### Figure out whether this is the replicate .1, or has one.
    sequence <- unlist(strsplit(as.character(Pot[1]), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) replicate <- TRUE
    else replicate <- FALSE

### Get name of replicate
    if(replicate) other <- paste(head(sequence, n=n-2), collapse = "")
    else other <- paste(c(sequence, ".", "1"), collapse = "")

    H1 <- max(Time)
    H2 <- max((TAf %>% filter(Pot==other, Period=="growth1"))$Time)
    H1_t <- HarvestData[as.character(Pot)[1], "harvest_date"]
    H2_t <- HarvestData[other, "harvest_date"]

    if(H1<=H2) { H <- H1; Ht <- H1_t }
    else { H <- H2; Ht <- H2_t }

    in1 <- HarvestData[as.character(Pot)[1], "inoculation_date"]
    in2 <- HarvestData[other, "inoculation_date"]

    start_date <- as.Date("2017-05-24", format="%Y-%m-%d")

    ### Find common regions, which are comparable.
    ## Find biggest common region, with min(HarvestDate)-max(Inoculationdate)

    if(in1<in2){
        y <- in2 - in1
        p_start_date <- start_date

        #cat(as.character(Pot[1]), Ht-H, y, "\n")

        if((Ht-H)>y){
            cat("This is true \n")
            x <- H - p_start_date
        } else {
            x <- H - p_start_date - y
        }

        p_end_date <- p_start_date + x
    } else {
        y <- in1 - in2
        p_start_date <- start_date+y

        if((Ht-H)>y){
            cat("This is also true \n")
            x <- H - start_date
            #x <- x + y
            #print(x)
        } else {
            cat("Why is this one true? \n")
            x <- H - p_start_date
        }

        p_end_date <- p_start_date + x
    }

    v <- ((Time>=p_start_date)+(Time<=p_end_date))>1
    v
}

testGrowthCurveMini <- function(Plant, filter=TRUE){
    Plants <- c(Plant, paste(Plant, ".1", sep=""))
    print(Plants)
    if(filter){
        TAf %>% filter(any(as.character(Pot)==Plants), Period=="growth1") %>%
            group_by(Pot) %>% filter(select_time_frame_filter(Pot, TotalArea, Time)) -> TAtemp
    } else {
        TAf %>% filter(any(as.character(Pot)==Plants), Period=="growth1") -> TAtemp
    }

    cat("Incoluation date of:", Plants[1], "\n")
    print(HarvestData[Plants[1], "inoculation_date"])
    cat("Incoluation date of:", Plants[2], "\n")
    print(HarvestData[Plants[2], "inoculation_date"])
    cat("Harvest date of:", Plants[1], "\n")
    print(HarvestData[Plants[1], "harvest_date"])
    cat("Harvest date of:", Plants[2], "\n")
    print(HarvestData[Plants[2], "harvest_date"])

    print(TAarea[TAarea$Pot==Plants[1],])
    print(TAarea[TAarea$Pot==Plants[2],])

    TAplot <- TAtemp

    figure <- ggplot(TAplot, aes(x=Time, y=TotalArea, color=Pot)) + geom_line(size=1.05) +
        geom_errorbar(data=TAplot, aes(ymin=TotalArea-Std, ymax=TotalArea+Std),
                      alpha=0.6, width=0.1) + geom_point(aes(color=Predicted)) +
        labs(x="", title=paste(pickone(TAplot$Pot), pickone(TAplot$Camera),
                               pickone(TAplot$Placement), collapse=" ")) +
        scale_x_date(date_labels = "%b %d", date_breaks = "2 days") + theme_bw() +
        scale_y_continuous(limits=c(0, 1500)) +
        scale_color_manual(values=c("red", "blue", "gray", "black"))

    print(figure)

}

showtime <- function(pot){
    (TAf %>% filter(Pot==pot))$Time
}
```

Heritability
------------

The script containing all of the heritability estimatations, and prepartion of the data for the heritability estimates. In particular the harvest data is setup in this script.

``` r
library(lme4)

TAarea$weight <- HarvestData[TAarea$Pot, "weight"]
TAarea$growthday <- HarvestData[TAarea$Pot, "growth_day"]

TAarea$Camera <- unlist(convertToCamera(TAarea$Pot))
TAarea$Placement <- unlist(convertToPot(TAarea$Pot))

combination_split <- Vectorize(function(x, i) unlist(strsplit(x, split=".", fixed=TRUE))[i],
                               vectorize.args = "x")

is_replicate <- function(x) c(2, 1)[is.na(combination_split(x, 3))+1]

join <- Vectorize(function(x, y) paste(x, y, sep="."))

TAarea$Clover <- combination_split(TAarea$Pot, 1)
TAarea$Rhizobium <- combination_split(TAarea$Pot, 2)
TAarea$Combination <- join(TAarea$Clover, TAarea$Rhizobium)
TAarea$Replicate <- is_replicate(TAarea$Pot)

TAarea$Inoculationdate <- HarvestData[TAarea$Pot, "inoculation_date"]
TAarea$Harvestdate <- HarvestData[TAarea$Pot, "harvest_date"]


TAarea$Period <- NULL

Get_originals <- function(Name){
    sequence <- unlist(strsplit(as.character(Name), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) return(FALSE)
    else return(TRUE)
}

Get_originals <- Vectorize(Get_originals)

Get_replicate <- function(Name){
    paste(Name, ".1", sep="")
}

any_match <- function(x, y){
    v <- vector(length=length(x))
    j <-
    for(i in 1:length(x)){
        if(any(x[i]==y)) v[i] <- TRUE
        else v[i] <- FALSE
    }
    v
}

## Score of 1 perfect consitency, 1-0 (varying degrees of constistency)
## <0 inverse consistency, meaning the replicates are worse than the average.
replicate_consistency <- function(dataset, measure){
    x <- dataset[[measure]]
    x <- x/max(x)
    originals <- TAarea$Pot[Get_originals(as.character(dataset$Pot))]
    within_variance <- 0
    for(original in originals){
        replicate <- Get_replicate(as.character(original))
        pair <- c(original, replicate)
        new_variance <- var(x[any_match(dataset$Pot, pair)])
        if(!is.na(new_variance)) within_variance <- within_variance+new_variance
    }
    within_variance <- within_variance/length(originals)
    total_variance <- var(x)
    cat("Replicate variance:", within_variance, "\n")
    cat("Total variance:", total_variance, "\n")
    (total_variance-within_variance)/total_variance
}

replicate_consistency(TAarea, "weight")
replicate_consistency(TAarea, "AreaRatio")
replicate_consistency(TAarea, "growthday")
replicate_consistency(TAarea, "Area")

randint <- function(min, max) round(runif(1, min, max))


find_minimum_time_frame <- function(Harvestdates, Inoculationdates){
    as.numeric(min(Harvestdates-Inoculationdates))
}


normalize_weight <- function(growthrates, timeframe){
    growthrates*timeframe
}

# Just transforms it into the growthday measurement, scaled by time.
t <- find_minimum_time_frame(TAarea$Harvestdate, TAarea$Inoculationdate)
TAarea$normweight <- normalize_weight(TAarea$growthday, t)
replicate_consistency(TAarea, "normweight")

normalized_replicate_time <- (TAf %>% filter(is_or_has_replicate(Pot, all_names), Period=="growth1") %>% group_by(Pot) %>% summarise(Days=sum(select_time_frame_filter(Pot, TotalArea, Time))))

normalized_replicate_time <- as.data.frame(normalized_replicate_time)

rownames(normalized_replicate_time) <- as.character(normalized_replicate_time$Pot)

as.data.frame(normalized_replicate_time %>% arrange(as.character(Pot)))

t <- as.numeric(TAarea$Harvestdate-TAarea$Inoculationdate)
TAarea$normweight <- TAarea$weight*(normalized_replicate_time[TAarea$Pot, "Days"]/t)

replicate_consistency(TAarea, "weight")
replicate_consistency(TAarea, "AreaRatio")
replicate_consistency(TAarea, "growthday")
replicate_consistency(TAarea, "Area")
replicate_consistency(TAarea, "normweight")

multi_effective_strains <- as.character(read.delim("multi_effective_strains.csv", header=FALSE)$V1)
multi_effective_strains[1] <- unlist(strsplit(multi_effective_strains[1], split="\277"))[2]

TAarea %>% filter(any_match(Rhizobium, multi_effective_strains)) %>%
    filter(is_or_has_replicate(Pot, all_names)) -> MEarea

replicate_consistency(MEarea, "weight")
replicate_consistency(MEarea, "AreaRatio")
replicate_consistency(MEarea, "growthday")
replicate_consistency(MEarea, "Area")
replicate_consistency(MEarea, "normweight")

fit <- glm(Area ~ Rhizobium, data=MEarea)
fit2 <- glm(growthday ~ Rhizobium, data=MEarea)

qplot(summary(fit2)$coefficients[,1], summary(fit)$coefficients[,1], xlim=c(-0.1, 0.3))

HarvestData$rhizobium <- as.character(HarvestData$rhizobium)

HarvestData %>% filter(any_match(rhizobium, multi_effective_strains)) -> MEharvest


summary(glm(growth_day ~ rhizobium, data=MEharvest))


fit <- lm(Area ~ Rhizobium, data=MEarea)
fit2 <- lm(growth_day ~ rhizobium, data=MEharvest)

summary(fit)
summary(fit2)

qplot(summary(fit2)$coefficients[c(-1, -21),1], summary(fit)$coefficients[-1,1], xlim=c(-0.1, 0.3))


fit <- lm(Area ~ Rhizobium + Clover, data=MEarea)
fit2 <- lm(growth_day ~ rhizobium + clover, data=MEharvest)

summary(fit)
summary(fit2)

qplot(summary(fit2)$coefficients[c(2:20, 22:23),1], summary(fit)$coefficients[2:22,1])

AreaTable <- as.data.frame(summary(fit)$coefficients[2:22,])
AreaTable$Names <- rownames(summary(fit)$coefficients[2:22,])
AreaTable <- AreaTable[, c('Names', 'Estimate')]

GrowthTable <- as.data.frame(summary(fit2)$coefficients[2:22,])
GrowthTable$Names <- rownames(summary(fit2)$coefficients[2:22,])
GrowthTable <- GrowthTable[, c('Names', 'Estimate')]

MEharvest %>% filter(any_match(Combination, MEarea$Pot)) -> MEharvestsample

MEarea$Rhizobium <- as.factor(MEarea$Rhizobium)
MEarea$Clover <- as.factor(MEarea$Clover)

MEharvest$rhizobium <- as.factor(MEharvest$rhizobium)
MEharvest$clover <- as.factor(MEharvest$clover)

fit <- glm(Area ~ Rhizobium + Clover, data=MEarea)
fit2 <- glm(growth_day ~ rhizobium + clover, data=MEharvest)

summary(fit)
summary(fit2)

summary(glht(fit2, linfct=mcp(rhizobium="Tukey")))

qplot(summary(fit2)$coefficients[c(2:20, 22:23),1], summary(fit)$coefficients[2:22,1])

AreaTable <- as.data.frame(summary(fit)$coefficients[2:22,])
AreaTable$Names <- rownames(summary(fit)$coefficients[2:22,])
AreaTable <- AreaTable[, c('Names', 'Estimate')]

GrowthTable <- as.data.frame(summary(fit2)$coefficients[2:22,])
GrowthTable$Names <- rownames(summary(fit2)$coefficients[2:22,])
GrowthTable <- GrowthTable[, c('Names', 'Estimate')]


qplot(MEarea$growthday, MEarea$Area)

cor(MEarea$growthday, MEarea$Area, method="spearman")

cor(summary(fit2)$coefficients[c(2:20, 22:23),1], summary(fit)$coefficients[2:22,1], method="spearman")


fit <- lm(Area ~ Rhizobium, data=MEarea)

summary(fit)

summary(glht(fit, linfct=mcp(Rhizobium="Tukey")))


MEarea$Clover <- factor(MEarea$Clover)

fit <- lm(growthday ~ Rhizobium + Clover, data=MEarea)

summary(fit)

filtered <- (MEarea %>% filter(!any_match(Clover, c("Aalon_0718", "Aoost_0215",
                                                       "Banna_0733", "Clfin_0213",
                                                       "Banca_0947", "Clfin_0102",
                                                       "Aalon_0617", "Aalon_0512")),
                                  !any_match(Rhizobium, c("SM88"))))

fit <- lm(Area ~ Rhizobium + Clover,
          data=(MEarea %>% filter(!any_match(Clover, c("Aalon_0718", "Aoost_0215",
                                                       "Banna_0733", "Clfin_0213",
                                                       "Banca_0947", "Clfin_0102",
                                                       "Aalon_0617", "Aalon_0512")),
                                  !any_match(Rhizobium, c("SM88")))))


summary(fit)

significance_column <- Vectorize(function(pvalue){
    if(pvalue>0.1) return(' ')
    if(pvalue<=0.1 && pvalue>0.5) return('.')
    if(pvalue<=0.5 && pvalue>0.01) return('*')
    if(pvalue<=0.01 && pvalue>0.001) return('**')
    if(pvalue<=0.001) return('***')
})

rhivrhi <- summary(glht(fit, linfct=mcp(Rhizobium="Tukey")))

rvrtable <- data.frame(estimates=rhivrhi$test$coefficients)
rvrtable$std.error <- rhivrhi$test$sigma
rvrtable$tvalue <- rhivrhi$test$tstat
rvrtable$pvalue <- rhivrhi$test$pvalues
rvrtable[[' ']] <- significance_column(rvrtable$pvalue)
rvrtable <- cbind(rownames(rvrtable), rvrtable)
colnames(rvrtable)[1] <- "Test"

rvrtable %>% filter(pvalue<0.05)



arearc <- lm(Area ~ Rhizobium + Clover, data=MEarea)
areatable <- as.data.frame(summary(arearc)$coefficients)
colnames(areatable)[4] <- "p.value"
areatable <- cbind(rownames(areatable), areatable)
colnames(areatable)[1] <- "Genotypes"


areatable %>% filter(grepl("Rhizobium", Genotypes)) %>% filter(p.value<0.05) %>%
    arrange(-Estimate) %>% mutate(' '=significance_column(p.value))

gpdrc <- lm(growth_day ~ rhizobium + clover, data=MEharvest)
gpdtable <- as.data.frame(summary(gpdrc)$coefficients)
colnames(gpdtable)[4] <- "p.value"
gpdtable <- cbind(rownames(gpdtable), gpdtable)
colnames(gpdtable)[1] <- "Genotypes"

gpdtable %>% filter(grepl("rhizobium", Genotypes)) %>% filter(p.value<0.05) %>%
    arrange(-Estimate) %>% mutate(' '=significance_column(p.value))



## Try Bayz for hertitability

source("http://www.bayz.biz/Rbayz.R")

BAYZHOME = "/usr/local/bin/"

fit <- bayz.mm(data=MEarea, resp="Area", fixmod="fac.Replicate", ranmod="fac.Rhizobium+fac.Clover", chain=c(99900, 1000, 50))

bayz.summ(fit)

fit$post.mean[181]/(fit$post.mean[181]+fit$post.mean[117])

## Replicate versus replicates

(MEarea %>% filter(Replicate==1) %>% arrange(Pot))$Area
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$Area

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$Area,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$Area) + theme_classic() +
    labs(x="Area (Replicate 1)", y="Area (Replicate 2)")

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$growthday,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$growthday) + theme_classic() +
    labs(x="Growth per day (Replicate 1)", y="Growth per day (Replicate 2)")

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$weight,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$weight) + theme_classic() +
    labs(x="Weight (Replicate 1)", y="Weight (Replicate 2)")

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$AreaRatio,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$AreaRatio) + theme_classic() +
    labs(x="Area ratio (Replicate 1)", y="Area ratio (Replicate 2)")

cat("Area replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$Area,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$Area), "\n")

cat("Growth per day replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$growthday,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$growthday), "\n")

cat("Weight replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$weight,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$weight), "\n")

cat("Weight replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$AreaRatio,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$AreaRatio), "\n")

## Calculate heritiabilities


calculate_heritability <- function(lmerMod, variable){
    RE <- as.data.frame(VarCorr(lmerMod))
    if(!class(variable)=="character") stop("Variable must be a string of Characters!")
    if(length(variable)==1) {
        return(RE[RE$grp==variable, "vcov"]/sum(RE$vcov))
    }

print(RE)

    any_match <- function(x, y){
        v <- vector(length=length(x))
        j <-
            for(i in 1:length(x)){
                if(any(x[i]==y)) v[i] <- TRUE
                else v[i] <- FALSE
            }
        v
    }

    if(length(variable)>1){
        return(sum(RE[any_match(RE$grp,variable), "vcov"])/sum(RE$vcov))
    }
}

fit <- lmer(growthday ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))

fit <- lmer(Area ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))

fit <- lmer(weight ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))

fit <- lmer(AreaRatio ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))
```
