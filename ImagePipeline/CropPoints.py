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
