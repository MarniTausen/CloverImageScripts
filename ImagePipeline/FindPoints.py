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
