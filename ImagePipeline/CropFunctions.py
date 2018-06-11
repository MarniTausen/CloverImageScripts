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
