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
