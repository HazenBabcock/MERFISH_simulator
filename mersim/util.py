#!/usr/bin/env python
"""
Utility functions.
"""
import cv2
import math
import numpy as np
import shapely
import shapely.geometry
import scipy
import scipy.ndimage
import skimage
import skimage.draw


def all_fov(simParams):
    """
    Return all FOV and their union.
    """
    allFOV = []
    for fov in range(simParams.get_number_positions()):
        allFOV.append(simParams.get_fov_rect(fov))

    allUnion = shapely.ops.unary_union(allFOV)
    b = allUnion.bounds
    return [allFOV, allUnion, [math.ceil(b[2]-b[0]), math.ceil(b[3]-b[1])]]

    
def add_images(im1, im2, xo, yo):
    """
    Add im2 to im1, it is assumed that im2 is smaller than im1.

    https://stackoverflow.com/questions/9886303/adding-different-sized-shaped-displaced-numpy-matrices
    """
    xr1 = slice(max(0, xo), max(min(xo + im2.shape[0], im1.shape[0]), 0))
    yr1 = slice(max(0, yo), max(min(yo + im2.shape[1], im1.shape[1]), 0))
    
    xr2 = slice(max(0, -xo), min(-xo + im1.shape[0], im2.shape[0]))
    yr2 = slice(max(0, -yo), min(-yo + im1.shape[1], im2.shape[1]))

    im1[xr1, yr1] += im2[xr2, yr2]


def concat(list1, list2):
    if list1[0] is None:
        return list2
    else:
        for i in range(len(list1)):
            list1[i] = np.concatenate((list1[i], list2[i]))
        return list1


def convolve(image, kernel):
    return cv2.filter2D(image, -1, kernel,
                        borderType = cv2.BORDER_REPLICATE)


def gradient_fill(polygons, bounds, fovShape, gradient = 0.1):
    """
    Return image with gradient fill. The 'on' values image
    will start at 1.0 and increase by gradient * pixel
    distance from the edge of the fill.
    """
    img = uniform_fill(polygons, bounds, fovShape)
    img_grad_fill = scipy.ndimage.morphology.distance_transform_edt(img)
    img_grad_fill = gradient * img_grad_fill + (1.0 - gradient)

    return img_grad_fill


def random_points_in_shape(poly, density):
    """
    Return XY points randomly distributed inside a polygon.
    """
    minx, miny, maxx, maxy = poly.bounds
    totalPnts = int(poly.area * density)

    pntX = np.zeros(totalPnts)
    pntY = np.zeros(totalPnts)
    cnt = 0
    while (cnt < totalPnts):
        # Choose random XY.
        pnt = shapely.geometry.Point(np.random.uniform(minx, maxx),
                                     np.random.uniform(miny, maxy))

        if poly.contains(pnt):
            pntX[cnt] = pnt.x
            pntY[cnt] = pnt.y
            cnt += 1

    return [pntX, pntY]


    
def uniform_fill(polygons, bounds, fovShape):
    """
    Return binary image of the polygons in poly.
    """
    minx, miny, maxx, maxy = list(map(int, bounds))
    img = np.zeros(fovShape, dtype = np.uint8)
    for poly in polygons:
        r = np.array(poly.exterior.coords.xy[0]) - minx
        c = np.array(poly.exterior.coords.xy[1]) - miny
        [rr, cc] = skimage.draw.polygon(r, c, img.shape)
        img[rr, cc] = 1

    return img

    
def uniform_points_in_shape(poly, spacing, dither = 0.0):
    """
    Return XY points uniformly distributed (on a grid) inside
    a polygon.
    """
    minx, miny, maxx, maxy = poly.bounds
    totalPnts = int(2 * poly.area/(spacing * spacing))

    xv = np.arange(minx, maxx, spacing)
    yv = np.arange(miny, maxy, spacing)
    
    pntX = np.zeros(totalPnts)
    pntY = np.zeros(totalPnts)
    
    cnt = 0
    for xi in range(xv.size):
        for yi in range(yv.size):
            if (dither > 0):
                pnt = shapely.geometry.Point(xv[xi] + np.random.uniform(-dither, dither),
                                             yv[yi] + np.random.uniform(-dither, dither))
            else:
                pnt = shapely.geometry.Point(xv[xi], yv[yi])

            if poly.contains(pnt):
                pntX[cnt] = pnt.x
                pntY[cnt] = pnt.y
                cnt += 1

    return [pntX[:cnt], pntY[:cnt]]
