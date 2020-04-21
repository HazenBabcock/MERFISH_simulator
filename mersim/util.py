#!/usr/bin/env python
"""
Utility functions.
"""
import numpy as np


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
