#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:42:00 2024

@author: economist
"""

import matplotlib.image as mplimg

def avg_grayscale_matrix(image):
    """
    from input image data read from a PNG file, having 3 color channels,
    generate an output matrix having grayscales as average of channel
    luminosities.

    Parameters
    ----------
    image : array
        a tensor with PNG image data.

    Returns
    -------
    image2d : array
        a 2d-Matrix of luminosities.

    """
    # image dimensions y-axis, x-axis, 3x color 1x ?
    image2d = image[:,:,0] + image[:,:,1] + image[:,:,2]
    image2d = image2d / 3
    return image2d

def read_grayscale_2d(pngfilename):    
    """
    creates a 2d grayscale matrix from input PNG file

    Parameters
    ----------
    pngfilename : str
        path to png file.

    Returns
    -------
    image2d : array
        2d luminosity array of image.

    """
    image = mplimg.imread(pngfilename)
    # the downscaling averaging is now a simple reshape and mean
    # job. However, the image must first be cut to be evenly
    # divisible at both dimentions by the downscale factor.
    image2d = avg_grayscale_matrix(image)    
    return image2d


def matrix_to_pix(image2d, downscale=10):
    """
    downscale image by taking avery downscale stepped pixel.
    avoids introducing "new" materials at boundaries
    but note that png needs be exported w/o aliasing.

    Parameters
    ----------
    pngfilename : array
        image matrix.
    downscale : int, optional
        Factor to downscale input image. The default is 10.

    Returns
    -------
    outimage : array
        downscaled image matrix.

    """

    dso = int(downscale/2)
    image_xclip = image2d[:,dso::downscale]
    image_yclip = image_xclip[dso::downscale,:]
    return image_yclip

    
    