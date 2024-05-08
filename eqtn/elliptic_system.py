#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 10:56:07 2024

@author: economist
"""

import numpy as np


def fill_system_weights(C):
    """
    Fills a 2D simulation system matrix from heat conductivity map
    with computations of conductivities between squares of the map.
    
    So the conductivities of the two halves of adjacent simulation squares
    are combined into a connection conductivity.
    
    u12 = (u1 * u2) / (u1 + u2)
    
    Parameters
    ----------
    C : array
        a matrix 1. dimension y 2. x conductivity of material.

    Returns
    -------
    A : array
        Square matrix with system weights (x*y)x(x*y) with x,y dim of C.

    """
    shape_A = C.size
    (cy,cx) = C.shape
    A = np.zeros((shape_A, shape_A))
    di = 0
    
    ucomb = lambda u1, u2: 2 * u1 * u2 / (u1 + u2)
    
    for i in range(cy):
        for j in range(cx):
            cij = C[i,j]
            if j > 0:
                A[di,di-1] = ucomb(cij, C[i,j-1])
            if j < cx - 1:
                A[di,di+1] = ucomb(cij, C[i,j+1])
            if i > 0:
                A[di,di-cx] = ucomb(cij, C[i-1,j])
            if i < cy - 1:
                A[di,di+cx] = ucomb(cij, C[i+1,j])
            cw = np.sum(A[di,:])
            A[di,di] = -cw
            di += 1  
    return A


def set_identity(A, size2):
    """
    sets a portion of the equation system array to identity
    for edge condition

    Parameters
    ----------
    A : array
        equation system coefficients array.
    size2 : int
        portion to set to identity.

    Returns
    -------
    None.

    """
    idm = np.identity(A.shape[0])
    A[0:size2,:] = idm[0:size2,:]
    A[-size2:,:] = idm[-size2:,:]
    
def create_system_results(shape_A, cx, v0, v1):
    """
    creates the equation system results vector

    Parameters
    ----------
    shape_A : int
        edge length of system coefficients array.
    cx : x dimension of heat conductivity matrix
        DESCRIPTION.
    v0 : float
        edge value 1.
    v1 : float
        edge value 2.

    Returns
    -------
    b : array
        result vector of equation system.

    """
    b = np.zeros(shape_A)
    # first row temperatures
    frt = v0 * np.ones(cx)
    b[0:cx] = frt
    # last row temperatures
    lrt = v1 * np.ones(cx)
    b[-cx:] = lrt
    return b
    