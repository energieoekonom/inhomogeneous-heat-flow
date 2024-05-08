#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:05:10 2024

@author: economist
"""

import numpy as np

def heat_flow_magnitude(T, C):
    """
    computes sum of incoming and outgoing heat flows for grid nodes
    omitting the first and last rows, as these had the 1e9 near infinite
    conductivities and don't interest here.

    Parameters
    ----------
    T : array
        temperatures.
    C : array
        heat conductivities.

    Raises
    ------
    ValueError
        for shape mismatch.

    Returns
    -------
    array
        heat flow magnitudes array.

    """
    if T.shape != C.shape:
        raise ValueError(f"shapes T{T.shape} != C{C.shape}")
    
    F = np.zeros(T.shape)
    
    ucomb = lambda u1, u2: 2 * u1 * u2 / (u1 + u2)
    
    # iteration won't touch first and last row of newly allocated
    # array. Easier to index, and afterwards unfilled rows can be removed
    for i in range(1, T.shape[0]-1):
        for j in range(T.shape[1]):
            hf = 0
            Cij = C[i,j]
            Tij = T[i,j]

            def calc_hf(i_,j_):
                u = ucomb(Cij, C[i_,j_])
                dT = Tij - T[i_,j_]
                return abs(dT*u)
            
            if j > 0:
                hf += calc_hf(i, j-1)
            if j < T.shape[1] -1:
                hf += calc_hf(i,j+1)
            
            hf += calc_hf(i-1,j)
            hf += calc_hf(i+1,j)
            F[i,j] = hf
    F = F[1:-2,:]
    return F
                
    
def heat_flow_gradient(T, C):
    """
    computes sum of incoming and outgoing heat flows for grid nodes
    omitting the first and last rows, as these had the 1e9 near infinite
    conductivities and don't interest here.

    Parameters
    ----------
    T : array
        temperatures.
    C : array
        heat conductivities.

    Raises
    ------
    ValueError
        for shape mismatch.

    Returns
    -------
    array
        heat flow magnitudes array.

    """
    if T.shape != C.shape:
        raise ValueError(f"shapes T{T.shape} != C{C.shape}")
    
    X = np.zeros(T.shape)
    Y = np.zeros(T.shape)
    ucomb = lambda u1, u2: 2 * u1 * u2 / (u1 + u2)
    
    # iteration won't touch first and last row of newly allocated
    # array. Easier to index, and afterwards unfilled rows can be removed
    for i in range(1, T.shape[0]-1):
        for j in range(T.shape[1]):
            Cij = C[i,j]
            Tij = T[i,j]
            
            def calc_hf_y():
                u_u = ucomb(Cij, C[i-1,j])
                u_d = ucomb(Cij, C[i+1,j])
                d_u = T[i-1,j] - Tij
                d_d = Tij - T[i+1,j]
                return (u_u * d_u + u_d * d_d)/2
            
            def calc_hf_x():
                if j > 0:
                    u_l = ucomb(Cij, C[i,j-1])
                    d_l = Tij - T[i,j-1]
                else:
                    u_l, d_l = 0.0,0.0
                if j < T.shape[1] -1:
                    u_r = ucomb(Cij, C[i,j+1])
                    d_r = T[i,j+1] - Tij
                else:
                    u_r, d_r = 0.0,0.0
                return (u_l * d_l + u_r * d_r)/2
            
            X[i,j] = calc_hf_x()
            Y[i,j] = calc_hf_y()
    X = X[1:-2,:]
    Y = Y[1:-2,:]
    return X, Y

                
    