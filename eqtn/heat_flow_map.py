#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:05:10 2024

@author: energieoekonom
"""

import numpy as np
def heat_flow_magnitude(T, C):
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
                
    
    