#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 07:02:46 2024

@author: economist
"""
import sys

# import libraries for numerical functions and plotting
import numpy as np
import matplotlib.pyplot as plt
import eqtn.elliptic_system as els

# these lines are only for helping improve the display
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'png')
plt.rcParams['figure.dpi']= 300
plt.rcParams['savefig.dpi'] = 300

def flip_vertical(A):
    """
    flips matrix along axis 0

    Parameters
    ----------
    A : array
        input matrix.

    Returns
    -------
    A_flip : array
        flipped matrix.

    """
    cy = A.shape[0]
    A_flip = A[np.arange(cy)[::-1],:]
    return A_flip

def minimal_example():
    C = np.array([
                            [1.0, 2.0, 3.0],
                            [4.0, 5.0, 6.0],
                            [7.0, 8.0, 9.0]
                            ])
    shape_A = C.size
    (cy,cx) = C.shape
    A = els.fill_system_weights(C)
    els.set_identity(A, cx)
    
    # boundary conditions - per default equations resolve to zero
    # but can set temperature of chosen nodes
    # first row temperatures
    T_i, T_o = 30.0, 20.0
    b = els.create_system_results(shape_A, cx, T_i, T_o)

    u = np.linalg.solve(A, b)

    u_square = u.reshape((cy,cx))
    
    # align image axes with matrix
    u_square_flip = flip_vertical(u_square)
    
    plt.contourf(u_square_flip, levels=20)
    plt.colorbar()
    plt.show()
    
    

def main(args):
    minimal_example()
    return 0

try:
    args = []
    v = main(args)
    sys.exit(v)
except Exception as ex:
    print(f"Exception: {ex}")
    print("terminating exit(1)")
    sys.exit(1)

