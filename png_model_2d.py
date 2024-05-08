#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:13:28 2024

@author: economist

Reads png image representing a cross section of material to model
accompanied by a CSV file having grayscale to heat conductivity mappings.

Sample command line:
    png_model_2d.py -pixel_step 20 cross_sections/cross-section2.png

"""

import sys
import argparse
from pathlib import Path

import img.png_to_matrix as pngm
import img.value_map as pnvm

# import libraries for numerical functions and plotting
import numpy as np
import matplotlib.pyplot as plt
import eqtn.elliptic_system as els
import eqtn.heat_flow_map as hfm

# these lines are only for helping improve the display
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'png')
plt.rcParams['figure.dpi']= 300
plt.rcParams['savefig.dpi'] = 300

argv = sys.argv[1:]

parser = argparse.ArgumentParser()
parser.add_argument('-pixel_step', type=int, default=10)
parser.add_argument('pngimage', type=str)

args = parser.parse_args(argv)

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

def print_pixel_portions(v,c,vm):
    """
    prints a breakdown of materials mapped in the input image

    Parameters
    ----------
    v : list
        list of luminosity values found in grayscale image.
    c : list
        list of counts of pixels having a luminosity.
    vm : dict
        value map from luminosity to material description.

    Returns
    -------
    None.

    """
    percent = c/sum(c)*100
    for zv,zc in zip(v,percent):
        model = vm[zv]
        fix = ''
        if model.isFix():
            fix = ' fix'
        print(f"{zv:.3f} -> {model.getUValue()}, {model.getDescription()}{fix}: {zc:.2f}%")

def print_wall_portions(v,c,vm, thickness):
    """
    similar to above print_pixel_portions, prints out wall constituents
    relative to their surface in the cross-section, but taking out
    the inside and environment air layers.

    Parameters
    ----------
    v : list
        list of luminosity values found in grayscale image.
    c : list
        list of counts of pixels having a luminosity.
    vm : dict
        value map from luminosity to material description.
    thickness : ModelMaterial
        magic value mapped ModelMaterial having the overall thickness.

    Returns
    -------
    None.

    """
    # filter the fixed layers (the in and outside air resistance)
    the_wall = [(zv,zc) for zv,zc in zip(v,c) if not vm[zv].isFix()]
    the_air = [(zv,zc) for zv,zc in zip(v,c) if vm[zv].isFix()]
    v,c = zip(*the_wall)
    va,ca = zip(*the_air)
    percent = c/sum(c)*100
    for zv,zc in zip(v,percent):
        model = vm[zv]
        print(f"{zv:.3f} -> {model.getUValue()}, {model.getDescription()}: {zc:.2f}%")
    
    # Now a naive modelling breakdown
    # 1) The constituents of the wall could either be side by side so that their
    # conductivites simply add up to the wall's conductivity. A shortcut to that is
    # to assume that the wall has the consituent's average conductivity.
    # 2) Alternatively, the wall could be built layer by layer, so that the resistances
    # of all layers add up.
    vv = np.array([vm[v_].getUValue() for v_ in v])
    cv = np.array(c)
    va = [vm[v_].getUValue() for v_ in va]
    U_avg = sum(vv * cv) / sum(cv) / thickness.getPixelValue()
    Us_case_1 = va + [U_avg]
    R_case_1 = sum([1/u_ for u_ in Us_case_1])
    U_case_1 = 1/R_case_1
    print(f"case 1 mixed compound U = {U_case_1:.2f}")
    # for 2), first get the thickness of layers
    print("")
    print("case 2 layer thicknesses")
    thicknesses = c / sum(c) * thickness.getPixelValue()
    for v_,t_ in zip(v,thicknesses):
        m = vm[v_]
        print(f"{m.getDescription()}: thickness {t_*100:.2f} cm")
    Us_case_2 = [vm[v_].getUValue()/t_ for v_,t_ in zip(v,thicknesses)]
    Us_case_2 = Us_case_2 + va
    R_case_2 = sum([1/u_ for u_ in Us_case_2])
    U_case_2 = 1/R_case_2
    print(f"case 2 layered wall U = {U_case_2:.2f}")

def png_model2d(args):
    """
    runs the 2d finite difference model with parameters from command line

    Parameters
    ----------
    args : namespace
        parameters passed from command line.

    Raises
    ------
    ValueError
        unspecified error condition.

    Returns
    -------
    None.

    """
    grayimage2d = pngm.read_grayscale_2d(args.pngimage)
    pixelpicked2d = pngm.matrix_to_pix(grayimage2d, args.pixel_step)
    
    v_grayimage2d, c_grayimage2d = np.unique(grayimage2d, return_counts=True)
    v_pixelpicked2d, c_pixelpicked2d = np.unique(pixelpicked2d, return_counts=True)

    if len(v_grayimage2d) != len(v_pixelpicked2d):
        raise ValueError(f"error in image or pixel picking."
                         f" Have {len(v_grayimage2d)} values in original,"
                         f" {len(v_pixelpicked2d)} in reduced x {args.pixel_step}")
    
    path = Path(args.pngimage)
    imgBaseName = path.stem
    
    wallThickness = pnvm.get_wall_thickness(f"{args.pngimage}.csv")
    
    print(f"Wall thickness [m]: {wallThickness.getPixelValue()}")
    print("")    
    value_map = pnvm.generate_value_map(v_grayimage2d, f"{args.pngimage}.csv")

        
        
        
        

    print("grayimage wall portions")
    print_pixel_portions(v_grayimage2d, c_grayimage2d, value_map)
    print("")
    print("pixelpicked downscale in comparison")
    print_pixel_portions(v_pixelpicked2d, c_pixelpicked2d, value_map)
    
    print("")
    print("taking outer air layers out of the balance gives wall volume components")
    print_wall_portions(v_grayimage2d, c_grayimage2d, value_map, wallThickness)

    ###########################################################
    # now the fun bit: the numerical approximation :)
    ###########################################################
    
    # first, reduce the image size by giving the air layers a single row
    # in the image matrix - since the air layers have fixed U value that's
    # acceptable so long as the conductivity is adapted to the equivalent
    # layer's thickness
    p2d = pixelpicked2d # handier alias please!
    voc = p2d[0,0]
    vlc = p2d[-1,-1]
    # now vertical vectors of having this value
    vvo = np.arange(p2d.shape[0])[p2d[:,0] == voc]
    vvl = np.arange(p2d.shape[0])[p2d[:,0] == vlc]
    # the area of interest is now 2 rows above and below the plaster layer
    # of the wall.
    p2d = p2d[(max(vvo)-1):(min(vvl)+2),:]
    # as a modelling hack, assume an environment of infinite conductivity
    # separated from the wall by a 1 pixel wide air insulation layer.
    U_env = 1e9
    p2d[0,:] = U_env
    p2d[-1,:] = U_env
    # now replace the luminosities from the image by their mapped
    # heat conductivities.
    p2dm = np.array(p2d)
    for v in value_map.keys():
        mask = p2d == v
        p2dm[mask] = value_map[v].getUValue()
    # scale the heat conductivities in the environment air layers to 
    # the edge length of one pixel square
    pixels_wall = p2dm.shape[0] - 4
    pixel_edge_length = wallThickness.getPixelValue() / pixels_wall
    p2dm[1,:] = p2dm[1,:] * pixel_edge_length
    p2dm[-2,:] = p2dm[-2,:] * pixel_edge_length
    
    # now set up the equation system
    shape_A = p2dm.size
    (cy,cx) = p2dm.shape

    A = els.fill_system_weights(p2dm)
    els.set_identity(A, cx)
    
    # boundary conditions - per default equations resolve to zero
    # but can set temperature of chosen nodes
    # first row temperatures
    T_i, T_o = 20.0, -15.0
    b = els.create_system_results(shape_A, cx, T_i, T_o)
    
    u = np.linalg.solve(A, b)
    u_square = u.reshape((cy,cx))
    
    u_square_flip = flip_vertical(u_square)
    
    plt.contourf(u_square_flip, levels=40)
    plt.colorbar()
    plt.savefig(f"{imgBaseName}.temp_contour.png")
    plt.show()
    
    # compute the inner and outer wall surface temperatures
    delta_t_inner_air = (u_square[0,:] - u_square[1,:]) * 2
    delta_t_outer_air = (u_square[-2,:] - u_square[-1,:]) * 2
    T_inner_wall = u_square[0,:] - delta_t_inner_air
    T_outer_wall = u_square[-1,:] + delta_t_outer_air
    
    plt.plot(np.arange(len(T_inner_wall)), T_inner_wall)
    plt.savefig(f"{imgBaseName}.temp_inner_wall.png")
    plt.show()
    
    plt.plot(np.arange(len(T_outer_wall)), T_outer_wall)
    plt.savefig(f"{imgBaseName}.temp_outer_wall.png")
    plt.show()
    
    # temperature drops in environment air layers let us compute
    # the wall's U value from the cross sectional heat flow.
    delta_t_inner_air_avg = np.mean(delta_t_inner_air)
    R_inner_air = 1 / p2dm[1,0] * pixel_edge_length
    R_total = (T_i - T_o)/delta_t_inner_air_avg * R_inner_air
    U_total = 1/R_total
    print(f"simulated inner wall U = {U_total}")
    
    delta_t_outer_air_avg = np.mean(delta_t_outer_air)
    R_outer_air = 1 / p2dm[-2,0] * pixel_edge_length
    R_o_total = (T_i - T_o)/delta_t_outer_air_avg * R_outer_air
    U_o_total = 1/R_o_total
    print(f"simulated outer wall U = {U_o_total}")
    
    # heat flow gradient must be run with rougher mesh to obtain
    # meaningful arrows
    gX, gY = hfm.heat_flow_gradient(u_square, p2dm)
    gx, gy = np.meshgrid(np.arange(gX.shape[1]), np.arange(gX.shape[0]))
    gX = -flip_vertical(gX)
    gY = -flip_vertical(gY)
    plt.figure()
    plt.quiver(gx,gy, gX,gY)
    plt.savefig(f"{imgBaseName}.heat_flow_gradient.png")
    plt.show()
    
    # heat flow magnitude looks good in fine meshes
    heat_flow = hfm.heat_flow_magnitude(u_square, p2dm)
    heat_flow_flip = flip_vertical(heat_flow)
    plt.contourf(heat_flow_flip, levels=40)
    plt.colorbar()
    plt.savefig(f"{imgBaseName}.heat_flow.png")
    plt.show()    
    
    



def main(args):
    png_model2d(args)
    return 0

try:
    v = main(args)
    sys.exit(v)
except Exception as ex:
    print(f"Exception: {ex}")
    print("terminating exit(1)")
    sys.exit(1)

