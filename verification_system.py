#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 07:02:46 2024

@author: economist
"""

import sys

# import libraries for numerical functions
import numpy as np
import eqtn.elliptic_system as els

def verification_system():
    # create a wall of 10 simulation layers, 2 layers of inside/outside air
    # and 2 layers representing inside/outside environments.
    # thickness 10cm -> one cm**2 square size assumed. Wall consists of
    # An inner and outer layer of solid clay brick with an in between 
    # thermal insulation layer of polyurethane.
    kappa_brick = 0.62 # W/(mK)
    # polyurethane
    kappa_pur = 0.028 # W/(mK)
    
    # on the scaling of factors: grid is assumed 1cm. So 1cm wide and
    # 1cm thick. Just assume 1m in depth. Then all scaling is effected 
    # by computing matching thermal conductivities of imaginary
    # 1cm air layers on the outside and inside insulating the wall from
    # inner room and outside environment.
    simulation_width = 2
    simulation_depth = 14
    C = kappa_brick * np.ones((simulation_depth, simulation_width))
    # now put the air layers adding to insulation
    # external air resistance
    # Rse in m² • K/W ISO 6946:2017
    R_se = 0.04
    # want a 1cm layer in simulation to have R_se
    # hence kappa_se = 1/R_se / 100
    kappa_se = 1/R_se/100
    # internal air resistance
    R_si = 0.13
    # want a 1cm layer in simulation to have R_si
    # hence kappa_se = 1/R_se / 100
    kappa_si = 1/R_si/100
    # and in this simulation variant, the first layer sits between
    # environment and brick wall, full size.
    kappa_env = 1e9

    # put conductivites different from brick by overwriting respective
    # rows of matrix C
    C[0,:] = kappa_env
    C[-1,:] = kappa_env
    C[1,:] = kappa_si
    C[-2,:] = kappa_se
    
    n_pur = 2
    start_pur = int(simulation_depth/2 - n_pur/2)
    C[(start_pur):(start_pur + n_pur),:] = kappa_pur
    
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
    print("Result matrix of temperatures:")
    print(u_square)
            
    # pick the wall inner and outer temperatures
    t_i = u_square[2,0]
    t_o = u_square[-3,0]
    
    # simulation computes temperatures in middle of layers. Hence the above
    # wall temperatures are half a centimeter into the wall on the 
    # inside and outside. To fix this, compute the temperature drop for
    # that half centimeter. As brick layers are several centimeters, can 
    # compute this from temperature drop between brick simulation nodes.
    half_brick_drop = -(u_square[3,0] - u_square[2,0]) / 2
    t_i += half_brick_drop
    t_o -= half_brick_drop
    
    # and the theoretical values from resistances and thickness of layers
    R_wall_clay = (simulation_depth-4-n_pur) / 100 / kappa_brick
    R_wall_pur = n_pur/ 100 / kappa_pur
    R_total = R_se + R_si + R_wall_clay + R_wall_pur
    
    delta_t_air_e_th = (T_i - T_o) * R_se / R_total
    delta_t_air_i_th = (T_i - T_o) * R_si / R_total
    delta_t_air_e = t_o - T_o
    delta_t_air_i = T_i - t_i
    print(f"Simulated/theortical inner air layer temperature drop: "
          f"{delta_t_air_i:.3f}/{delta_t_air_i_th:.3f}")
    print(f"Simulated/theortical outer air layer temperature drop: "
          f"{delta_t_air_e:.3f}/{delta_t_air_e_th:.3f}")
    
    # get the simulated and theoretical temperature drops across the 
    # layer of polyurethane
    # first, indices of brick bordering the polyurethane
    idx_brick_pur_i = start_pur -1
    idx_brick_pur_e = idx_brick_pur_i + n_pur + 1
    t_pur_i = u_square[idx_brick_pur_i,0] - half_brick_drop
    t_pur_e = u_square[idx_brick_pur_e,0] + half_brick_drop
    delta_t_pur_th = (T_i - T_o) * R_wall_pur / R_total
    print(f"simulated/theoretical temperature drop across polyurethane: "
          f"{t_pur_i - t_pur_e:.3f}/{delta_t_pur_th:.3f}")
    
    print(f"theoretical pur layered brick {simulation_depth-4} cm "
          f"coefficient of heat conductivity:\n"
          f"U [W/(m² K)] = {1/R_total:.3f}")

def main(args):
    verification_system()
    return 0

try:
    args = []
    v = main(args)
    sys.exit(v)
except Exception as ex:
    print(f"Exception: {ex}")
    print("terminating exit(1)")
    sys.exit(1)

