#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:27:06 2024

@author: economist
"""

import pandas as pd
import numpy as np

class ModelMaterial(object):
    """
    Class instance represents a material mapped from the 
    heat bridge cross-section.
    """
    col_pixel_value = 'Pixel_Value'
    col_u_value = 'U_Value'
    col_fix = 'Fix'
    col_description = 'Description'
    
    # magic pixel luminosity to map the overall thickness of
    # represented structure
    magic_wallThickness = 1000
    
    def __init__(self, pixelValue, materialData):
        self.pixelValue = pixelValue
        self.materialData = materialData

    def isFix(self):
        return self.materialData[self.col_fix] == 1
    
    def getDescription(self):
        return self.materialData[self.col_description]
    
    def getUValue(self):
        return self.materialData[self.col_u_value]
    
    def getPixelValue(self):
        return self.pixelValue

def generate_value_map(values_to_map, map_csv):
    """
    takes a list of grayscale values from a cross-section image,
    and maps these to heat conductivities listed in a mapping CSV file.

    Parameters
    ----------
    values_to_map : list
        list of grayscale steps from a cross section image.
    map_csv : str
        path to a CSV file with mapping information.

    Returns
    -------
    vm : dict
        a dictionary mapping input values to heat conductivities.

    """
    csv = pd.read_csv(map_csv)
    values = csv[ModelMaterial.col_pixel_value]
    vm = {}
    for value in values_to_map:
        diffv = values - value
        mi = np.argmin(np.abs(diffv))
        vm[value] = ModelMaterial(value, csv.loc[mi])
    return vm
        
def get_wall_thickness(map_csv):
    """
    Magic value extractor to read the overall wall thickness from 
    heat conductivity map.

    Parameters
    ----------
    map_csv : str
        path to a CSV file with mapping information.

    Returns
    -------
    ModelMaterial
        a material instance having the magic key mapped thickness.

    """
    csv = pd.read_csv(map_csv)
    values = csv[ModelMaterial.col_pixel_value]
    diffv = values - ModelMaterial.magic_wallThickness
    mi = np.argmin(np.abs(diffv))
    vr = csv.loc[mi]
    return ModelMaterial(vr[ModelMaterial.col_u_value],csv.loc[mi])    
