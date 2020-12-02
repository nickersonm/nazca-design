#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------
# This file is part of Nazca.
#
# Nazca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# Nazca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
#
#@Authors: Katarzyna Lawniczuk, Ronald Broeke
#@email: nazca@brightphotonics.eu
#2017(c)

"""
**DEMO PDK technology settings.**

This module is intended to define the mask layers and xsections of a technology,
which can be accomplished by loading the coresponding csv tables with those
settings and/or use the nazca function that add layers and sections.

The **Xsection objects** are defined in this module. A xsection corresponds typically
to a waveguide or metal line. The following properies of the xsection can be set

    * os: straight-to-bend waveguide offset
    * width: default waveguide width
    * radius: default bend radius
    * taper: taper length
"""



import os
import nazca as nd
from math import copysign

# full path of the (demo) directory that has __init__ file.
dir_path = os.path.dirname(__file__)
gds_path = os.path.join(dir_path, 'gdsBB/')


#==============================================================================
# Load Demo foundry layers from tables store as csv files.
#==============================================================================
layer_file           = os.path.join(dir_path, 'table_layers.csv')
xsection_layer_file  = os.path.join(dir_path, 'table_xsection_layer_map.csv')
xsection_file        = os.path.join(dir_path, 'table_xsections.csv')
layercolor_file      = os.path.join(dir_path, 'table_colors_Light.csv')

#load the foundry tables
nd.load_layers(filename=layer_file)
nd.load_xsections(filename=xsection_file)
nd.load_xsection_layer_map(filename=xsection_layer_file)
nd.load_layercolors(filename=layercolor_file)


nazca_logo_layers = {'ring':1, 'box':None, 'bird':11}
nazca_logo = nd.nazca_logo(nazca_logo_layers)

#==============================================================================
# Define xsection related functions as needed.
#==============================================================================
def os_shallow(width, radius):
    """Offset straight to bend for shallow waveguides.

    Args:
        width (float): waveguide width
        radius (float): waveguide radius

    Returns:
        float: offset value
    """
    if radius == 0:
        offset = 0
    else:
        offset = min(0.1*(200/radius)**2, xsShallow.width*0.4)
    return copysign(offset, radius)

def os_deep(width, radius):
    """Offset straight to bend for deep waveguides.

    Args:
        width (float): waveguide width
        radius (float): waveguide radius

    Returns:
        float: offset value
    """
    if radius == 0:
        offset = 0
    else:
        offset = min(0.05*(30/radius)**2, xsDeep.width*0.4)
    return copysign(offset, radius)


#==============================================================================
# xsection attributes
#==============================================================================
xsShallow = nd.add_xsection('Shallow')
xsShallow.os = os_shallow
xsShallow.width = 3.0
xsShallow.radius = 200.0
xsShallow.minimum_radius = 200.0
xsShallow.taper = 50.0

xsDeep = nd.add_xsection('Deep')
xsDeep.os = os_deep
xsDeep.width = 1.5
xsDeep.radius = 75.0
xsDeep.minimum_radius = 75.0
xsDeep.taper = 50.0

xsMetalDC = nd.add_xsection('MetalDC')
xsMetalDC.os = 0.0
xsMetalDC.width = 30.0
xsMetalDC.radius = 30.0
xsMetalDC.taper = 50.0

xsMetalRF = nd.add_xsection('MetalRF')
xsMetalRF.os = 0.0
xsMetalRF.width = 20.0
xsMetalRF.radius = 20.0
xsMetalRF.taper = 50.0


# =============================================================================
# DRC
# =============================================================================
# No angle DRC for the metal:
nd.get_xsection('MetalDC').drc_angle = False # Allow for all angles: skip the angle drc.
nd.get_xsection('MetalRF').drc_angle = False # Allow for all angles: skip the angle drc.
nd.get_xsection('MetalDC').drc_width = False # Allow for non-equal width connections.
nd.get_xsection('MetalRF').drc_width = False # Allow for non-equal width connections.

# Allow xsections to be connected to others than itself.
nd.cfg.drc_rule_xs = {
    'MetalDC': ['MetalDC', 'MetalRF'],
    'MetalRF': ['MetalDC', 'MetalRF'],
    }




