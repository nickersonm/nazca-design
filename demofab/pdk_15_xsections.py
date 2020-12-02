#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**DEMO PDK xsectional structures and epi.**

This module is intended to define the xsectional structure of the waveguides of a technology.
The structure can be defined using horizontal and vertical layers stacks of the materials defined in
pdk_10_materials.py.
"""


import nazca as nd
from . import pdk_10_materials as mat

#Introduce xsection objects. Provide a string name for future referencing.
xsShallow = nd.add_xsection('Shallow')
xsDeep = nd.add_xsection('Deep')

#==============================================================================
# 2D waveguide xsections
# Define waveguide xsections via the Xsection methods.
# A xsection description can be used in visualisation or input to mode solvers.
#==============================================================================
wguide = 3.0
hfilm = 0.6
hsub = 1.0
hclad = 1.0
epi = [(mat.InP,hsub), (mat.Q125,hfilm), (mat.InP,hclad)]

#Shallow guide
xsShallow.layers = epi
xsShallow.background = mat.air
Ls1 = xsShallow.add_vstack(name='backgrnd', etchdepth=hclad+0.2)
Ls2 = xsShallow.add_vstack(name='ridge')
xsShallow.add_hstack(name='guide', layers=[(Ls1,1.0), (Ls2,wguide), (Ls1,1.0)])

#Deep guide
xsDeep.layers = epi
xsDeep.background = mat.air
Ld1 = xsDeep.add_vstack(name='backgrnd', etchdepth=hclad+1.0)
Ld2 = xsDeep.add_vstack(name='ridge')
xsDeep.add_hstack(name='guide', layers=[(Ld1,1.0), (Ld2,wguide), (Ld1,1.0)])
