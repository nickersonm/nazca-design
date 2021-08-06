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
"""
**DEMO PDK, EPI and xsectional waveguide structures.**

This module is intended to define the xsectional structure of the waveguides in a technology.
The structure can be defined using horizontal and vertical layers stacks of the materials defined in
pdk_10_materials.py.

Note that for layout purposes this module is not required.
"""


import nazca as nd
from . import pdk_10_materials as mat10

# Create xsection objects. Provide a string name for future referencing.
xsShallow = nd.add_xsection('Shallow')
xsDeep = nd.add_xsection('Deep')


#==============================================================================
# 2D waveguide xsections
# Define waveguide xsections via the Xsection methods.
# A xsection description can be used in visualisation or input to mode solvers.
#==============================================================================
wguide = 3.0  # waveguide width in um
hfilm = 0.6  # thickness of the guiding layer in um
hsub = 1.0  # (optical) substrate thickness in um
hclad = 1.0  # cladding thickness in um
epi = [(mat10.InP, hsub), (mat10.Q125, hfilm), (mat10.InP, hclad)]


# Define the Shallow guide cross sectional geometry
xsShallow.layers = epi
xsShallow.background = mat10.air
Ls1 = xsShallow.add_vstack(name='backgrnd', etchdepth=hclad + 0.2)
Ls2 = xsShallow.add_vstack(name='ridge')
xsShallow.add_hstack(name='shallow_guide', layers=[(Ls1, 1.0), (Ls2, wguide), (Ls1, 1.0)])


# Define the Deep guide cross sectional geometry
xsDeep.layers = epi
xsDeep.background = mat10.air
Ld1 = xsDeep.add_vstack(name='backgrnd', etchdepth=hclad + 1.0)
Ld2 = xsDeep.add_vstack(name='ridge')
xsDeep.add_hstack(name='deep_guide', layers=[(Ld1, 1.0), (Ld2, wguide), (Ld1, 1.0)])
