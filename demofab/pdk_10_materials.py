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
#==============================================================================
# (c) 2017 Katarzyna Lawniczuk, Ronald Broeke
#==============================================================================

"""
**DEMO PDK material settings.**

This module is intended to define the material models and materials used in the DEMOfab.
"""

import nazca as nd


#==============================================================================
# Materials
#==============================================================================
#define material refractive functions:
def N(wl, pol, **kwargs):
    return 3.5+0.1*wl+pol

#Define materials. Use value or functions to set the refractive index:
InP  = nd.xsection.Material(Nmat=3.2, name='InP',  rgb=(0.0, 0.4, 0.9))
Q125 = nd.xsection.Material(Nmat=N,   name='Q125', rgb=(0.0, 0.8, 0.3))
air  = nd.xsection.Material(Nmat=1.0, name='air',  rgb=(0.95, 0.95, 1.0))

