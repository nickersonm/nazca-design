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
"""
Created on Mon May 29 10:11:41 2017

@author: katarzyna
"""

#==============================================================================
# (c) 2016-2017 Katarzyna Lawniczuk, Ronald Broeke
#==============================================================================

from math import tan, radians
import nazca as nd
import nazca.geometries as geom
import nazca.demofab.pdk as bbbb


instant = True # TODO: False displaces BBBB!!!


nd.add_layer(name='asd1', layer=1, accuracy=0.001)
nd.add_layer(layer=200)
nd.add_layer(layer=54)


nd.add_xs('wb_soa', 1, growx=0)
nd.add_xs('wb_soa', 2, growx=8)
nd.add_xs('wb_eopm', 1, growx=0)
nd.add_xs('wb_eopm', 2, growx=3)
nd.add_xs('wb_eopm', 3, growx=20)

def soa(length=100):
    """Parametrized SOA whitebox generation."""
    angle = 70.0
    width_active = 30.0
    width = 2.0
    name = 'wb_shallow.soa_{}'.format(length)

    with nd.Cell(name=name, instantiate=False) as C:
        bb = bbbb.soa(length=length).put(0,0,0)
        bb.raise_pins(namesin=['a0'])

        #active
        paral = geom.parallelogram2(
            length=length, height=width_active,
            angle=angle, position=2,
            shift=(0+0.5*width_active/tan(radians(angle)), 0))
        active = nd.Polygon(layer=50, points=paral).put(0,0,0)

        #pad
        rect = geom.rectangle(length=length, height=200,
             position=2, shift=(0, 30))
        metal = nd.Polygon(layer=10, points=rect).put(0,0,0)

        #waveguide
        nd.strt(width=3, length=length, xs='wb_soa').put(0, 0, 0)

    return C


def eopmdp(length=500, pad=1):
    """Parametrized EOPM whitebox generation."""
    width = 55.0
    padl = 100.0
    padh = 100.0
    name = 'wb_deep.eopm_{}_{}'.format(length, pad)

    with nd.Cell(name=name, instantiate=False) as C:
        bb = bbbb.eopm(length=length, pad=pad).put(0,0,0)
        bb.raise_pins(namesin=['a0'])

        if pad == 0:
            outline = [
                (0, 0.5*width), (length, 0.5*width),
                (length, -0.5*width), (0, -0.5*width)
                ]
            nd.Polygon(layer=11, points=outline).put(0,0,0)
            nd.Polygon(layer=52, points=outline).put(0,0,0)
            nd.strt(width=1.5, length=length, xs='wb_eopm').put(0, 0, 0)

        elif pad==1:
            outline = [(0, -0.5*width), (length, -0.5*width),
                       (length, 0.5*width), (length/2-padl/2, 0.5*width),
                       (length/2-padl/2, 0.5*width+padh), (length/2+padl/2, 0.5*width+padh),
                       (length/2+padl/2, 0.5*width), (0, 0.5*width)]
            nd.Polygon(layer=11, points=outline).put(0,0,0)
            nd.Polygon(layer=52, points=outline).put(0,0,0)
            nd.strt(width=1.5, length=length, xs='wb_eopm').put(0, 0, 0)

        elif pad==2:
            offset=10.0
            outline = [
                (0, -0.5*width), (length, -0.5*width),
                (length, 0.5*width), (length+offset, 0.5*width),
                (length+offset, 0.5*width+padh), (length+offset-padl, 0.5*width+padh),
                (length+offset-padl, 0.5*width), (padl-offset, 0.5*width),
                (padl-offset, 0.5*width+padh), (-offset, 0.5*width+padh),
                (-offset, 0.5*width), (0, 0.5*width)
                ]
            nd.Polygon(layer=11, points=outline).put(0,0,0)
            nd.Polygon(layer=52, points=outline).put(0,0,0)
            nd.strt(width=1.5, length=length, xs='wb_eopm').put(0, 0, 0)

    return C


#map black box BBs to white box functions

black2whiteMap= {
    'shallow.soa': soa,
    'deep.eopm': eopmdp
    }


#cells_to_remove = ('stmetal', 'stshallow')
