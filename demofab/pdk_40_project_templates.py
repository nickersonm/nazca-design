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
**DEMO PDK project definition.**

This module is intended to define the cell templates of a technology,
e.g. the cells for an MPW.
It may contains cells that define **optical or electrical IO** functionality of
the die with the external world, e.g. fiber-coupling or wire-bonding.
"""

import nazca as nd
from . import pdk_20_technology as tech
from . import pdk_30_BB_library as bblib
import nazca.geometries as geom
import nazca.cfg as cfg
import nazca.bb_util as bbu


package_pin_text = 'package_pin_text' #annotation layer in cell
#==============================================================================
# Die template
#==============================================================================
class DesignArea():
    """
    Class containing a Die - Design Area template for Demofab.
    """


    def __init__(self, name='DesignArea', length=4000, height=4000, cleave=100,
                 coatingL='NO', coatingR='NO', pitch=50):
        """Construct a DesignArea object, i.e. a foundry cell on a wafer

        The design area indicates the area that is diced/cleaved later into a die.
        This area can also be referred to as "project" in case of an MPW
        multi-project wafer.

        Args:
            name (str): name of the cell containing the design area
            length (float): physical length of the die after dicing in um
            width (float): physical width of the die after dicingin um
            cleave (float): width of the cleave/dice area in um
            coatingL (str): coating type on left side of the die: AR | HR | NO | DC: (default = 'NO')
            coatingR (str): coating type on right side of the die: AR | HR | NO | DC: (default = 'NO')
            pitch (float): pitch of optical io (default = 50)

        Returns:
            None
        """

        self.name   = name
        self.length = length # die size
        self.height = height # die size
        self.cleave = cleave
        self.coatingR = coatingR
        self.coatingL = coatingL
        self.pitch = pitch
        self.arrow = bbu.make_pincell(layer='package_pin')


    def size(self, length=4000, height=4000):
        """Set size of the design area.

        Args:
            length (float): physical length of the die after dicing in um
            width (float): physical width of the die after dicingin um

        Returns:
            None
        """
        self.length = length
        self.height = height


    def coating(self, coatingL='NO', coatingR='NO'):
        """Set coating at the left and at the right side of the design area.

        Args:
            coatingL (str): coating type on left side of the die: AR | HR | NO | DC: (default = 'NO')
            coatingR (str): coating type on right side of the die: AR | HR | NO | DC: (default = 'NO')

        Returns:
            None
        """
        self.coatingL = coatingL
        self.coatingR = coatingR


    def cell(self):
        """Create a Project Cell with IO postition pins.

        Returns:
            Cell: cell containing design area and cleave boundaries
        """
        with nd.Cell(name=self.name) as C:
            # put the dice/cleave border
            frame = geom.frame(sizew=self.cleave, sizel=self.length, sizeh=self.height)
            nd.Polygon(layer='DiceArea', points=frame).put(0)

            nd.Pin('left', show=False).put(-0.5*(self.cleave), -0.5*self.cleave+0.5*self.height, 0)
            bbu.put_boundingbox('left', self.length, self.height)

            # put the pins
            IOcountmax = round((self.height-self.cleave-self.pitch)/self.pitch)
            pin_indent = 0
            ioy0 = 50

            for no in range(0, IOcountmax):
                pinIDL = 'ioL{:03d}'.format(no)
                pinIDR = 'ioR{:03d}'.format(no)
                if no % 2 == 0:
                    angle = 7
                else:
                    angle = 0

                #IOs left
                p = nd.Pin(name=pinIDL, xs='Shallow', type=angle).put(-self.cleave+pin_indent, ioy0+no*self.pitch, angle)
                self.arrow.put(p)
                nd.text(pinIDL, layer=package_pin_text, height=0.15, align='rc').put(p.move(-0.1))

                #IOs right
                p = nd.Pin(name=pinIDR, xs='Shallow', type=angle).put(-pin_indent+self.length, ioy0+no*self.pitch, 180+angle)
                self.arrow.put(p)
                nd.text(pinIDR, layer=package_pin_text, height=0.15, align='rc').put(p.move(-0.575, 0, 180))

            # Adding coating
            options = {'AR': 'Coating_AR', 'HR': 'Coating_HR',
                       'NO': 'Coating_NO', 'DC': 'Coating_DC'}
            coating_on_chip = 100

            box_coatno = [
               (0, 0), (self.cleave/2, 0),
               (self.cleave/2, self.height), (0, self.height)
            ]

            box_coatother = [
                (0, 0), (coating_on_chip+0.5*self.cleave, 0),
                (coating_on_chip+0.5*self.cleave, self.height), (0, self.height)
            ]

            pinL = nd.Pin().put(-0.5*self.cleave, -0.5*self.cleave)
            pinR = nd.Pin().put(self.length-0.5*self.cleave, self.height-0.5*self.cleave, 180)

            coatings = {'L': (self.coatingL, pinL), 'R': (self.coatingR, pinR) }

            for coat, pin in coatings.values():
                if coat in options.keys():
                    layer = options[coat]
                else:
                    print('Available coating options: AR | HR | NO | DC. Default is NO.')

                if coat == 'NO':
                    box = box_coatno
                else:
                    box = box_coatother
                nd.Polygon(layer=layer, points=box).put(pin)

            return C
        return cell


    def get_height(self):
        """
        Returns:
            float: height of the design area
        """
        return self.height


    def get_length(self):
        """
        Returns:
            float: length of the design area
            """
        return self.length


    def get_cleave(self):
        """
        Returns:
            float: cleave of the design area
        """
        return self.cleave


    def get_NIO(self):
        """
        Returns:
            int: number of optical IOs
        """
        return round((self.height-self.cleave-self.pitch)/self.pitch)


#==============================================================================
# Standard IOs
#==============================================================================



def presetio(pin, bend=False, deep=False):
    """Returns a preset IO cell based on the attributes of the provided pin.

    The cp will be updated to cp = pin.

    Args:
        bend (bool): Add a bend to horizontal alignment in case of an angled input guide
        transition (bool): Add a transition element to the deep waveguide xsection

    Returns:
        Cell: cell containing the io element
    """

    if pin.xs == 'Shallow':
        with nd.Cell('io', instantiate=False) as C:
            if pin.type == 7:
                bend = True
                shape = 'angled'
            else:
                bend = False
                shape = 'tapered'
            bb = bblib.io(bend=bend, shape=shape, deep=deep).put()
            bb.raise_pins()
    elif pin.xs == 'MetalDC':
        with nd.Cell('io', instantiate=False) as C:
            bb = bblib.bb_io_electrical.pad_dc().put()
            bb.raise_pins(['c0', 'c0'], ['a0', 'b0'])
    elif pin.xs == 'MetalRF':
        with nd.Cell('io', instantiate=False) as C:
            bb = bblib.bb_io_electrical.pad_rf().put()
            bb.raise_pins(['c0', 'c0'], ['a0', 'b0'])
    else:
        with nd.Cell('io', instantiate=False) as C:
           pass # empty cell.

    nd.cfg.cp = pin
    return C


