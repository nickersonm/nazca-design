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
DEMO Package template
"""

import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom


#==============================================================================
# Package template
#==============================================================================
annotationlayer = 'package_pin_text'
fiberIOlayer = 'FiberIO'
textheight = 0.15
textmove = -0.08
RFname = 'rf{:03d}'
DCtopname = 'dcT{:03d}'
DCbotname = 'dcB{:03d}'

pin_naming = {
    'leftDC':'dcL{:03d}',
    'topDC':'dcT{:03d}',
    'rightDC':'dcR{:03d}',
    'bottomDC':'dcB{:03d}',
    'leftRF':'rfL{:03d}',
    'topRF':'rfT{:03d}',
    'rightRF':'rfR{:03d}',
    'bottomRF':'rfB{:03d}',
    }


class Package1():
    """
    Class containing a Demo package template.
    """

    def __init__(self, name='Package-Demo',
            die_length=4000, die_height=4000, cleave=100,
            DCside=250, DCedge=20, DCpitch=180, DCcount=None, DCcenter=False,
            RFside=500, RFedge=25, RFpitch=200, RFcount=None, RFcenter=False,
            show_fiberarea=True, fiberarea=1000,
            double_row_DC=False, DCx_doublerow=None, DCy_doublerow=None,
            double_row_RF=False, RFx_doublerow=None, RFy_doublerow=None,
            textlayer=None):
        """Construct a Package1 object.

        Args:
            name (str): name of the package cell
            lenght (float): length of the package
            width (float): width of the package
            cleave (float): thickness of cleave/dice area

            DCside (float): distance from side for DC pads
            DCedge (float): distance from the bond-edge for DC pads

            RFside (float): distance from side for RF pads
            RFedge (float): distance from the bond-edge for RF pads

        Returns :
            None
        """
        self.name = name
        self.die_length = die_length
        self.die_height = die_height
        self.cleave = cleave

        self.DCside = DCside
        self.DCedge = DCedge
        self.DCpitch = DCpitch
        self.DCcenter = DCcenter
        self.DCcount = DCcount

        self.RFedge = RFedge
        self.RFside = RFside
        self.RFpitch  = RFpitch
        self.RFcenter = RFcenter
        self.RFcount = RFcount

        self.show_fiberarea = show_fiberarea
        self.fiberarea = fiberarea
        self.double_row_DC = double_row_DC
        self.DCx_doublerow = DCx_doublerow
        self.DCy_doublerow = DCy_doublerow
        self.double_row_RF = double_row_RF
        self.RFx_doublerow = RFx_doublerow
        self.RFy_doublerow = RFy_doublerow
        self.textlayer = textlayer

        if self.textlayer is None:
            self.textlayer = annotationlayer
        else:
            self.textlayer = textlayer

        self.arrow = bbu.make_pincell(layer='package_pin')


    def die_size(self, die_length=4000, die_height=4000):
        """Set size of the package area.

        Args:
            length (float): length of the package in um
            heigth (float): height of the package in um

        Returns:
            None
        """
        self.die_length = die_length
        self.die_height = die_height


    def cell(self):
        """Create a Cell with DC, RF, Optical fiber postitions.

        Returns:
            Cell: package cell
        """
        with nd.Cell(name=self.name) as C:

            #Placing PAD IOs
            IOcountmax = round((self.die_length-2*self.DCside)/self.DCpitch)+1
            if self.DCcount is not None and self.DCcount < IOcountmax:
                IOcount = self.DCcount
            else:
                IOcount = IOcountmax

            #DC PAD positions
            for n in range(0, IOcount):
                pinIDT = DCtopname.format(n)
                pinIDB = DCbotname.format(n)

                if self.DCcenter:
                    dist_edgex = 0.5*(self.die_length -(IOcount-1)*self.DCpitch)
                else: dist_edgex = self.DCside

                #Top DC ports
                p = nd.Pin(pinIDT).put(
                    dist_edgex + n*self.DCpitch - 0.5*self.cleave,
                    self.die_height-self.DCedge-0.5*self.cleave,
                    -90)
                self.arrow.put(p)
                text = nd.text(pinIDT, layer=annotationlayer,
                    height=textheight, align='rc')
                text.put(p.move(textmove))

                #Bottom DC ports
                p = nd.Pin(pinIDB).put(
                    dist_edgex + n*self.DCpitch - 0.5*self.cleave,
                    self.DCedge-0.5*self.cleave,
                    90)
                self.arrow.put(p)
                text = nd.text(pinIDB, layer=annotationlayer,
                    height=textheight, align='rc')
                text.put(p.move(textmove))

            #DC PAD positions double row
            if self.double_row_DC:
                for ndr in range(0, IOcount-1):
                    n = ndr+IOcountmax
                    pinIDT = DCtopname.format(n)
                    pinIDB = DCbotname.format(n)

                    if self.DCcenter == True:
                        dist_edgex = self.die_length - self.cleave -\
                            (IOcount-1)*self.DCpitch-2*self.DCside
                    else: dist_edgex = 0

                    if self.DCx_doublerow == None:
                        self.DCx_doublerow = 0.5*self.DCpitch

                    if self.DCy_doublerow == None:
                        self.DCy_doublerow = 1.0*self.DCpitch

                    #Top DC ports
                    p = nd.Pin(pinIDT).put(
                        ndr*self.DCpitch+self.DCside+0.5*dist_edgex+self.DCx_doublerow,
                        self.die_height-self.cleave-self.DCedge-self.DCy_doublerow,
                        -90)
                    self.arrow.put(p)
                    text = nd.text(pinIDT, layer=annotationlayer,
                        height=textheight, align='rc')
                    text.put(p.move(textmove))

                    #Bottom DC ports
                    p = nd.Pin(pinIDB).put(
                        ndr*self.DCpitch+self.DCside+0.5*dist_edgex+self.DCx_doublerow,
                        self.DCedge+self.DCy_doublerow,
                        90)
                    self.arrow.put(p)
                    text = nd.text(pinIDB, layer=annotationlayer,
                        height=textheight, align='rc')
                    text.put(p.move(textmove))

            #RF PAD positions
            IOcountmax = round((self.die_height-2*self.RFside)/self.RFpitch)
            if self.RFcount is not None and self.RFcount < IOcountmax:
                IOcount = self.RFcount
            else:
                IOcount = IOcountmax

            for n in range(0, IOcount):
                piRFn = RFname.format(n)

                if self.RFcenter:
                    dist_edgey = self.die_height-self.cleave -\
                        (IOcount-1)*self.RFpitch-2*self.RFside
                else: dist_edgey = 0

                p = nd.Pin(piRFn).put(
                    self.die_length-self.cleave-self.RFedge,
                    n*self.RFpitch+self.RFside+0.5*dist_edgey,
                    180)
                self.arrow.put(p)
                text = nd.text(piRFn, layer=annotationlayer,
                    height=textheight, align='rc')
                text.put(p.move(textmove))

            #RF PAD positions double row
            if self.double_row_RF:
                for ndr in range(0, IOcount-1):
                    n = ndr+IOcountmax
                    piRFn = RFname.format(n)

                    if self.RFcenter == True:
                        dist_edgey = self.die_height-self.cleave -\
                            (IOcount-1)*self.RFpitch-2*self.RFside
                    else: dist_edgey = 0

                    if self.RFx_doublerow == None:
                        self.RFx_doublerow = 1.0*self.RFpitch

                    if self.RFy_doublerow == None:
                        self.RFy_doublerow = 0.5*self.RFpitch

                    p = nd.Pin(piRFn).put(\
                        self.die_length-self.cleave-self.RFedge-self.RFx_doublerow,
                        ndr*self.RFpitch+self.RFside+0.5*dist_edgey+self.RFy_doublerow,
                        180)
                    self.arrow.put(p)
                    text = nd.text(piRFn, layer=annotationlayer,
                        height=textheight, align='rc')
                    text.put(p.move(textmove))

            #Fiber positions
            dist_edgey = self.die_height-self.cleave-self.fiberarea
            faa_len = self.cleave+50

            if self.show_fiberarea:
                farea = geom.box(width=self.fiberarea, length=faa_len)
                poly = nd.Polygon(layer=fiberIOlayer, points=farea)
                poly.put(-faa_len-0.5*self.cleave, 0.5*dist_edgey+0.5*self.fiberarea, 0)
            else: print('Set show_fiberarea == True to see the fiber allowed area for this package.')

        return C


    @property
    def maxDCcount(self):
        """Get max number of DC ports that fit in a package

        Returns:
            int: maximum number of DC pads
        """
        return round((self.die_length-2*self.DCside)/self.DCpitch)


    @property
    def maxRFcount(self):
        """Get max number of RF ports that fit in a package.

        Returns:
            int: maximum number of RF pads
        """
        return round((self.die_length-2*self.RFside)/self.RFpitch)



class Package2():
    """Class containing a Demo package template.

    This package provides flexible placement of RF and/or DC pins location around
    the all edges of the DIE.
    """

    def __init__(self, name='Package-Demo',
            die_length=4000, die_height=4000, cleave=100,
            pads=None,
            show_fiberarea=True, fiberarea=1000):
        """Construct a Package2 object.

        Args:
            name (str): name of the package cell
            die_length (float): length of the package
            die_height (float): width of the package
            cleave (float): thickness of cleave/dice area
            pads (list of dict): list of arrays decribed by a dictionary
            show_fiberarea (bool): display the fiber area (default = True)
            fiberarea (float): height of the fiber area
            textlayer=None):

        Returns:
            None
        """

        if pads is None:
            DC1 = {
                'edge': 'top',
                'type': 'DC',
                'count': 25,
                'count0': 0,
                'pitch': 200,
                'center': True,
                'edge_sep_side': 150,
                'edge_sep_front': 100}
            DC2 = {
                'edge': 'bottom',
                'type': 'DC',
                'count': 25,
                'count0': 0,
                'pitch': 200,
                'center': True,
                'edge_sep_side': 100,
                'edge_sep_front': 100}
            RF1 = {
                'edge': 'left',
                'type': 'RF',
                'count': 16,
                'count0': 0,
                'pitch': 400,
                'center': True,
                'edge_sep_side': 150,
                'edge_sep_front': 100}

            self.pads = [DC1, DC2, RF1]
        else:
            self.pads = pads

        self.name = name
        self.die_length = die_length
        self.die_height = die_height
        self.cleave = cleave

        self.show_fiberarea = show_fiberarea
        self.fiberarea = fiberarea
        self.arrow = bbu.make_pincell(layer='package_pin')


    def die_size(self, die_length=4000, die_height=4000):
        """Set size of the package area.

        Args:
            length (float): length of the package in um
            heigth (float): height of the package in um

        Returns:
            None
        """
        self.die_length = die_length
        self.die_height = die_height


    def cell(self):
        """Create a Cell with DC, RF, Optical fiber postitions.

        Returns:
            Cell: package cell
        """
        with nd.Cell(name=self.name) as C:

            for arr in self.pads:
                if 'count0' not in arr.keys():
                    arr['count0'] = 0
                count = arr['count']
                #Placing PAD IOs
                if arr['center'] == True:
                    IOcountmax = round(\
                        (self.die_length - 2*arr['edge_sep_side']) / arr['pitch'])
                else:
                    IOcountmax = round(\
                        (self.die_length - 1*arr['edge_sep_side']) / arr['pitch'])


                if count is not None and count < IOcountmax:
                    IOcount = count
                else:
                    IOcount = IOcountmax

                pin = arr['edge'] + arr['type']
                if arr['edge'] == 'left':
                    x0, y0, a0 = 0, 0, 90
                elif arr['edge'] == 'top':
                    x0, y0, a0 = 0, self.die_height, 0
                elif arr['edge'] == 'right':
                    x0, y0, a0 = self.die_length, self.die_height, -90
                elif arr['edge'] == 'bottom':
                    x0, y0, a0 = self.die_length, 0, 180
                else:
                    print("Edge not recognized in Package2: '{}'".format(arr['edge']))

                #PAD positions
                for n in range(0, IOcount):
                    pinID = pin_naming[pin].format(n+arr['count0'])
                    if arr['center']:
                        if arr['edge'] in ['top', 'bottom']:
                            dist_edge = self.die_length -\
                                (IOcount-1)*arr['pitch'] - 2*arr['edge_sep_side']
                        elif arr['edge'] in ['left', 'right']:
                            dist_edge = self.die_height -\
                                (IOcount-1)*arr['pitch'] - 2*arr['edge_sep_side']
                    else: dist_edge = 0

                    #place metal io pins
                    start = nd.Pin().put(x0, y0, a0)
                    p = nd.Pin(pinID).put(start.move(
                        n*arr['pitch'] + arr['edge_sep_side'] + 0.5*dist_edge,
                        -arr['edge_sep_front'],
                        -90))
                    self.arrow.put(p)
                    text = nd.text(pinID, layer=annotationlayer,
                        height=textheight, align='rc')
                    text.put(p.move(textmove))

            #Fiber positions
            dist_edgey = self.die_height-self.cleave-self.fiberarea
            faa_len = self.cleave+50

            if self.show_fiberarea:
                farea = geom.box(width=self.fiberarea, length=faa_len)
                poly = nd.Polygon(layer=fiberIOlayer, points=farea)
                poly.put(-faa_len-0.5*self.cleave, 0.5*dist_edgey+0.5*self.fiberarea, 0)
            #else: print('Set show_fiberarea == True to see the fiber allowed area for this package.')

        return C


    @property
    def maxDCcount(self):
        """Get max number of DC ports that fit in a package

        Returns:
            int: maximum number of DC pads
        """
        return round((self.die_length-2*self.DCside)/self.DCpitch)


    @property
    def maxRFcount(self):
        """Get max number of RF ports that fit in a package.

        Returns:
            int: maximum number of RF pads
        """
        return round((self.die_length-2*self.RFside)/self.RFpitch)

