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
# (c) 2017 Katarzyna Lawniczuk, Ronald Broeke, Xaveer Leijtens
#==============================================================================
"""
Pixapp Package template
"""

from itertools import count
import nazca as nd
import nazca.bb_util as bbu
#import nazca.pdk_template as pdk
import nazca.geometries as geom
import nazca.cfg as cfg
import smart as sp


#==============================================================================
# Package template
#==============================================================================
#annotationlayer = 'AnnotationIO'
annotationlayer = 'package_pin_text'
fiberIOlayer = 'FiberIO'
textheight = 0.15
textmove = -0.08
padtext = 30

# RF pads can be on all sides, but not staggered.
RFname = ['rfN{:02d}', 'rfS{:02d}', 'rfW{:02d}', 'rfE{:02d}']

# DC pads can be on all sides, possibly staggered (in two rows).
# First row
DC0name = ['dc0S{:02d}', 'dc0N{:02d}', 'dc0W{:02d}', 'dc0E{:02d}']
# Second row
DC1name = ['dc1S{:02d}', 'dc1N{:02d}', 'dc1W{:02d}', 'dc1E{:02d}']

class Package1():
    """Class containing a Pixapp Die package template.

    Pad specification values:
        False: no arrows (no placeholder pins)
        True:  show arrows
    """
    instnum = count()

    def __init__(self, name='Pixapp-Package1',
            length=4000, height=4000, cleave=100,
            DCminside=0, RFminside=0,
            DC0edge=0, DC1edge=150, DCpitch=150,
            RFedge=0, RFpitch=150,
            fiberareaL=None, fiberareaR=1000,
            DC0S=True,  DC0N=True,  DC0W=False, DC0E=False,
            DC1S=False, DC1N=False, DC1W=False, DC1E=False,
            RFS=False, RFN=False, RFW=False, RFE=False,
            DCpadcell=None, RFpadcell=None,
            textlayer=annotationlayer,
            ):
        """Construct a package object.

        Args:
            length (float): die length (after cleaving)
            height (float): die height (after cleaving)
            cleave (float): thickness of cleave border
            DCminside (int): minimum distance between pad and side of die
            RFminside (int): minimum distance between pad and side of die
            DC0edge (float): distance between pad and edge of die
            DC1edge (float): distance between pad and edge of die
            DCpitch (float): pitch for DC-pad placement
            RFedge (float): distance between pad and edge of die
            RFpitch (float): pitch for RF-pad placement
            fiberareaL (float): length of position for waveguides to fiber
            fiberareaR (float): length of position for waveguides to fiber
            DC0S (bool): place 0th-row Bottom dc pads
            DC0N (bool): place 0th-row Top dc pads
            DC0W (bool): place 0th-row Left dc pads
            DC0E (bool): place 0th-row Right dc pads
            DC1S (bool): place 1st-row Bottom dc pads
            DC1N (bool): place 1st-row Top dc pads
            DC1W (bool): place 1st-row Left dc pads
            DC1E (bool): place 1st-row Right dc pads
            RFS (x): place Bottom rf pads
            RFN (x): place Top rf pads
            RFW (x): place Left rf pads
            RFE (x): place Right rf pads
            DCpadcell (cell): pad BB to place
            RFpadcell (cell): pad BB to place
            textlayer (str | int): mask layer for text output

        Returns:
            None

        Example:
            Create a package object, create a cell with its method and put it::

            import nazca as nd
            from nazca.demopackager.packages import Package1

            PACKAGE = Package1(name='Package-A1', length=5000, height=4000)
            package.cell().put(0)

            nd.export_plt()
        """

        self.name   = name
        self.length = length
        self.height = height
        self.cleave = cleave

        # Distance of pad bounding box to edge
        self.DC0edge = DC0edge
        self.DC1edge = DC1edge
        self.DCpitch = DCpitch

        # Distance of pad bounding box to edge
        self.RFedge = RFedge
        self.RFpitch = RFpitch

        # Specify which DC pads to put.
        self.DC0 = [DC0S, DC0N, DC0W, DC0E]
        self.DC1 = [DC1S, DC1N, DC1W, DC1E]

        # Specify which RF pads to put.
        self.RF = [RFS, RFN, RFW, RFE]

        # minimum distance to sides of chip (pads are placed from the
        # center of the sides of the die)
        self.DCminside = DCminside
        self.RFminside = RFminside

        self.fiberareaL = fiberareaL
        self.fiberareaR = fiberareaR
        self.textlayer = textlayer

        if DCpadcell is None:
            self.DCpadcell = sp.pad_dc_lw(width=70,length=90)
        else:
            self.DCpadcell = DCpadcell

        self.RFpadcell = RFpadcell

        if self.textlayer is None:
            self.textlayer = 'AnnotationIO'
        else:
            self.textlayer = textlayer

        self.arrow = bbu.make_pincell(layer='package_pin')


    def die_size(self, length=4000, height=4000):
        """Set size of the package area.

        Args:
            length (float): length of the package in um
            heigth (float): height of the package in um

        Returns:
            None
        """
        self.length = length
        self.height = height


    def cell(self):
        """Create a Cell with DC and RF postitions.

        Returns:
            Cell: cell for packaging"""
        # Optical positions come from foundry.
        l = self.length - self.cleave
        h = self.height - self.cleave

        with nd.Cell(name='{}_{}'.format(self.name, next(self.instnum))) as C:

            # Number of DC pads
            nNS = int((l-2*self.DCminside)/self.DCpitch) + 1
            nEW = int((h-2*self.DCminside)/self.DCpitch) + 1
            nDC = [nNS, nNS, nEW, nEW]

            # tuple in pins: start position of dc pins on all four sides.
            dc0 = (
                # Bottom
                nd.Pin('p0B').\
                    put(l/2-(nNS/2-0.5)*self.DCpitch, self.DC0edge, 90),
                # Top
                nd.Pin('p0T').\
                    put(l/2+(nNS/2-0.5)*self.DCpitch,  h-self.DC0edge, 270),
                # Left
                nd.Pin('p0L').\
                    put(self.DC0edge, h/2+(nEW/2-0.5)*self.DCpitch, 0),
                # Right
                nd.Pin('p0R').\
                    put(l-self.DC0edge, h/2-(nEW/2-0.5)*self.DCpitch, 180),
            )

            # Number of RF pads
            nNS = int((l-2*self.RFminside)/self.RFpitch)
            nEW = int((h-2*self.RFminside)/self.RFpitch)
            nRF = [nNS, nNS, nEW, nEW]

            # tuple in pins: start position of dc pins on all four sides.
            rf0 = (
                # Bottom
                nd.Pin('rfS').\
                    put(l/2-(nNS/2-0.5)*self.RFpitch, self.RFedge, 90),
                # Top
                nd.Pin('rfN').\
                    put(l/2+(nNS/2-0.5)*self.RFpitch, h-self.RFedge, 270),
                # Left
                nd.Pin('rfW').\
                    put(self.RFedge, h/2+(nEW/2-0.5)*self.RFpitch, 0),
                # Right
                nd.Pin('rfE').\
                    put(l-self.RFedge, h/2-(nEW/2-0.5)*self.RFpitch, 180),
            )

            # Place DC0 pins
            for pin, n, name, arrow in zip(dc0, nDC, DC0name, self.DC0):
                if not arrow:
                    continue # Don't generate the pins
                for i in range(0, n):
                    loc = nd.Pin(name.format(i)).\
                        put(pin.move(0,-i*self.DCpitch,0))
                    self._put_arrow(loc, name, i)

            # Place DC1 pins
            for pin, n, name, arrow in zip(dc0, nDC, DC1name, self.DC1):
                if not arrow:
                    continue # Don't generate the pins
                for i in range(0, n-1):
                    loc = nd.Pin(name.format(i)).\
                        put(pin.move(self.DC1edge-self.DC0edge,
                            -(i+0.5)*self.DCpitch, 0))
                    self._put_arrow(loc, name, i)

            # Place RF pins
            for p, n, name, arrow in zip(rf0, nRF, RFname, self.RF):
                if not arrow:
                    continue # Don't generate the pins
                for i in range(0, n):
                    loc = nd.Pin(name.format(i)).\
                        put(p.move(0,-i*self.RFpitch,0))
                    self._put_arrow(loc, name, i)

            #Fiber positions
            faalen = self.cleave+50
            if self.fiberareaR:
                dist_edgey = self.height-self.cleave-self.fiberareaR
                farea = geom.box(width=self.fiberareaR, length=faalen)
                poly = nd.Polygon(layer=fiberIOlayer, points=farea)
                poly.put(self.length-0.5*self.cleave,
                    0.5*dist_edgey+0.5*self.fiberareaR, 0)
            if self.fiberareaL:
                dist_edgey = self.height-self.cleave-self.fiberareaL
                farea = geom.box(width=self.fiberareaL, length=faalen)
                poly = nd.Polygon(layer=fiberIOlayer, points=farea)
                poly.put(-faalen-0.5*self.cleave,
                    0.5*dist_edgey+0.5*self.fiberareaL, 0)

        return C


    @property
    def maxDCcount(self):
        """Get max number of DC ports that fit in a package

        Returns:
            int: maximum number of DC pads
        """
        return round((self.length-2*self.DCx)/self.DCpitch)


    @property
    def maxRFcount(self):
        """Get max number of RF ports that fit in a package

        Returns:
            int: maximum number of RF pads
        """
        return round((self.length-2*self.RFy)/self.RFpitch)


    def _put_arrow(self, loc, name, i):
        self.arrow.put(loc)
        t = nd.text(name.format(i), layer=annotationlayer,
            height=textheight, align='rc')
        t.put(loc.move(textmove))

    def DC_pad(self, pin):
        """Standard DC pad for pixapp layout.

        The current pin will be set at rc of the pad.

        Args:
            pin (Node): pin to connect standard DC pad to.

        Returns:
            Cell: cell with DC pad and name of the pad
        """
        cfg.cp = pin
        pnpos = pin.name.find('_dc')
        padname = pin.name[pnpos+1:pnpos+7]
        with nd.Cell('DC-pad-'+padname, instantiate=False) as C:
            pad = self.DCpadcell.put('rc')
            pad.raise_pins()
            C.default_pins('rc', 'c0')
            txt = '{}\n{}'.format(padname[0:3], padname[3:])
            nd.text(txt, height=30, align='cc', layer=self.textlayer).\
                    put(pad.pin['c0'].move(0,0,-90))
        return C


    def RF_pad(self, pin):
        """Standard RF pad for pixapp layout.

        The current pin will be set at rc of the pad.

        Args:
            pin (Node): pin to connect standard RF pad to.

        Returns:
            Cell: cell with RF pad and name of the pad
        """
        cfg.cp = pin
        pnpos = pin.name.find('_rf')
        padname = pin.name[pnpos+1:pnpos+6]
        with nd.Cell('RF-pad-'+padname) as C:
            pad = sp.gsg_pad(length_pad=90, width_pad_sig=70,
                    width_pad_gnd=70, gap2=80+30).put('rc')
            pad.raise_pins()
            C.default_pins('rc', 'c1')
            txt = '{}\n{}'.format(padname[0:2], padname[2:])
            nd.text(txt, height=30, align='ct', layer=self.textlayer).\
                    put(pad.pin['c1'].move(-7,0,-90))
        return C


