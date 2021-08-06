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
#
# 2017 (c) Ronald Broeke
#-----------------------------------------------------------------------
# -*- coding: utf-8 -*-
"""
Nazca module for interconnecting guides.
"""

import numpy as np
from numpy import sign, sin, cos, radians, sqrt
import math as m

import nazca as nd
import nazca.cfg as cfg
from nazca.logging import logger
import nazca.cp as cp
import nazca.geometries as geom
import nazca.bb_util as bbu
import nazca.trace as trace
from nazca.netlist import interconnect_logger


cfg.interconnects = {} # store all interconnects for reference/documentation
cfg.interconnect_errcnt = 0  # cnt warnings and errors to locate and catch them better.
num = -1


def posRad(a):
    """Clip angle to [0, 2pi>.

    Args:
        a (float): angle in radians

    Returns:
        float: angle in radians
    """
    return a % (2 * m.pi)


def negRad(a):
    """Clip angle to <-2pi, 0]

    Args:
        a (float): angle in radians

    Returns:
        float: angle in radians
    """
    return 0 if a == 0 else a % (2 * m.pi) - 2 * m.pi


def zeroRad(a):
    """Clip angle to <-pi, pi].

    Args:
        a (float): angle in radians

    Returns:
        float: angle in radians
    """
    a = posRad(a)
    if a > m.pi:
        a -= 2 * m.pi
    return a


def posDeg(a):
    """Clip angle to [0, 360>.

    Args:
        a (float): angle in degrees

    Returns:
        float: angle in degrees
    """
    return a % 360


def negDeg(a):
    """Clip angle to <-360, 0].

    Args:
        a (float): angle in degrees

    Returns:
        float: angle in degrees
    """
    return 0 if a == 0 else a % 360 - 360


def zeroDeg(a):
    """Clip angle in degrees between <-180, 180].

    Args:
        a (float): angle in degrees

    Returns:
        float: angle in degrees
    """
    a = a % 360
    return a if a <= 180 else a - 360


def ic_exception(msg=''):
    """Raise interconnect exception if conditions apply."""
    if cfg.interconnect_raise:
        if cfg.interconnect_msgnum is None\
        or cfg.interconnect_msgcnt == cfg.interconnect_msgnum:
            raise Exception(msg)
            # TIP: use the traceback (2 steps up) to find the cause of the exception.


class Interconnect():
    """Interconnect class for drawing waveguides and metal lines.

    An Interconnect object can be configured to match a specifc foundry.
    This includes properties like width, xsection, straight-bend offset and others.

    Example:
        Create two Interconnect objects, each for a different kind of waveguide::

            import nazca as nd
            import nazca.interconnects as IC

            ic1 = IC.Interconnect(width=2.0, radius=20)
            ic2 = IC.Interconnect(width=1.0, radius=50)

            # use the interconnect objects to create layout:
            ic1.strt(length=10).put(0)
            ic1.bend(angle=45).put()

            ic2.strt(length=20).put(20)
            ic2.bend(angle=45).put()

            nd.export_plt()
    """
    num = 0

    def __init__(
        self,
        radius=None,
        width=None,
        angle=90,
        xs=None,
        layer=None,
        adapt_width=False,
        adapt_xs=False,
        instantiate=False,
        pinstyle=None,
        offset=None,
        varname=None,
        doc='',
        PCB=False,
        modes=None,
    ):
        """Construct an Interconnect object.

        If a xsection value <xs> is provided for 'xs' then values for
        'radius' and 'width' are copied from the <xs> attributes if present.

        Any values radius and/or width that are explicitly set in __init__
        take priority over values in <xs>.

        If xs nor radius and/or width are set an interconnect with
        the name 'nazca' is created and the values are set in via the
        the values in the cfg module.

        Inside interconnects an arrow ('pinshape') is placed at the beginning
        and the end pin of the interconnect.

        Args:
            radius (float): default radius in um
            width (float): default waveguide width im um
            angle (float): default angle of a bend (default=90 degrees)
            xs (str): waveguide xsection (default='nazca')
            layer (str): layer to draw interconnect in.
                It is preferred to use a xs rather than a layer for Interconnects
                to store additional information like offset, index, etc.
                (default=None)
            adapt_width (bool): adapt interconnect width to the pin it connects to
                (default=False)
            adapt_xs (bool): adapt interconnect width to the pin it connects to
                (default=False)
            instantiate (bool): instantiate the interconnect cell (default=False)
            pinstyle: visualisation style of the interconnect pin. default=None
                uses the default pinstyle from cfg
            offset (value | func): if not None then overide the xsection's straight-bend
                offset with the provided offset value or function (default=None).
            PCB (bool): if true, set maximum angle in bends equal to 45 deg (default=False).
            modes (list): list of integers, containing the label of netlist modes for this interconnects

        Returns:
            None
        """
        # Optional: store info about the interconnect:
        if varname is None:
            varname = 'noname_{}'.format(Interconnect.num)
        else:
            # register the interconnect if it got a name
            cfg.interconnects[varname] = self
        self.varname = varname
        self.doc = doc
        Interconnect.num += 1

        if radius is None:
            try:
                R = nd.get_xsection(xs).radius
                if R is None:
                    self.radius = cfg.default_xs_radius
                else:
                    self.radius = R
            except:
                self.radius = cfg.default_xs_radius
        else:
            self.radius = radius

        if width is None:
            try:
                W = nd.get_xsection(xs).width
                if W is None:
                    self.width = cfg.default_xs_width
                else:
                    self.width = W
            except:
                self.width = cfg.default_xs_width
        else:
            self.width = width

        if layer is None and xs is None:
            self.xs = cfg.default_xs_name
        else:
            self.xs = xs
        if self.xs is None:
            interconnect_logger("Created an Interconnect object with layer {} "\
                "but no xsection.".format(layer), 'warning')

        if layer is None:
            self.layer = None
        else:
            self.layer = nd.get_layer(layer)

        self.offset = offset

        if pinstyle is None: # then use xsection pin setting as default
            try:
                xsstyle = nd.get_xsection(self.xs).pinstyle
                pinstyle = xsstyle
            except:
                pass # no attribute set in xscetion: do nothing

        if modes is None:
            self.modes = nd.sim.modes
        else:
            self.modes = modes

        self.pinstyle = pinstyle
        self.angle = angle
        self.length = 10
        self.instantiate = instantiate
        self.xya = (100, 100, 10)
        self.line = nd.Tp_straight(xs=self.xs, layer=self.layer, modes=self.modes)
        self.sinecurve = nd.Tp_sinecurve(xs=self.xs, layer=self.layer)
        self.cobra = nd.Tp_cobra(xya=self.xya, width1=None, width2=None, radius1=0,
            radius2=0, offset1=None, offset2=None, xs=self.xs,
            layer=self.layer)
        self.euler_base = nd.Tp_euler(width1=None, width2=None, radius=self.radius,
            xs=self.xs, layer=self.layer, angle=90)
        if PCB:
            self.Farc = nd.Tp_arc2(xs=self.xs, layer=self.layer)
        else:
            self.Farc = nd.Tp_arc(xs=self.xs, layer=self.layer, modes = self.modes)
        self._ptaper = nd.Tp_ptaper(xs=self.xs, layer=self.layer)
        self._taper = nd.Tp_taper(xs=self.xs, layer=self.layer)
        self.adapt_width = adapt_width
        self.adapt_xs = adapt_xs
        self.arrow = bbu.make_pincell(style=self.pinstyle)
        self.max_length = 1e5 # maximum line length in p2l.
        self.pinflip = True

        self.gridpatch = 0.0  # Experimental: add taper the size fo gridpatch to strt.


    def copy(self, ic=None):
        """Create a copy of the Interconnect object.

        The copy is created by initializing a new Interconnect with __init__
        using only the parameters send to init.

        Returns:
            Interconnect: copy of self

        Example:
            Create a new interconnect based on an existing interconnect
            with a different bend radius but otherwise unchanged settings::

                import nazca as nd
                import nazca.interconnects as IC

                ic1 = IC.Interconnect(width=2.0, radius=20)
                ic2 = ic1.copy()
                ic2.radius = 100
        """
        if ic is None:
            ic = self
        return Interconnect(
            radius=ic.radius,
            width=ic.width,
            angle=ic.angle,
            xs=ic.xs,
            layer=ic.layer,
            modes=ic.mode,
            adapt_width=ic.adapt_width,
            adapt_xs=ic.adapt_xs,     
            instantiate=ic.instantiate,
            pinstyle=ic.pinstyle,
            offset=ic.offset,
            varname=ic.varname,
            doc=f"copy: {ic.doc}",
            PCB=ic.PCB,
        )


    def _arc(self, radius=None, width=None, angle=None, xs=None,
            layer=None, offset=None, name=None):
        """Return an arc function which is overloadable in derived classes.

        Interconnect bends will use the arc function returned by this method.
         If a keyword value is not provided explicitly
        the Interconnect default is applied.

        Note that a derived function can't overload a base attribute. Hence
        the arc is here as defined as an Interconnect method.

        Arguments:
            radius (float): bend radius in um
            width (float): waveguide width in um
            angle (float): default bend angle in deg
            xs (str): waveguide xsection
            layer (str): waveguide layer name
            offset (float): straight-bend-offset in um
            name (str): name of the mask element generated

        Returns:
            function: arc function
        """
        if layer is None:
            layer = self.layer
        if xs is None:
            xs = self.xs
        if radius is None:
            radius = self.radius
        if width is None:
            width = self.width
        if angle is None:
            angle = self.angle
        if offset is None:
            offset = self.offset
        return self.Farc(radius=radius, width=width, angle=angle, xs=xs,
            layer=layer, offset=offset, name=name)


    def _getpinin(self, pin):
        """Return the pin as is or the default_in pin if an instance is provided."""
        if isinstance(pin, nd.Instance):
            pin = pin.pin[pin.cnode.cell.default_in]
        return pin


    def _getpinout(self, pin):
        """Return the pin as is or the default_out pin if an instance is provided."""
        if isinstance(pin, nd.Instance):
            pin = pin.pin[pin.cnode.cell.default_out]
        return pin


    def _getwidth(self, pin=None, width=None, xs=None, adapt=False):
        """Return width based on interconnect rules.

        Args:
            pin (Node): pin to connect
            width (float): waveguide width
            xs (str): xsection
            end (bool): True if output output side of taper

        Returns:
            width
        """
        if width is not None:
            return width
        if adapt or self.adapt_width:
            if pin is None:
                pin = cp.here()
            if pin.width is not None:
                width = pin.width
        if width is None:
            return self.width
            #try:
            #    width = nd.get_xsection(xs).width
            #except:
            #    width = self.width
        return width


    def _getradius(self, pin, radius, xs):
        """Obtain the bend radius based on interconnect rules.

        Args:
            pin (Node):
            radius (float):
            xs (str): xsection

        Returns:
            float: radius
        """
        if self.adapt_xs:
            if pin is None:
                pin = cp.here()
            if xs is None:
                try:
                    xs = pin.xs
                except:
                    xs = self.xs
        if radius is None:
            try:
                radius = nd.get_section(xs).radius
            except:
                radius = self.radius
        return radius


    def _getxs(self, pin, xs):
        """Obtain the xs based on interconnect rules.

        Returns:
            str: xsection
        """
        if xs is not None:
            return xs
        if self.adapt_xs:
            if pin is None:
                pin = cp.here()
            try:
                if pin.xs is not None:
                    xs = pin.xs
            except:
                xs = self.xs
                # raise Exception('No xsection defined in pin. Add a xs attribute.')
        if xs is None:
            xs = self.xs
        return xs


    def _p2p_parse(self, pin1, pin2=None, xs=None, width1=None, width2=None,
            radius1=None, radius2=None):
        """Parse point-to-points parameters.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            width1 (float): optional waveguide width in um
            width2 (float): optional waveguide width in um at end
            xs (str): optional xsection
            radius1 (float): optional first bend radius in um
            radius2 (float): optional second bend radius im um

        Returns:
            Node, Node, string, float, float, float:
                pin1, pin2, xs, width1, width2, radius1, radius2
        """
        if pin1 is None:
            pin1 = cp.here()
        if pin2 is None:
            pin1, pin2 = cp.here(), pin1

        # Keep the original p2p pins (identity) where possible to keep netlist
        # connectivity in the right cell scope.
        pin1b, T = nd.parse_pin(pin1)
        if T != [0, 0, 0] or not isinstance(pin1, nd.Node):
            pin1 = pin1b.move(*T)
        pin2b, T = nd.parse_pin(pin2, rot=True, default='out')
        if T != [0, 0, 0] or not isinstance(pin2, nd.Node):
            pin2 = pin2b.move(*T)

        xs = self._getxs(pin1, xs)
        width1 = self._getwidth(pin1, width1, xs)
        width2 = self._getwidth(pin2, width2, xs)
        radius1 = self._getradius(pin1, radius1, xs)
        radius2 = self._getradius(pin1, radius2, xs)
        return pin1, pin2, xs, width1, width2, radius1, radius2


    def _strt_solve():
        """To be implemented."""
        pass


    # TODO: add strt_solve
    def strt(self, length=None, width=None, pin=None, xs=None, edge1=None,
            edge2=None, edgepoints=50, name=None, arrow=True, gridpatch=None):
        """Create a straight waveguide.

        Args:
            length (float): length of guide in um
            width (float): width of guide in um
            pin (Node): optional Node for modeling info
            xs (str): optionals xsection of guide
            layer (int | str): layer number or layername
            edge1 (function): optional function F(t) describing edge1 of the waveguide
            edge2 (function): optional function G(t) describing edge2 of the waveguide
            edgepoints (int): optional number of points for edge1 and edge2 (default=50)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            gridpatch (float): patch gridsnap jumps at grid disconnect 
                of cells with chamfers of size gridpatch. Default=0 is no patch.

        Returns:
            Cell: waveguide element

        Example:
            Create and place a straight waveguide::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.strt(length=20)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        pinflip = not nd.get_xsection(xs).symmetry
        if gridpatch is None:
            gridpatch = self.gridpatch
        if length is None:
            length = self.length
        if name is None:
            name = 'ic_strt'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self.line(length=length, width=width, xs=xs, edge1=edge1,
                edge2=edge2, edgepoints=edgepoints, gridpatch=gridpatch).put(0)
            nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def _bend_solve():
        """To be implemented."""
        pass


    #TODO: add bend_solve()
    def bend(self, radius=None, angle=None, width=None, pin=None, xs=None,
            length1=0, length2=0, name=None, arrow=True, offset=None):
        """Create a bent waveguide (circular arc) with optinal in/out strt sections.

        Args:
            radius (float): radius at the center line of the arc in um
            width (float): width of the arc in um
            angle (float): angle of arc in degree (default=90)
            pin (Node): optional Node for modeling info
            xs (str): optiinal xsection of bend
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            offset (float): optional new offset for this bend only.
            length1 (float): length of an optional straight section before the bend
            length2 (float): length of an optional straight section after the bend

        Returns:
            Cell: circularly-bent waveguide element

        Example:
            Create and place a bend::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.bend(angle=45)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        radius = self._getradius(pin, radius, xs)
        pinflip = not nd.get_xsection(xs).symmetry

        if offset is None:
            offset = self.offset

        if angle is None:
            angle = self.angle
        if name is None:
            name = 'ic_bend'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            if length1 > 0:
                e1 = self.line(length=length1, width=width, xs=xs).put(0)
                self._arc(radius=radius, width=width, angle=angle, xs=xs, offset=offset).put()
            else:
                e1 =self._arc(radius=radius, width=width, angle=angle, xs=xs, offset=offset).put()
            if length2 > 0:
                self.line(length=length2, width=width, xs=xs).put()
            nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()

        if pin is not None:
            cfg.cp = pin
        return ICcell


    #TODO: add bend_solve()
    def bend_p2l(self, radius=None, angle=None, width=None, pin=None, xs=None,
            length1=0, ref=None, name=None, arrow=True, offset=None,
            max_length=None):
        """Create a bent waveguide (circular arc) with optional strt input and ending at a line.

        Args:
            radius (float): radius at the center line of the arc in um
            width (float): width of the arc in um
            angle (float): angle of arc in degree (default=90)
            pin (Node): optional Node for modeling info
            xs (str): optiinal xsection of bend
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            offset (float): optional new offset for this bend only.
            length1 (float): length of an optional straight section before the bend
            length2 (float): length of an optional straight section after the bend
            ref (Node | Instance)| tuple(x, y, a)): the reference line to intersect

        Returns:
            Cell: circularly-bent waveguide element

        Example:
            Create and place a bend::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.bend(angle=45)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        radius = self._getradius(pin, radius, xs)
        pinflip = not nd.get_xsection(xs).symmetry

        pinb, T = nd.parse_pin(pin)
        pin = pinb.move(*T)
        refb, T = nd.parse_pin(ref)
        ref = refb.move(*T)
        ref2b = nd.diff(pin, ref)
        ref2, T = nd.parse_pin(ref2b)

        if offset is None:
            offset = self.offset

        if angle is None:
            angle = self.angle
        if name is None:
            name = 'ic_bend'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            refrel = nd.Pin().put(*T)
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self.line(length=length1, width=width, xs=xs).put(0)
            e2 = self._arc(radius=radius, width=width, angle=angle, xs=xs,
                offset=offset).put()
            self.strt_p2l(pin=e2.pin['b0'], ref=refrel, width=width, xs=xs,
                name=None, arrow=False, max_length=max_length).put()
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
        if pin is not None:
            cfg.cp = pin
        return ICcell


    #TODO: add ptaper_solve()
    def ptaper(self, length=None, width1=None, width2=None, pin=None, xs=None,
            name=None, arrow=True):
        """Create a parabolic taper.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            pin (Node): optional Node for modeling info
            xs (str): optional xsection of taper
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: parabolic taper element

        Example:
            Create and place a parabolic taper::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.ptaper(length=10, width1=2.0, width2=5.0)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width1, xs, adapt=False)
        width2 = self._getwidth(pin, width2, xs, adapt=False)
        pinflip = not nd.get_xsection(xs).symmetry

        if length is None:
            length = self.length

        if name is None:
            name = 'ic_ptaper'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self._ptaper(length=length, width1=width1, width2=width2,
                name=name, xs=xs).put(0)
            nd.Pin('a0', io=0, width=width1, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width2, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()

        if pin is not None:
            cfg.cp = pin
        return ICcell


    #TODO: add taper_solve()
    def taper(self, length=None, width1=None, width2=None, shift=0, xs=None, pin=None,
            name=None, arrow=True):
        """Create a linear taper.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            shift (float): lateral shift of taper end
            xs (str): optional xsection of taper
            pin (Node): optional Node for modeling info
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: linear taper element

        Example:
            Create and place a linear taper::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.taper(length=10, width1=2.0, width2=5.0)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width1, xs, adapt=False)
        width2 = self._getwidth(pin, width2, xs, adapt=False)
        pinflip = not nd.get_xsection(xs).symmetry

        if length is None:
            length = self.length

        with nd.Cell('taper_{}'.format(xs), instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self._taper(length=length, width1=width1, width2=width2,
                shift=shift, xs=xs).put(0)
            nd.Pin('a0', io=0, width=width1, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width2, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def _strt_p2l_solve():
        """To be implemented."""
        pass


    #TODO: p2l_solve
    def strt_p2l(self, pin=None, ref=None, width=None, xs=None,
            name=None, arrow=True, max_length=None):
        """Create a straight guide to intersect a reference line.

        p2l: point-to-line. Note there is no solution for a reference line
        parallel to the pointer in pin. To avoid huge (near parallel) lines,
        a max-length can be specified.

        Args:
            pin (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            ref (Node | Instance)| tuple(x, y, a)): the reference line to intersect
            width (float): width of the interconnect in um
            xs (str): optional xsection of the strt
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            max_length (float): maximum length of the guide

         Example:
            Create and place a straigth waveguide to intersect with a
            reference line::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.bend(angle=45)
                guide.put()
                nd.export_plt()

        Returns:
            Cell: straight waveguide element
        """
        if max_length is None:
            max_length = self.max_length
        if pin is None:
            pin = cp.here()
        if ref is None:
            cfg.interconnect_errcnt += 1
            msg = "IC-{} strt_p2l needs a reference line via ref= keyword.".format(cfg.interconnect_errcnt)
            interconnect_logger(msg, 'error')
            ic_exception(msg)
            xs = 'error'

        pinb, T = nd.parse_pin(pin)
        pin = pinb.move(*T)
        refb, T = nd.parse_pin(ref)
        ref = refb.move(*T)

        if xs != 'error':
            xs = self._getxs(pin, xs)
            pinflip = not nd.get_xsection(xs).symmetry
        else:
            pinflip = False
        width = self._getwidth(pin, width, xs)

        x, y, a = nd.diff(pin, ref)
        L = x + y * m.tan(m.radians(a-90))

        rot = 0
        if abs(L) > max_length:
            L = sign(L) * max_length
            msg = "Solution for strt_p2l too large: >{}.".format(max_length)
            interconnect_logger(msg, 'warning')
            ic_exception(msg)
        if L < 0:
            rot = 180
            L = -L
            xs = 'error'
            cfg.interconnect_errcnt += 1
            msg = "IC-{} negative length for strt_p2l.".format(cfg.interconnect_errcnt)
            interconnect_logger(msg, 'error')
            ic_exception(msg)

        if name is None:
            name = 'ic_strt_p2l'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self.line(length=L, width=width, xs=xs).put(0)
            nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'].rot(rot))
            nd.Pin('b0', io=1, width=width, xs=xs).put()
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
        cfg.cp = pin
        return ICcell


    def _strt_p2p_solve(self, pin1=None, pin2=None, width=None, xs=None):
        """Find strt point-to-point solution.

        Returns:
            properties, parameters"""
        parse = self._p2p_parse(pin1=pin1, pin2=pin2, xs=xs, width1=width)
        pin1, pin2, xs, width1, width2, radius1, radius2 = parse

        x, y, a = nd.diff(pin1, pin2)
        length = m.hypot(x, y)
        b = m.degrees(m.atan2(y, x))
        #print(x, y, length, b)
        return parse, (length, b)


    def strt_p2p(self, pin1=None, pin2=None, width=None, xs=None, name=None,
            arrow=True):
        """Create point-to-point straight interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start of waveguide
            pin2 (Node | Instance | tuple(x, y, a)): end of waveguide
            width (float): width of waveguide
            xs (str): optional xsection of waveguide
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: straight waveguide element

        Example:
            Create and place a straight guide between points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.strt_p2p(pin1=(0, 0), pin2=(10, 10))
                guide.put()
                nd.export_plt()
        """
        parse, (length, b) = self._strt_p2p_solve(pin1=pin1, pin2=pin2,
            width=width, xs=xs)
        pin1, pin2, xs, width, _, radius1, radius2 = parse
        if name is None:
            name = 'ic_strt_p2p'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self.line(length=length, width=width, xs=xs).put(0, 0, b)
            p1 = nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'].rot(-b))
            p2 = nd.Pin('b0', io=1, width=width, xs=xs).put(e1.pin['b0'].rot(-b-pin1.a+pin2.a+180))
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            ICcell.pin2 = pin2
            # connect_opt needed because the extra rot in p1 and p2 does hide the line optical connection:
            nd.connect_path(p1, p2, trace.trace_length())
        cfg.cp = pin1
        return ICcell

    def taper_p2p(
        self,
        pin1=None,
        pin2=None,
        width1=None,
        width2=None,
        xs=None,
        name=None,
        arrow=True
    ):
        """Create point-to-point (angled) taper interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            width1 (float): width at start (taken from pin1 if None)
            width2 (float): width at end (taken from pin2 if None)
            xs (str): optional xsection of taper
            name (str): optional new name for the component
            BendEndFlag (int): (default=1)
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: taper element

        Example:
            Create and place a taper to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2, radius=10)

                guide = ic.taper_p2p(pin1=(0), pin2=(40, 10), width2=4)
                guide.put()
                nd.export_plt()
        """
        save_adapt = self.adapt_width
        self.adapt_width = True  # Adapt taper width to what it connects to.
        parse = self._p2p_parse(
            pin1=pin1,
            pin2=pin2,
            width1=width1,
            width2=width2,
            xs=xs
        )
        pin1, pin2, xs, width1, width2, radius1, radius2 = parse
        length, shift, _ = nd.diff(pin1, pin2)
        if name is None:
            name = 'ic_taper_p2p'
        with nd.Cell(f'{name}_{xs}', instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self._taper(
                length=length,
                width1=width1,
                width2=width2,
                shift=shift,
                xs=xs
            ).put(0, 0, pin1.xya()[2])
            p1 = nd.Pin('a0', io=0, width=width1, xs=xs).put(e1.pin['a0'])
            p2 = nd.Pin('b0', io=1, width=width2, xs=xs).put(e1.pin['b0'])
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            ICcell.pin2 = pin2
            # connect_opt needed because the extra rot in p1 and p2 does hide the line
            # optical connection:
            nd.connect_path(p1, p2, trace.trace_length())
        cfg.cp = pin1
        self.adapt_width = save_adapt  # restore previous value
        return ICcell

    def rot2ref_solve(self, pin=None, ref=None, angle=0, cw=None,
            length1=0, length2=0, width=None, radius=None, xs=None):
        """Calculate and return the angle to rotate from <pin> to reference direction <ref>.

        Note that only the angle part of <ref> is used in the calcuation.

        Args:
            pin (Node): starting pin (default=cp)
            ref (Node): reference pin (default=org)
            angle (float): rotation with repect to ref in [Degrees] (default=0)
            cw (bool): angle direction clockwise or counter clockwise (default is shortest bend)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            tuple, float: parse-info, angle to rotate from <pin> to <ref> + <a>
        """
        nd.cp.push()

        # TODO: handle ref as tuple
        if ref is None:
            ref = cfg.cells[-1].pin['org']
        else:
            ref = self._getpinout(ref)

        parse = self._p2p_parse(pin1=pin, xs=xs, width1=width, radius1=radius)
        pin1, pin2, xs, width, _, radius1, radius2 = parse

        if pin2 is None:
            raise Exception('Source pin not specified in rot2ref.')

        x, y, a = nd.diff(ref.rot(angle), pin2)
        if a >= 180:
            a -= 360

        if cw is True:
            if a < 0:
                a += 360
        elif cw is False:
            if a > 0:
                a -= 360

        nd.cp.pop()

        xya = (length1 - radius1 * m.sin(m.radians(a)), radius1 * (1 - m.cos(m.radians(a))), a)
        result = {
            'geo': [
                ('s', (length1, width)),
                ('b', (-m.degrees(a), radius1, width)),
                ('s', (length2, width))],
            'xya': xya,
            'params': {
                'lenght1': length1,
                'lenght2': length2,
                'angle': -a}
            }
        return parse, result


    def rot2ref(self, pin=None, ref=None, angle=0, length1=0, length2=0,
            cw=None, width=None, xs=None, radius=None, name=None, arrow=True):
        """Rotate a waveguide from <pin> to a specific <angle> with respect to reference <ref>.

        Args:
            pin (Node): starting pin (default=cp)
            ref (Node): reference pin (default=org)
            angle (float): rotation with repect to ref (default=0)
            length1 (float): length of straight waveguide starting section (default=0)
            length2 (float): length of straight waveguide ending section (default=0)
            cw (bool): angle direction clockwise or counter clockwise (default is shortest bend)
            width (float): width of waveguide
            xs (str): optional xsection of waveguide
            radius (float): radius at the center line of the bend in um
            name (str): optional new name for the component (default=rot2ref)
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: waveguide element rotating from <pin> into to the desired direction

        Example:
            Rotate to angle 125 degree w.r.t. 'org' after a straight guide of 100 um::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.rot2ref(angle=125, length1=10)
                guide.put()
                nd.export_plt()
        """
        parse, result = self.rot2ref_solve(
            pin=pin, ref=ref, angle=angle, cw=cw,
            length1=length1, length2=length2, width=width, radius=radius, xs=xs)
        pin1, pin2, xs, width, _, radius1, _ = parse
        pinflip = not nd.get_xsection(xs).symmetry
        a = result['params']['angle']

        if name is None:
            name = 'ic_rot2ref'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=self.instantiate, cnt=True) as ICcell:
            trace.trace_start()
            ICcell.group = 'interconnect'
            e1 = self.line(length=length1, width=width, xs=xs).put(0)
            self._arc(angle=a, radius=radius1, width=width, xs=xs).put()
            self.line(length=length2, width=width, xs=xs).put()
            nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, xs=xs).put()
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)

        if pin2 is not None:
            cfg.cp = pin2
        return ICcell


    def _sbend_solve(self, radius=None, width=None, pin=None, xs=None, offset=20,
            Ltot=None, length1=None, length2=None, Amax=90.0):
        """Calculate sbend solution.

        Returns:
            tuple, dict: parse-info, solution-parameters
        """
        if Ltot is not None:
            if length1 is not None and length2 is not None:
                raise Exception("Can not use Ltot together with length1 and length2.")
        else:
            if length1 is None:
                length1 = 0
            if length2 is None:
                length2 = 0
            Ltot = 0

        xs = self._getxs(pin, xs)
        width  = self._getwidth(pin, width, xs)
        radius  = self._getradius(pin, radius, xs)

        L, La, Lb = 0, 0, 0
        Amax = m.radians(Amax)

        if 2*radius*(1-m.cos(Amax)) > abs(offset):
            A = m.acos(1-abs(offset)/(2*radius))
            if offset > 0:
                A = -A
        else:
            L = (abs(offset) - abs(2*radius*(1-m.cos(Amax))))/m.sin(Amax)
            A = sign(-offset)*Amax

        Lx = 2*radius * m.sin(abs(A)) + L * m.cos(abs(A))
        dLx = abs(Ltot)-Lx

        if Ltot == 0:
            La = length1
            Lb = length2
        else:
            if length1 is not None:
                if La > dLx:
                    raise Exception("Ltot too short for length1={}".
                        format(length1))
                else:
                    La = length1
                    Lb = dLx - La
                    Ltot = abs(Ltot)
            elif length2 is not None:
                if Lb > dLx:
                    raise Exception("Ltot too short for length2={}".
                        format(length2))
                else:
                    Lb = length2
                    La = dLx - Lb
                    Ltot = abs(Ltot)
            elif Ltot < 0 and dLx > 0:
                La = dLx
                Lb = 0
            elif Ltot > 0 and dLx > 0:
                La = 0
                Lb = dLx
        parse = xs, width, radius

        xya = (La+Lb+2*abs(radius*m.sin(A))+L*abs(radius*(1-m.cos(A)))), offset, 0

        result = {
            'geo': [
                ('s', (La, width)),
                ('b', (-m.degrees(A), radius, width)),
                ('s', (L, width)),
                ('b', (m.degrees(A), radius, width)),
                ('s', (Lb, width))],
            'xya': xya,
            'params': {
                'La': La,
                'L': L,
                'Lb': Lb,
                'A': A}
            }

        return parse, result


    def sbend(self, radius=None, width=None, pin=None, xs=None, offset=20,
            Ltot=None, length1=None, length2=None, name=None, arrow=True, Amax=90.0):
        """Create an s-bend interconnect.

        Args:
            radius (float): bend radius at the center line of the arc in um
            width (float): width of the interconnect in um
            pin (Node): optional Node for modeling info
            xs (str): xsection of sbend
            offset (float): lateral offset of the sbend in um
            Ltot (float): optional total forward length of the sbend in um.
                When positive, additional length is added at the
                end of the s-bend, when positive it is added at the end,
                all provided the forward length of the
                s-bend itself is shorter than abs(Ltot).
            length1 (float):
            length2 (float):
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            Amax: maximum angle of bends. default is 90

        Returns:
            Cell: sbend element

        Example:
            Create and place a sbend waveguide::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.sbend(offset=20)
                guide.put()
                nd.export_plt()
        """
        parse, result = self._sbend_solve(xs=xs, width=width, radius=radius,
            pin=pin, offset=offset, Ltot=Ltot, length1=length1, length2=length2, Amax=Amax)
        xs, width, radius = parse
        pinflip = not nd.get_xsection(xs).symmetry

        if name is None:
            name = 'ic_sbend'
        C = self.tube(geo=result['geo'], name='{}_{}'.format(name, xs), xs=xs, arrow=arrow)
        return C




    #TODO: new tube option?
    def sinebend(self, width=None, pin=None, xs=None, distance=200, offset=20,
            name=None, arrow=True):
        """Create a sine-bend interconnect.

        This interconnect has zero curvature at both ends and can be
        connected without offsets to straight waveguides.

        Args:
            width (float): width of the interconnect in um
            pin (Node): optional Node for modeling info
            xs (str): xsection of sinebend
            distance (float): total forward length of the sinebend in um
            offset (float): lateral offset of the sinebend in um
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: sinebend element

        Example:
            Create and place a sinebend waveguide::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.sinebend(distance=100, offset=50)
                guide.put()
                nd.export_plt()
        """
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        pinflip = not nd.get_xsection(xs).symmetry

        if name is None:
            name = 'ic_sinebend'
        with nd.Cell('{}_{}_{}_{}'.format(name, xs, int(distance), int(offset)),
                instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            if abs(offset) < 1e-6:
                # Straight line will do just fine.
                e1 = self.strt(length=distance).put(0)
            else:
                e1 = self.sinecurve(width=width, distance=distance, offset=offset,
                    xs=xs).put(0)
            nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, xs=xs).put(distance, offset, 0)
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)

        if pin is not None:
            cfg.cp = pin
        return ICcell

    # TODO: new tube option?

    def _sbend_p2p_solve(self, pin1=None, pin2=None, width=None, radius=None,
            xs=None, length1=0, doStrFirst=1, ortho=True, ref=None, Amax=90.0):
        """Calculate an sbend_p2p solution.

        Args:
            pin1 (Node): start point
            pin2 (Node): end point
            width (float): optional width to overrule xs
            radius (float): optional radius to overrule xs
            xs (str): xsection
            length1 (float):
            doStrFirst (bool)
            ortho (bool): if True, ensure smooth sbend connections by adding
                start and end bends (default=True)
            ref (Node | xya | a): optional direction of the sbend io (and straight section).
                Example xya = (0, 0, 20), example a = 20.
                Use with ortho=True.
            Amax: maximum angle of the bends

        Returns:
            tuple, dict: parse-info, solution-parameters
        """
        parse = self._p2p_parse(pin1, pin2, xs, width1=width, radius1=radius)
        pin1, pin2, xs, width, _, radius1, radius2 = parse

        message = ''
        found = True
        A = 0
        Ltap = 0
        Amax = m.radians(Amax)

        # rotate to sbend-ref
        geo1 = []
        geo2 = []
        da1 = da2 = 0
        if ortho and ref is None:
            ref = pin1
        if ref is not None:
            # make a pin out of ref:
            if isinstance(ref, (int, float)):
                ref = (0, 0, ref)
            ref, T = nd.parse_pin(ref)
            ref =  ref.move(*T)

            if ref is not pin1:
                dx1, dy1, da1 = nd.diff(pin1, ref)
                da1 = zeroDeg(da1)
                pin1 = pin1.move(abs(radius1*m.sin(m.radians(da1))),
                    sign(da1)*radius1*(1-m.cos(m.radians(da1))), da1)
                if da1 != 0:
                    geo1 = [('b', (da1, radius1, width))]
            if ref is not pin2:
                dx2, dy2, da2 = nd.diff(pin2, ref.rot(180))
                da2 = zeroDeg(da2)
                pin2 = pin2.move(abs(radius1*m.sin(m.radians(da2))),
                    sign(da2)*radius1*(1-m.cos(m.radians(da2))), da2)
                if da2 != 0:
                    geo2 = [('b', (-da2, radius1, width))]

        xya = nd.diff(pin1, pin2)
        dx, dy, da = xya

        if length1 < 0:
            length1 = abs(length1)
            doStrFirst = 0

        # get dy from start and end position
        if length1 < Ltap:
            length1 = Ltap

        H = 2*radius1*(1-m.cos(Amax))+2*Ltap*m.sin(Amax)
        if H > abs(dy):
            # not enough offset for Amax--> no vertical straight guide section
            if Ltap == 0:
                A = m.acos(1-abs(dy)/(2*radius1))
            else:
                tel=0
                A = 0
                da=0.2
                damin=1e-8
                while abs(da) > damin and tel < 100:
                    if Ltap*m.sin(A)-radius1*m.cos(A) < (abs(dy)-2*radius1)/2.0:
                        A += da
                    else:
                        A -= da
                        da /= 2.0
                        A += da
                    tel += 1
                if A < 10*damin:
                    A=0

            A = sign(dy)*A
            Lstr = 0
        else: # use Amax angle
            A = sign(dy)*Amax
            Lstr = (abs(dy) - abs(2*radius1*(1-m.cos(Amax))) -2*Ltap)/abs(m.sin(Amax))

        Lfit = dx -2*radius1*abs(m.sin(A)) - (Lstr+2*Ltap)*abs(m.cos(A)) - length1
        La, Lb = 0, 0
        if doStrFirst == 1:
            La = length1-Ltap
            Lb = Lfit-Ltap
        else:
            Lb = length1-Ltap
            La = Lfit-Ltap

        if La < 0 or Lb < 0:
            found = False
            message = "(interconnect): No solution for sbend_p2p."

        else:
            if A == 0:
                Lstr += 4*Ltap
        geo = geo1 + [
            ('s', (La, width)),
            ('b', (m.degrees(A), radius1, width)),
            ('s', (Lstr, width)),
            ('b', (-m.degrees(A), radius1, width)),
            ('s', (Lb, width))] + geo2
        params = {
            'angle1': da1,
            'angle2': -da2,
            'angle': A,
            'length1': La,
            'length2': Lstr,
            'length3': Lb,
            'Lfit': Lfit,  # length left to reach end point.
            'message': message,
            'Amax': Amax,
            'Ltap': Ltap}
        result = {
            'found': found,
            'message': message,
            'solution': params,
            'geo': geo,
            'xya': xya}
        return parse, result

    def sbend_p2p(self, pin1=None, pin2=None, width=None, radius=None, Amax=90,
                xs=None, doStrFirst=1, Lstart=0, BendEndFlag=1, ref=None,
                name=None, arrow=True, bsb=True):
        """Create point-to-point s-bend interconnect.

        The direction of the end pin is ignored.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            width (float): width of the interconnect in um
            radius (float): bend radius of the interconnect in um
            xs (str): optional xsection of sbend
            doFirst (int): (default=1)
            Amax (float): maximum bend angle (default=90)
            Lstart (float): straight waveguide length at beginning (positive
                value) or end (negative value) of sbend
            ref (Node): reference direction for the sbend (default=pin1).
            name (str): optional new name for the component
            BendEndFlag (int): (default=1)
            arrow (bool): draw connection arrows (default=True)
            bsb (bool): If True, use bend_straight_bend_p2p() as fallback (default=True)

        Returns:
            Cell: sbend element

        Example:
            Create and place a sbend to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.sbend_p2p(pin1=(0), pin2=(40, 20))
                guide.put()
                nd.export_plt()
        """
        parse, params = self._sbend_p2p_solve(
            pin1=pin1, pin2=pin2, width=width, radius=radius, length1=Lstart,
            doStrFirst=doStrFirst, ref=ref, Amax=Amax)
        pin1, pin2, xs, width, _, radius1, radius2 = parse
        pinflip = not nd.get_xsection(xs).symmetry

        if not params['found'] or abs(params['solution']['angle']) < 0*m.pi:
            interconnect_logger(params['message'], 'info')
            #TODO: alterative back up if La is so long that Lb is negative?
            if bsb:
                msg = "replacing sbend with bend_strt_bend."
                interconnect_logger(msg, 'info')
                return self.bend_strt_bend_p2p(
                    pin1=pin1, pin2=pin2, radius=radius, width=width, xs=xs,
                    name=name, arrow=arrow)
            else:
                return self.strt_p2p(pin1, pin2, xs='error')

        else:
            if name is None:
                name = 'ic_sbend'
            cfg.cp = pin1
            C = self.tube(geo=params['geo'], name='{}_{}'.format(name, xs), xs=xs, arrow=arrow)
            C.pin2 = pin2
            return C


    def _bend_strt_bend_p2p_solve(self, pin1=None, pin2=None, radius1=None,
             radius2=None, xs=None, width=None, ictype='shortest',
             length1=0, length2=0):
        """Calculate geometry for a bend_strt_bend interconnect.

        Returns:
            tuple, dict: parse-info, solution-parameters
        """
        parse = self._p2p_parse(pin1, pin2, xs=xs, width1=width,
            radius1=radius1, radius2=radius2)
        pin1, pin2, xs, width, _, radius1, radius2 = parse

        Ltap = 0
        A =  pin1.move(Ltap+length1, 0, 0) # to calculate the geometry with Ltap
        B =  pin2.move(Ltap+length2, 0, 180)
        xya = nd.diff(A, B)
        dx, dy, da = xya

        radius1 -= 1e-8
        radius2 -= 1e-8
        # Calculate circle centers. Note that pin1 is put at (0,0,0)
        c1Lx, c1Ly = 0, radius1
        c1Rx, c1Ry = 0, -radius1
        c2Lx, c2Ly = dx-radius2*m.sin(m.radians(da)), dy+radius2*m.cos(m.radians(da))
        c2Rx, c2Ry = dx+radius2*m.sin(m.radians(da)), dy-radius2*m.cos(m.radians(da))

        message = ''
        solutions = {}
        if ictype in ['rr', 'rl', 'lr' 'll']:
            shapes = [ictype]
        else:
            shapes = ['rr', 'rl', 'lr', 'll']

        options = ['rr', 'rl', 'lr', 'll', 'shortest', 'all']
        if ictype not in options:
            msg = "ictype '{}' not defined in cell '{}', switching to 'shortest'.".\
                format(ictype, cfg.cells[-1].cell_name)
            interconnect_logger(msg, 'warning')
            interconnect_logger("Valid options are {}.".format(options), 'warning')
            ictype = 'shortest'
            ic_exception(msg)

        for shape in shapes:
            if shape == 'rr':
                sx, sy = c2Rx-c1Rx, c2Ry-c1Ry
                d1, d2 = 1, -1
            if shape == 'rl':
                sx, sy = c2Lx-c1Rx, c2Ly-c1Ry
                d1, d2 = 1, 1
            if shape == 'lr':
                sx, sy = c2Rx-c1Rx, c2Ry-c1Ly
                d1, d2 = -1, -1
            if shape == 'll':
                sx, sy = c2Lx-c1Lx, c2Ly-c1Ly
                d1, d2 = -1, 1

            rr = d1*radius1 + d2*radius2
            s = m.sqrt(sx**2+sy**2)  # (sx,sy): (x,y)-distance between circle centers
            if rr > s+1e-6:
                result = {
                    'found': False,
                    'message': "Radii too large (points to close) in bend_strt_bend in cell '{}' for shape '{}'. Maximum radius={:0.3f} or radius1+radius2<{:0.3f}. (radius1={:0.3f}, radius2={:0.3f})".\
                        format(cfg.cells[-1].cell_name, shape, 0.5*s, s, radius1, radius2)}
                #if len(shapes) == 1:
                interconnect_logger(result['message'], 'warning')
                found = False
                #return parse, result #found, 0, 0, 0, 0
            else:
                gs = m.atan2(sy, sx) # angle through the circle centres at the start
                if abs(rr/s) <= 1:
                    found = True
                    gb = m.asin(rr/s) # angle through the circle centers after placing connecting straight horizontal
                else:
                    found = False
                    result = {
                        'found': False,
                        'message': "No solution Found in bend_strt_bend in cell '{}'.".\
                            format(cfg.cells[-1].cell_name)}
                    interconnect_logger(result['message'], 'warning')
                    #return parse, result #found, 0, 0, 0, 0

                gt = -gs+gb # angle of rotation of axis through circle centers from to put connection between circles horizontal.
                t1 = m.radians(0)+gt # angle of spoke that points to start-point bsb on the circle
                t2 = m.radians(da)+gt # angle of spoke that points to end-point bsb on the circle
                if d1 == 1:
                    b = negRad(-t1) # angle of drawn self._arc on start self._arc
                else:
                    b = posRad(-t1) # angle of drawn self._arc on start self._arc

                if d2 == -1:
                    e = negRad(t2) # angle of drawn self._arc on end self._arc
                else:
                    e = posRad(t2) # angle of drawn self._arc on start self._arc

                L = s*m.cos(gb) # length of straight self.line
                Ltot = L + radius1*abs(b)+radius2*abs(e) # total connection length

                if b == 0:
                    L += Ltap
                else:
                    L -= Ltap
                if e == 0:
                    L += Ltap
                else:
                    L -= Ltap
                if L < 0:
                    found = True

            solutions[shape] = {'solution': {
                'found': found,
                'Ltot': Ltot,
                'L': L,
                'angle1': b,
                'angle2': e}}
            solutions[shape]['geo'] = [
                ('b', (m.degrees(b), radius1, width)),
                ('s', (L, width)),
                ('b', (m.degrees(e), radius2, width))]

        if ictype == 'shortest':
            shortest = None
            LMin = 1e10
            for shape, info in solutions.items():
                #if info['solution']['found']:
                if info['solution']['Ltot'] < LMin and info['solution']['found']:
                    LMin = info['solution']['Ltot']
                    shortest = shape
            if shortest is not None:
                result = {
                    'found': True,
                    'type': shortest,
                    'solution': solutions[shortest]['solution'],
                    'geo': solutions[shortest]['geo'],
                    'xya': xya,
                    'message': message}
            else:
                raise Exception("No shortest solution found for bend_strt_bend.")
            return parse, result
        elif ictype == 'all':
            result = {'type': 'all'}
            if len(solutions) > 0:
                result['found'] = True
            curves = []
            for variation in solutions:
                curves.append({
                    'found': True,
                    'message': message,
                    'solution': solutions[variation]['solution'],
                    'geo': solutions[variation]['geo'],
                    'xya': xya,
                    'type': variation})
            result['variations'] = curves
            return parse, result
        elif ictype in ['rr', 'rl', 'lr', 'll']:
            result = {
                 'found': True,
                 'type': ictype,
                 'solution': solutions[ictype]['solution'],
                 'geo': solutions[ictype]['geo'],
                 'xya': xya,
                 'message': message}
            return parse, result
        else:
            result = {
                'found': False,
                'message': "No solution found in bend_strt_bend in cell '{}'.".\
                    format(cfg.cells[-1].cell_name)}
            return parse, result


    def bend_strt_bend_p2p(self, pin1=None, pin2=None, radius=None,
            radius1=None, radius2=None, width=None, xs=None,
            length1=0, length2=0,
            ictype='shortest', name=None, arrow=True):
        """Generate a point-to-point bend-straight-bend interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius1 (float): optional first bend radius in um
            radius2 (float): optional second bend radius im um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            ictype (str): interconnection type (default='shortest')
                options: 'shortest', 'll', 'lr', 'rl', rr', 'all'
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: bend_strt_bend element

        Example:
            Create and place a bend-straight-bend guide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.bend_strt_bend_p2p(pin1=(0, 0, 0), pin2=(40, 20, 90))
                guide.put()
                nd.export_plt()
        """
        if radius is not None:
            if radius1 is None:
                radius1 = radius
            if radius2 is None:
                radius2 = radius
        parse, curves = self._bend_strt_bend_p2p_solve(pin1, pin2,
            xs=xs, width=width, radius1=radius1, radius2=radius2,
            length1=length1, length2=length2, ictype=ictype)
        pin1, pin2, xs, width, _, radius1, radius2 = parse
        pinflip = not nd.get_xsection(xs).symmetry

        instantiate = self.instantiate
        if ictype == 'all':
            instantiate = True
            variations = curves['variations']
        else:
            variations = [curves]

        cells = []
        if not curves['found']:
            interconnect_logger(curves['message'], 'error')
            return self.strt_p2p(pin1, pin2, xs='error')

        else:
            for curve in variations:
                par   = curve['solution']
                shape = curve['type']
                param   = curve['solution']
                Ltot  = param['Ltot']
                L     = param['L']
                b     = param['angle1']
                e     = param['angle2']
                if name is None:
                    name = 'ic_bend_strt_bend_p2p'
                with nd.Cell("{}_{}".format(name, shape),
                        instantiate=instantiate, cnt=True) as ICcell:
                    ICcell.group = 'interconnect'
                    trace.trace_start()
                    if length1 > 0:
                        e1 = self.line(length1, width, xs=xs).put()
                        self._arc(radius1, angle=m.degrees(b), width=width, xs=xs).put()
                    else:
                        e1 = self._arc(radius1, angle=m.degrees(b), width=width, xs=xs).put()
                    self.line(L, width, xs=xs).put()
                    self._arc(radius2, angle=m.degrees(e), width=width, xs=xs).put()
                    if length2 > 0:
                        self.line(length2, width, xs=xs).put()
                    nd.Pin('a0', io=0, width=width, xs=xs).put(e1.pin['a0'])
                    nd.Pin('b0', io=1, width=width, xs=xs).put()
                    if arrow:
                        self.arrow.put(ICcell.pin['a0'])
                        self.arrow.put(ICcell.pin['b0'], flip=pinflip)
                    trace.trace_stop()
                    ICcell.length_geo = trace.trace_length()
                    ICcell.pin2 = pin2
                cfg.cp = pin1
                cells.append(ICcell)

        if not cells:
            msg = "No solution in bend_strt_bend_p2p."
            interconnect_logger(msg, 'warning')
            ic_exception(msg)

        elif len(cells) == 1:
            return cells[0]
        elif len(cells) > 1:
            with nd.Cell('{}_{}'.format(name, xs), cnt=True) as ICgroup:
                ICcell.group = 'interconnect'
                trace.trace_start()
                for cell in cells:
                    cell.put(180)
                trace.trace_stop()
                ICgroup.length_geo = trace.trace_length()
            cfg.cp = pin1
            return ICgroup


    def bend_strt_bend(self, pin=None, radius=None, radius1=None, radius2=None,
        width=None, xs=None, ictype='shortest', name=None, arrow=True):
        """Generate a bend-straight-bend connection starting at the current pointer.

        This is the same connection as 'bend_strt_bend_p2p' with pin1 = cp.
        """
        return self.bend_strt_bend_p2p(
            pin1=cp.here(), pin2=pin,
            radius=radius, radius1=radius1, radius2=radius2,
            width=width, xs=xs, ictype=ictype, name=name, arrow=arrow)


    def _strt_bend_strt_p2p_solve(self, pin1=None, pin2=None, radius=None, xs=None, width=None):
        """Solve geometry for a strt_bend_strt interconnect.

        Returns:
            tuple, dict: parse-info, solution-parameters
        """
        parse = self._p2p_parse(pin1, pin2, xs=xs, width1=width, radius1=radius)
        pin1, pin2, xs, width, _, radius1, _ = parse

        xya  = nd.diff(pin1, pin2.rot(180))
        dx, dy, da = xya
        if da >= 180:
            da -= 360
        g = (180 - da) / 2.0

        msg = ''
        L1 = 0
        L2 = 0

        if abs(180 - da) > 1e-5 and abs(180 + da) > 1e-5:
            x0 = dy / m.tan(m.radians(180 - da)) + dx
            found = True

            dx1 = abs(radius1 / m.tan(m.radians(g)))

            if sign(dy) * da < 0 or sign(dy) * da >= 180:
                found = False
                cfg.interconnect_errcnt += 1
                msg = "strt_bend_strt: Wrong direction of end-point. Switching to bend_strt_bend in cell '{}' between {} and {}".\
                    format(cfg.cells[-1].cell_name, pin1.fxya(), pin2.fxya())
            elif  abs(g - 90) < 1e-4:
                found = False
                cfg.interconnect_errcnt += 1
                msg = "strt_bend_strt: pointers in-line. Switching to bend_strt_bend in cell '{}' between {} and {}".\
                    format(cfg.cells[-1].cell_name, pin1.fxya(), pin2.fxya())
            else:
                Ltap = 0
                s = m.sqrt((x0 - dx)**2 + dy**2)
                L1 = x0 - dx1
                L2 = s - dx1
                L1 -= Ltap
                L2 -= Ltap
                L1 = round(L1, 6)
                L2 = round(L2, 6)
                if L1 >= 0 and L2 >= 0:
                    found = True
                else:
                    found = False
                    cfg.interconnect_errcnt += 1
                    msg = "strt_bend_strt: negative straight line. Switching to bend_strt_bend in cell '{}' between {} and {}".\
                    format(cfg.cells[-1].cell_name, pin1.fxya(), pin2.fxya())
        else:
            found = False
            cfg.interconnect_errcnt += 1
            msg = "strt_bent_strt: angle not possible. Switching to bend_strt_bendin cell '{}' between {} and {}".\
                format(cfg.cells[-1].cell_name, pin1.fxya(), pin2.fxya())
        geo = [
            ('s', (L1, width)),
            ('b', (da, radius1, width)),
            ('s', (L2, width))]
        params = {
            'length1': L1,
            'length2': L2,
            'da': da,
            'message': msg
            }
        result = {
            'found': found,
            'solution': params,
            'geo': geo,
            'xya': xya,
            'message': msg}
        return parse, result


    def strt_bend_strt_p2p(self, pin1=None, pin2=None, radius=None, width=None,
            xs=None, name=None, arrow=True):
        """Create point-to-point straight-bend-straight interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius (float): optional bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: strt_bend_strt element

        Example:
            Create and place a straight-bend-straight guide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.strt_bend_strt_p2p(pin1=(0, 0, 0), pin2=(40, 20, 90))
                guide.put()
                nd.export_plt()
        """
        # For the calculation always rotate the coordinate system into start
        # coordinate in direction of positive x-axis.

        parse, params = self._strt_bend_strt_p2p_solve(pin1=pin1, pin2=pin2,
            xs=xs, width=width, radius=radius)
        pin1, pin2, xs, width, _, radius1, radius2 = parse
        pinflip = not nd.get_xsection(xs).symmetry

        if params['found']:
            solution = params['solution']
            L1 = solution['length1']
            L2 = solution['length2']
            da = solution['da']

            if name is None:
                name = 'ic_strt_bend_strt_p2p'
            cfg.cp = pin1
            C = self.tube(geo=params['geo'], name='{}_{}'.format(name, xs), xs=xs, arrow=arrow)
            C.pin2 = pin2
            return C

        else:
            interconnect_logger(params['message'], 'warning')
            ic_exception(params['message'])

        # goto bend_strt_bend:
        return self.bend_strt_bend_p2p(pin1, pin2, radius1=radius,
            radius2=radius, width=width, xs=xs, arrow=arrow, name=name)


    def _ubend_p2p_solve(self, pin1, pin2, length=0, xs=None, width=None,
            radius=None, balance=0, end_angle=False):
        """Calculate a ubend geometry between two pins.

        Returns:
            tuple, dict: parse-info, solution-parameters
        """
        epsilon = 1e-6
        parse = self._p2p_parse(pin1, pin2, xs, width1=width, radius1=radius)
        pin1, pin2, xs, width, _, radius1, radius2 = parse

        geo = []
        xya = nd.diff(pin1, pin2)
        dx, dy, da = xya
        if end_angle:
            if da > 180:
                da -= 360
            ddx = radius1 * abs(m.sin(m.radians(da)))
            ddy = radius1 * m.copysign(1-m.cos(m.radians(da)), -da)
            dx += ddx
            dy -= ddy
            p2 = pin2.move(ddx, ddy, -da)
            b_end = [('b', (da, radius1, width))]
        else:
            p2 = pin2.rot(-da)
            b_end = []
        da = 0

        if dx < 0:
            L2 = length - dx
            L1 = length
        else:
            L1 = length + dx
            L2 = length

        if abs(dy) < 2*radius1:
            sign = -np.sign(dy)
            sign = sign if sign != 0 else 1
            d = sign * (epsilon + radius1 - 0.5 * abs(dy))
            d1 = d * (1 + balance) # right swing
            d2 = d * (1 - balance) # left swing
        else:
            d1, d2 = 0, 0

        geo += [('s', (L1, width))]
        p1 = pin1.move(L1)
        p2 = p2.move(L2) # p2 is adjusted end pin for angle_end
        if abs(d1) > epsilon:
            d1 = np.sign(d1)*epsilon + d1
            sbend = self._sbend_solve(offset=d1, radius=radius, width=width, xs=xs)
            geo += sbend[1]['geo']
            p1 = p1.move(*sbend[1]['xya'])
        geo_sb_out = []
        if abs(d2) > epsilon:
            d2 = np.sign(d2)*epsilon + d2
            sbend = self._sbend_solve(offset=d2, radius=radius, width=width, xs=xs)
            geo_sb_out = sbend[1]['geo'][::-1]
            dxs, dys, das = sbend[1]['xya']
            p2 = p2.move(dxs, -dys, -das)
        bsb = self._bend_strt_bend_p2p_solve(
            pin1=p1, pin2=p2,
            radius1=radius1, radius2=radius1, width=width, xs=xs)
        geo += bsb[1]['geo'] + geo_sb_out + [('s', (L2, width))] + b_end

        params = {
            'length1': L1,
            'length2': L2,
            'offset1': d1,
            'offset2': d2
        }
        result = {
            'found': True,
            'message': '',
            'solution': params,
            'geo': geo,
            'xya': xya
        }
        return parse, result


    def ubend_p2p(self, pin1=None, pin2=None, radius=None, width=None, xs=None,
            length=0, name=None, arrow=True, balance=0, end_angle=False):
        """Create point-to-point u-bend interconnect.

        An extra straight length can be added to the ubend with <length>.

        If the sideways translation needed in the ubend is <2*radius, then
        the ubend automatically introduces a 'horseshoe' shape. The horseshoe
        can be made sidelobed by a <balance> parameter between -1 and 1, where
        0 results in a symmetric shape.

        The orientation of the output pin does not matter unless
        end_angle=True is set (default=False). If True an extra bend is
        introduced on pin2 to align its direction with pin1 before drawing an
        orthogonal ubend.


        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius (float): optional bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection of ubend
            length (float): extra straight section for longer ubend (default=0)
            balance (float): for a ubend <2*radius sidewyas, shift the horseshoe shape (default=0)
            end_angle (bool): Take pin2 angle into account when connecting if True (default=False)

        Returns:
            Cell: ubend element

        Example:
            Create and place a ubend to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.ubend_p2p(pin1=(0, 0, 0), pin2=(10, 20, 90), length=10)
                guide.put()
                nd.export_plt()
        """
        parse, params = self._ubend_p2p_solve(
            pin1,
            pin2,
            length=length,
            width=width,
            radius=radius,
            balance=balance,
            end_angle=end_angle
        )
        pin1, pin2, xs, width, _, radius1, radius2 = parse

        if name is None:
            name = 'ic_ubend_p2p'
        cfg.cp = pin1
        C = self.tube(geo=params['geo'], name='{}_{}'.format(name, xs), xs=xs, arrow=arrow)
        C.pin2 = pin2
        return C


    def ubend(self, pin=None, offset=20.0, radius=None, width=None, xs=None,
            length=0, name=None, arrow=True, balance=0, end_angle=False):
        """Create u-bend interconnect from an offset.

        An extra straight length can be added to the ubend with <length>.

        If the sideways translation needed in the ubend is <2*radius, then
        the ubend automatically introduces a 'horseshoe' shape. The horseshoe
        can be made sidelobed by a <balance> parameter between -1 and 1, where
        0 results in a symmetric shape.

        The orientation of the output pin does not matter unless
        end_angle=True is set (default=False). If True an extra bend is
        introduced on pin2 to align its direction with pin1 before drawing an
        orthogonal ubend.


        Args:
            pin (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            offset (float): offset of the u-bend
            radius (float): optional bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection of ubend
            length (float): extra straight section for longer ubend (default=0)
            balance (float): for a ubend <2*radius sidewyas, shift the horseshoe shape (default=0)
            end_angle (bool): Take pin2 angle into account when connecting if True (default=False)

        Returns:
            Cell: ubend element
        """

        if pin is None:
            pin = cp.here()
        pin2 = pin.move(0, offset, 0)
        return self.ubend_p2p(pin1=pin, pin2=pin2, radius=radius, width=width, xs=xs,
            length=length, name=name, arrow=arrow, balance=balance, end_angle=end_angle)


    print_warning = True
    def pcurve_p2p(self, pin1=None, pin2=None, width=None, radius1=0, radius2=0,
            offset1=None, offset2=None, xs=None, name=None, arrow=True):
        if self.print_warning:
            print("""WARNING: the function 'pcurve_p2p' is now obsolete. You can obtain the
same functionality by calling the function 'cobra_p2p' in a similar way.
And the new cobra_p2p() interconnect function has more functionality.
In your code, please do now replace:
> pcurve_p2p(pin1, pin2, width, radius1, radius2, offset1, offset2, xs, name, arrow)
with
> cobra_p2p(pin1, pin2, width1, radius1, radius2, offset1, offset2, xs, name, arrow)
(not all arguments need to be present).
For more information, check the documentation of 'cobra_p2p'.
""")
            self.print_warning = False
        # Call the cobra equivalent.
        return self.cobra_p2p(pin1=pin1, pin2=pin2, width1=width, width2=width,
                radius1=radius1, radius2=radius2, offset1=offset1,
                offset2=offset2, xs=xs, name=name, arrow=arrow)


    def cobra_p2p(self, pin1=None, pin2=None, width1=None, width2=None,
            radius1=0, radius2=0, offset1=None, offset2=None, xs=None,
            name=None, arrow=True):
        """Create point-to-point cobra interconnect with a smooth curve and
        option to change the width along the curve.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            width1 (float|function): optional waveguide width in um. This can
                be a function, w(t). In that case width2 is not used.
            width2 (float): optional waveguide width in um at end
            xs (str): optional xsection
            radius1 (float): radius at start of the cobra (default=0 -> inf)
            radius2 (float): radius at start of the cobra (default=0 -> inf)
            offset1 (float): lateral offset at pin1
            offset2 (float): lateral offset at pin2
            name (str): optional new cell name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: cobra interconnect element

        Example:
            Create and place a cobra waveguide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.cobra_p2p(pin1=(0, 0, 0), pin2=(40, 20, 90))
                guide.put()
                nd.export_plt()
        """
        if width1 is None and pin1 is not None:
            try:
                width1 = pin1.width
            except AttributeError:
                width1 = None
        if width2 is None and pin2 is not None:
            try:
                width2 = pin2.width
            except AttributeError:
                width2 = None
        # For the calculation always rotate coordinate system to point start
        # coordinate in positive x-axis.
        parse = self._p2p_parse(pin1, pin2, xs=xs, width1=width1,
            width2=width2, radius1=radius1, radius2=radius2)
        pin1, pin2, xs, width1, width2, radius1, radius2 = parse
        pinflip = not nd.get_xsection(xs).symmetry

        dx, dy, da = nd.diff(pin1, pin2.rotate(180))
        xya = (dx, dy, da)

        if name is None:
            name = 'ic_cobra_p2p'
        with nd.Cell('{}_{}_{}_{}_{}'.format(name, xs, int(dx), int(dy), int(da)),
                instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            if abs(dy) < 1e-6 and abs(da) < 1e-6 and \
                    abs(radius1) < 1e-6 and abs(radius2) < 1e-6:
                # Straight line (or taper) will do just fine.
                e1 = self.ptaper(length=dx, width1=width1, width2=width2).put(0)
                ICcell.Rmin = 0
            else:
                e1 = self.cobra(xya, width1=width1, width2=width2, radius1=radius1,
                    radius2=radius2, offset1=offset1, offset2=offset2,
                    xs=xs).put(0)
                ICcell.Rmin = e1.cell.properties['Rmin']
            nd.Pin('a0', io=0, width=width1, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width2, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            ICcell.pin2 = pin2
        cfg.cp = pin1
        ICcell.length_geo = trace.trace_length()
        return ICcell


    def euler(self, pin=None, width=None, width2=None,
            radius=None, angle=90, xs=None,
            name=None, arrow=True):
        """Create an Euler bend from a straight guide to a curvature of radius at angle.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            width1 (float|function): optional waveguide width in um. This can
                be a function, w(t). In that case width2 is not used.
            width2 (float): optional waveguide width in um at end
            angle (float): end angle
            radius (float): end radius
            xs (str): optional xsection
            radius (float): radius at start of the cobra (default=0 -> inf)
            name (str): optional new cell name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: Euler interconnect element

        Example:
            Create and place a cobra waveguide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.euler(angle=45)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        if width2 is None:
            width2 = width
        radius = self._getradius(pin, radius, xs)
        pinflip = not nd.get_xsection(xs).symmetry

        if name is None:
            name = 'ic_euler'
        with nd.Cell('{}_{}'.format(name, xs),
                instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()

            e1 = self.euler_base(width1=width, width2=width2, radius=radius,
                angle=angle, xs=xs).put(0)
            ICcell.Rmin = radius
            nd.Pin('a0', io=0, width=width, xs=xs, radius=0).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width2, xs=xs, radius=radius).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])  # TODO: nopinflip here?
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            #ICcell.pin2 = pin2
        #cfg.cp = pin
        ICcell.length_geo = trace.trace_length()
        return ICcell


    def euler2(
        self,
        pin=None,
        width=None,
        width2=None,
        radius=None,
        angle=45,
        xs=None,
        name=None,
        arrow=True,
    ):
        """Create a symmetric Euler bend with no curacture at the start and minimum curvature radius.
        """
        halfeuler = self.euler(
            pin=pin,
            width=None,
            width2=None,
            radius=radius,
            angle=angle/2.0,
            xs=None,
            name=None,
            arrow=False,
        )
        xs = self._getxs(pin, xs)
        pinflip = not nd.get_xsection(xs).symmetry
        if name is None:
            name = 'ic_euler2'
        with nd.Cell(
            '{}_{}'.format(name, xs),
            instantiate=self.instantiate,
            cnt=True
        ) as ICcell:
            ICcell.group = 'interconnect'
            e1 = halfeuler.put(0)
            e2 = halfeuler.put('b0', e1.pin['b0'], flip=True)
            ICcell.Rmin = radius
            nd.Pin('a0', io=0, width=width, xs=xs, radius=0).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, xs=xs, radius=0).put(e2.pin['a0'])
            if arrow:
                self.arrow.put(ICcell.pin['a0'])  # TODO: nopinflip here?
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
        ICcell.length_geo = 2 * halfeuler.length_geo
        return ICcell


    def mamba(self, points, radius=None, width=None, pin=None, xs=None,
            N=1, pitch=10, offset=0, polyline=True, showpins=False,
            name=None, arrow=True):
        """Create a snake-like interconnect guided by a list of points (x, y).

        Start and end of a mamba are set as pins 'a0' and 'b0', respectively.
        To put the Mamba in the layout in absolute cell coordinates, i.e.
        on the literal coordinates as provided in <points> use:
        mamba(...).put('org', 0).

        Args:
            points: list of (x, y) positions to guide the mamba
            radius (float): optional waveguide radius (default self.radius)
            width (float): optional waveguide width (default self.width)
            pin (Node): optional Node for modeling info
            xs (str): optional xsection of mamba
            N (int): number of parallel guides in the mamba
            pitch (float): pitch of the guide if N>1
            offset (float): lateral offset in the position of all guides
            polyline (bool): boolean determining if the mamba is also drawn as polyline
                (default=True)
            showpins (bool): show the points as dots in the layout (default=False)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)

        Returns:
            Cell: Mamba element based on the provided <points>

        Example:
            Create a Mamba and attach the first point to the current pin::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.mamba(points=[(10, 10), (20, 20), (20, 50), (10, 50)])
                guide.put(0) # put first mamba point 'a0' on a pin
                guide.put('org', 0) # put mamba 'org' in 0 for absolute coordinates
                nd.export_plt()

            Hence to put a mamba at absolute coordinates of <points> in the cell::

                guide.put('org', 0)
        """
        nd.cp.push()

        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        radius = self._getradius(pin, radius, xs)
        pinflip = not nd.get_xsection(xs).symmetry

        ring = nd.Polygon(points=geom.circle(radius=width/2), layer=self.layer)

        #create points along the mamba for interconnects:
        p1, p2 = [], []
        size = len(points)
        for i in range(size-1):
            dx = points[i+1][0] - points[i][0]
            dy = points[i+1][1] - points[i][1]
            a = np.degrees(m.atan2(dy, dx))
            p1.append((points[i][0], points[i][1], a))
            p2.append((points[i+1][0], points[i+1][1], a+180))

        #print('p1[0]:', p1[0])
        if name is None:
            name = 'ic_mamba'
        with nd.Cell('{}_{}'.format(name, xs), instantiate=False, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            ICcell.default_pins('pla0', 'plb0')
            start = nd.Pin(width=width, xs=xs).put(p1[0])
            nd.Pin('pla0', width=width, xs=xs).put(start.rot(180))

            #loop over guides
            for num, dis in enumerate([pitch*(n-0.5*(N-1))+offset for n in range(N)]):
                trace.trace_start()
                last = start.move(0, -dis)
                nd.Pin('a'+str(num), io=0, xs=xs).put(last.rot(180))
                if showpins:
                    ring.put(last)
                for i in range(size-2): #loop over points
                    if showpins:
                        pin0 = nd.Pin().put(p1[i+1])
                        ring.put(pin0.move(0, -dis))

                    pin1 = last
                    pin2 = nd.Pin(width=width, xs=xs).put(p2[i+1]).move(0, dis)

                    cp.push()
                    parse, params = self._strt_bend_strt_p2p_solve(pin1, pin2, radius)
                    found    = params['found']
                    solution = params['solution']
                    L1       = solution['length1']
                    L2       = solution['length2']
                    da       = solution['da']
                    message  = solution['message']
                    cp.pop()
                    if found is True:
                        self.strt(length=L1, pin=last).put(last)
                        self.bend(angle=da, radius=radius).put()
                        last = cp.here()
                        if i is size-3:
                            self.strt(length=L2, pin=last).put()
                    else:
                        if i < size-3:
                            self.bend(radius=radius, pin=last, angle=da).put()
                            last = cp.here()
                        else:
                            self.bend_strt_bend_p2p(pin1, pin2, radius=radius).put(last)

                    if showpins:
                        plast = nd.Pin().put(p2[-1])
                        ring.put(plast.move(0, -dis))

                nd.Pin('b'+str(num), type=1).put(cp.here())
                trace.trace_stop()
                ICcell.length_geo = trace.trace_length()
            end = nd.Pin(width=width, xs=xs).put(p2[-1])
            nd.Pin('plb0').put(end.rot(180))
            if arrow:
                self.arrow.put(ICcell.pin['pla0'])
                self.arrow.put(ICcell.pin['plb0'], flip=pinflip)

            if polyline is True:
                if pin is None:
                    nd.Polyline(points=points, width=2, layer=1111).put(0)
                else:
                    nd.Polyline(points=points, width=2, layer=1111).put(p1[0][0], p1[0][1])
        return ICcell


    def tube(self, geo, showpins=False, name=None, xs=None, arrow=True):
        """draw interconnect based on symbols.

        Returns:
            Cell: Tube element based on the provided elements
        """
        pinflip = not nd.get_xsection(xs).symmetry
        if name is None:
            name = 'ic_tube'
        with nd.Cell(name, instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            el_list=[]
            new_geo=[]
            for i, tube in enumerate(geo):
                if tube[0] == 's':
                    if abs(tube[1][0])<1.e-5:
                        continue
                    el_list.append(self.strt(length=tube[1][0], width=tube[1][1], arrow=False).put())
                    new_geo.append(tube)
                elif tube[0] == 'b':
                    if abs(tube[1][0])<1.e-5:
                        continue
                    el_list.append(self.bend(angle=tube[1][0], radius=tube[1][1], width=tube[1][2], arrow=False).put())
                    new_geo.append(tube)
                else:
                    cfg.interconnect_errcnt += 1
                    msg = "IC-{} Tube element not recognized {} in cell '{}'".\
                        format(cfg.interconnect_errcnt, tube, cfg.cells[-1].cell_name)
                    interconnect_logger(msg, 'error')
                    ic_exception(msg)
            if el_list:      
                nd.Pin('a0', io=0, width=el_list[0].pin['a0'].width, xs=xs).put(el_list[0].pin['a0'])
                nd.Pin('b0', io=1, width=el_list[-1].pin['b0'].width, xs=xs).put(el_list[-1].pin['b0'])
            else: # nothing to place. TODO: Handle in cell.put(), where empty cells can be skipped, e.g. with cell.void = True
                nd.Pin('a0', io=0, width=None, xs=xs).put(0, 0, 180)
                nd.Pin('b0', io=1, width=None, xs=xs).put(0, 0, 0)
         
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
        return ICcell


    def connect(self, length=10, width1=None, width2=None, ic=None, name=None):
        """Connect two different interconnects using tapered shapes.

        The xsections of the interconnects need to have the same order and
        number of guides in the same layers.

        Args:
            ic (Interconnect): mandatory interconnect to connect to
            length (float): lenght of taper
            width1 (float): width of interconnect 'self'
            width2 (float): width of interconnect 'ic'
            name (str): name of cell. Default='connect'

        Returns:
            Cell: taper between two different interconnects 'self' -> 'ic'
        """
        if width1 is None:
            width1 = self.width
        if width2 is None:
            width2 = ic.width
        if ic is None:
            raise Exception("Keyword 'ic' is mandatory.")
        LT1 = nd.get_xsection(self.xs).mask_layers.reset_index()
        LT2 = nd.get_xsection(ic.xs).mask_layers.reset_index()
        D1 = LT1.to_dict(orient='index')
        D2 = LT2.to_dict(orient='index')
        if D1.keys() != D2.keys():
            raise Exception("Not the same layer stacks. Can not connect interconnects.")
        if name is None:
            name = 'connect'
        with nd.Cell(name=name, instantiate=self.instantiate, cnt=True) as C:
            for k in D1.keys():
                al1 = D1[k]['leftedgefactor']
                bl1 = D1[k]['leftedgeoffset']
                al2 = D2[k]['leftedgefactor']
                bl2 = D2[k]['leftedgeoffset']
                ar1 = D1[k]['rightedgefactor']
                br1 = D1[k]['rightedgeoffset']
                ar2 = D2[k]['rightedgefactor']
                br2 = D2[k]['rightedgeoffset']
                lay = D2[k]['layer_name']
                if not D2[k]['polyline']:
                    nd.Polygon(points=[(0, al1*width1+bl1), (length, al2*width2+bl2),
                        (length, ar2*width2+br2), (0, ar1*width1+br1)], layer=lay).put()
                else:
                    wpoly = width1*(al1-ar1)+(bl1-br1)
                    centre = 0.5*(width1*(al1+ar1)+(bl1+br1))
                    centreline = [(0, centre), (length, centre)]
                    nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
            nd.Pin('a0', xs=self.xs).put(0, 0, 180)
            nd.Pin('b0', xs=ic.xs).put(length)
        return C


# TODO: add anglei and angleo as keyword option or even as a function of paramters.
    def Tp_viper(
        self, 
        x,
        y,
        w,
        width1=None, 
        width2=None,
        xs=None, 
        layer=None, 
        N=200, 
        epsilon=1e-6, 
        name='viper', 
        arrow=True,
        anglei=None,
        angleo=None,
        **kwargs,
    ):
        """Add a viper method to the interconnect's name space.
        
        Args:
            
        Returns:
            function: Viper
        """
        pin = None
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width1, xs)
        width2 = self._getwidth(pin, width2, xs)

        Tp = nd.Tp_viper(
            x=x, 
            y=y, 
            w=w, 
            width1=width1, 
            width2=width2,
            xs=xs, 
            layer=layer, 
            N=N,
            epsilon=epsilon, 
            name=name, 
            anglei=anglei,
            angleo=angleo,
            **kwargs,
        )

        #if arrow:
        #    self.arrow.put(Tp.pin['a0'])
        #    self.arrow.put(Tp.pin['b0'])

        return Tp
 
