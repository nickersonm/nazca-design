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
# 2017-2021 (c) Ronald Broeke
#-----------------------------------------------------------------------
# -*- coding: utf-8 -*-
"""
Nazca module for interconnecting guides.
"""
from __future__ import annotations
from functools import partial
import numpy as np
from numpy import sign, sin, cos, radians, sqrt
import math as m
from scipy.integrate import quad
import nazca as nd
import nazca.cfg as cfg
from nazca.logging import logger
import nazca.cp as cp
import nazca.geometries as geom
import nazca.bb_util as bbu
import nazca.trace as trace
from nazca.netlist import interconnect_logger, Node, Cell
from nazca.util import parabolicvarwidth, linvarwidth


cfg.interconnects = {} # store all interconnects for reference/documentation
cfg.interconnect_errcnt = 0  # cnt warnings and errors to locate and catch them better.
num = -1
MAMBALAYER = 1111

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


class RibIt():

    def __init__(self, pitch=None, width=None, gap=None):
        """"""
        if (pitch is None and width is None) or (gap is None and width is None):
            nd.main_logger("Incomplete data to create a rib.", "error")
        try:
            lenp = len(pitch)
        except:
            lenp = 1
        try:
            lenw = len(width)
        except:
            lenw = 1
        N = max(lenp, lenw)
        if N == 1:
            N = None

        if gap is not None:
            pitch = []
            p0 = 0
            pitch.append(p0)
            for w1, w2 in zip(width[0:N-1], width[1:N]):
                p0 += 0.5 * (w1 + w2) + gap
                pitch.append(p0)
        #self.pitch = pitch
        #self.width = width
        #self.N = N
        self._expand(pitch, width, N)
        return None

    def _expand(self, pitch, width, N):
        """Expand a constant width or pitch in a list for each guide.

        Returns:

        """
        if isinstance(pitch, (float, int)):
            pitch = [pitch]
        if isinstance(width, (float, int)):
            width = [width]
        N = max(len(pitch), len(width))
        if len(pitch) == 1:
            pitch = [pitch[0] * n for n in range(N)]
        if len(width) == 1:
            width *= N
        self.pitch = pitch
        self.width = width
        self.N = N
        return None

    def get(self):
        return {'pitch': self.pitch, 'width': self.width, 'N': self.N}

    def add(self, rib, offset=10):
        """
        Returns:
            RibIt
        """
        pitch = [ self.pitch[-1] + offset + p for p in rib.pitch ]
        return RibIt(pitch=self.pitch + pitch, width=self.width + rib.width)

    def reverse(self):
        """
        Returns:
            RibIt
        """
        pitchmax = self.pitch[-1]
        pitch = [ pitchmax - p for p in self.pitch[::-1] ]
        return RibIt(pitch=pitch, width=self.width[::-1])

    def __repr__(self):
        return f"RibIt(pitch={self.pitch}, width={self.width})"

    def __str__(self):
        return f"pitch={self.pitch}, width={self.width}, N={self.N}"


    def __getitem__(self, val):
        """Allows list slicing for ribbons."""
        _pitch = [p for p in self.pitch.__getitem__(val)]
        _width = [w for w in self.width.__getitem__(val)]
        return RibIt(pitch=_pitch, width=_width)

# Examples:
#   R1 = RibIt(pitch=[0, 10, 30], width=[1, 2, 3])
#   R2 = RibIt(pitch=[0, 50, 60], width=[4, 5, 6])
#   R1.add(R2)
#   print(R1)
#   print(R2)



        #return  'pitch': rib1['pitch'] + rib2b['pitch'],
        # 'width': rib1['width'] + rib2['width'],
        # 'N': rib1['N'] + rib2['N'],



# RB: Separate the handling of interconnect that are based on a mask_element from
# those that are a combination of mask_elements. The former can use the
# mask_element's without redefining those.
# Hence, composite interconnects <-> mask_elements interconnects


# RB: remove xs from component parameters?
# or will this mess with an "adaptable xs that adapts to the pins you connect to.


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
        angle=None,
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
        pitch=10,
        N=1,
        sectionangle=110,
        mambalayer=MAMBALAYER,
        euler_scale=None,
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
            radius (float): default radius in um. Also calibration radius for euler bends.
            width (float): default waveguide width im um
            angle (float): default angle of a bend and calibation angle for euler bends.
                Default to None. This case the angle is taken for the Xsection (if defined) or set to 90.0.
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
            pitch (float): pitch in um in ribbon mode, default = 10.
            N (int): number of interconnects in a ribbon (default=1)
            sectionangle (float): maximum bend angle in ribbons before dividing it in sections
            mambalayer (layer): gds layer for Mamba polyline
            varname (str): string for documentation purposes only, i.e. to store
               the interconnects variable name it is assigned to like
               ic = Interconnect(var="ic")).
            euler_scale (float): default scaling for Euler bends in the Interconnect.
                Default is obtained from xsection.euler_scale if exists,
                else from the default radius at the default angle. The interconnect euler scale
                can also be changed via class method euler_calibrate(). In addition, multiple Euler
                scales can be defined using class method add_scaled_euler().

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

        if layer is None and xs is None:
            self.xs = cfg.default_xs_name
        else:
            self.xs = xs
        if self.xs is None:
            interconnect_logger("Created an Interconnect object with layer {} "\
                "but no xsection.".format(layer), 'warning')
        XS = nd.get_xsection(self.xs)

        self.radius = radius
        if self.radius is None:
            if XS is not None:
                self.radius = getattr(XS, 'radius', cfg.default_xs_radius)
        if self.radius is None:
            self.radius = cfg.default_xs_radius

        self.angle = angle
        if self.angle is None:
            if XS is not None:
                self.angle = getattr(XS, 'angle', cfg.default_xs_angle)
        if self.angle is None:
            self.angle = cfg.default_xs_angle

        self.width = width
        if self.width is None:
            if XS is not None:
                self.width = getattr(XS, 'width', cfg.default_xs_width)
        if self.width is None:
            self.width = cfg.default_xs_width

        self.layer = None if layer is None else nd.get_layer(layer)
        self.offset = offset

        if pinstyle is None:  # then use xsection pin setting as default if defined
            if XS is not None:
                xsstyle = getattr(XS, 'pinstyle', None)
                pinstyle = xsstyle
        self.pinstyle = pinstyle
        self.arrow = bbu.make_pincell(style=self.pinstyle)
        self.ribbon0 = bbu.make_pincell(style='ribbon0')
        self.ribbonN = bbu.make_pincell(style='ribbonN')

        self.length = 10
        self.sectionangle = sectionangle
        self.instantiate = instantiate
        self.xya = (100, 100, 10)
        self._line = nd.Tp_straight(
            xs=self.xs,
            layer=self.layer,
        )
        self._sinecurve = nd.Tp_sinecurve(
            xs=self.xs,
            layer=self.layer
        )
        self._cobra = nd.Tp_cobra(
            xya=self.xya,
            width1=None,
            width2=None,
            radius1=0,
            radius2=0,
            offset1=None,
            offset2=None,
            xs=self.xs,
            layer=self.layer,
        )

        self.euler_scale = euler_scale
        if self.euler_scale is None:
            if XS is not None:
                self.euler_scale = getattr(XS, 'euler_scale', None)
        if self.euler_scale is None:
            self.euler_scale = nd.euler_scale(radius=self.radius, angle=self.angle)
        self._euler_base = nd.Tp_euler(
            width=None,
            width2=None,
            radius=self.radius,  # TODO: This should be scale based
            xs=self.xs,
            layer=self.layer,
            angle=self.angle,
            scale=self.euler_scale,
        )

        # store euler components for multiple scales
        self.scaled_euler = {}
        self.scaled_euler_arc = {}
        self.scaled_euler_arc_euler = {}

        self.pcb = PCB
        if PCB:
            self.Farc = nd.Tp_arc2(xs=self.xs, layer=self.layer)
        else:
            self.Farc = nd.Tp_arc(xs=self.xs, layer=self.layer)
        self.Farc_static = nd.Tp_arc(xs=self.xs, layer=self.layer)  # Keep a static version on case Farc gets updated
        self._ptaper = nd.Tp_ptaper(xs=self.xs, layer=self.layer)
        self._taper = nd.Tp_taper(xs=self.xs, layer=self.layer)
        self.adapt_width = adapt_width
        self.adapt_xs = adapt_xs
        self.max_length = 1e5 # maximum line length in p2l.
        self.length = 10
        self.pinflip = True

        self.gridpatch = 0.0  # Experimental: add taper the size fo gridpatch to strt.
        self.pitch = pitch
        self.N = N
        self.mambalayer = mambalayer

        self.ticks = None
        self.widths = None
        return None


    @property
    def PCB(self):
        return self.pcb

    @PCB.setter
    def PCB(self, val):
        self.pcb = val
        if val:
            self.Farc = nd.Tp_arc2(xs=self.xs, layer=self.layer)
        else:
            self.Farc = nd.Tp_arc(xs=self.xs, layer=self.layer)

    def copy(
            self,
            ic=None,
            xs=None,
            PCB=None,
            varname=None,
            width=None,
            radius=None,
            pitch=None,
            pinstyle=None,
            offset=None,
            doc=None,
            layer=None,
        ):
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
        ic = self if ic is None else ic
        if xs is None:
            xs = ic.xs
        return Interconnect(
            radius=ic.radius if radius is None else radius,
            width=ic.width if width is None else width,
            angle=ic.angle,
            xs=xs if xs is None else xs,
            layer=ic.layer if layer is None else layer,
            adapt_width=ic.adapt_width,
            adapt_xs=ic.adapt_xs,
            instantiate=ic.instantiate,
            pinstyle=ic.pinstyle if pinstyle is None else pinstyle,
            offset=ic.offset if offset is None else offset,
            varname=varname,  # Do not copy original 'varname'. This would overwrite the original ic.
            doc=f"copy: {ic.doc}" if doc is None else doc,
            PCB=ic.pcb if PCB is None else PCB,
            pitch=ic.pitch if pitch is None else pitch,
        )


    def _arc(
        self,
        radius=None,
        width=None,
        angle=None,
        xs=None,
        layer=None,
        offset=None,
        name=None,
    ):
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
        layer = self.layer if layer is None else layer
        xs = self.xs if xs is None else xs
        radius = self.radius if radius is None else radius
        width = self.width if width is None else width
        angle = self.angle if angle is None else angle
        offset = self.offset if offset is None else offset
        return self.Farc(
            radius=radius,
            width=width,
            angle=angle,
            xs=xs,
            layer=layer,
            offset=offset,
            name=name
        )


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
            pin = cp.here() if pin is None else pin
            width = pin.width if pin.width is not None else width
        if width is None:
            return self.width
            #try:
            #    width = nd.get_xsection(xs).width
            #except:
            #    width = self.width
        return width


    def _getwidths(self, pin=None, width=None, xs=None, widths=None, N=None, adapt=False):
        """Return width based on interconnect rules.

        Args:
            pin (Node): pin to connect
            width (float): waveguide width
            xs (str): xsection
            end (bool): True if output output side of taper

        Returns:
            width
        """
        if widths is not None:
            return widths
        if width is not None:
            return [width] * N
        if adapt or self.adapt_width:
            pin = cp.here() if pin is None else pin
            width = pin.width if pin.width is not None else width
        if width is None:
            return [self.width] * N
        return [width] * N


    def _getwidths2(self, N=None, pin=None, width=None, xs=None, config=None, adapt=False):
        """Return width based on interconnect rules.

        Args:
            pin (Node): pin to connect
            width (float): waveguide width
            xs (str): xsection
            end (bool): True if output output side of taper

        Returns:
            width
        """
        if config is not None:
            return config.width
        if isinstance(width, list):
            return width
        if width is not None:
            return [width] * N
        if adapt or self.adapt_width:
            pin = cp.here() if pin is None else pin
            width = pin.width if pin.width is not None else width
        if width is None:
            return [self.width] * N
        return [width] * N


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
            pin = cp.here() if pin is None else pin
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
            pin = cp.here() if pin is None else pin
            try:
                xs = pin.xs if pin.xs is not None else xs
            except:
                xs = self.xs
                # raise Exception('No xsection defined in pin. Add a xs attribute.')
        xs = self.xs if xs is None else xs
        return xs


    def _p2p_parse(
        self,
        pin1,
        pin2=None,
        xs=None,
        width1=None,
        width2=None,
        radius1=None,
        radius2=None
    ):
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


    def _p2p_parse2(
       self,
       N,
       pin1,
       pin2=None,
       xs=None,
       width1=None,
       width2=None,
       radius1=None,
       radius2=None,
       pitch=None,
       config=None,
    ):
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
       width1 = self._getwidths2(N, pin1, width1, xs, config)
       width2 = self._getwidths2(N, pin2, width2, xs, config)  # TODO: make specific for width2
       pitch = self._getpitch(N, pitch, config)
       radius1 = self._getradius(pin1, radius1, xs)
       radius2 = self._getradius(pin1, radius2, xs)
       ribwidth = 0
       ribradius = 0
       return pin1, pin2, xs, width1, width2, radius1, radius2, pitch, ribwidth, ribradius


    def _getribbon(self, N=None, pitch=None, width=None, radius=None):
        """
        """
        ribwidth = pitch[-1] + 0.5 * (width[0] + width[-1])
        if radius is None:
            ribradius = None
        else:
            ribradius = radius + 0.5 * pitch[-1]
        return ribwidth, ribradius


    def _getticks(self, N, pitch, ticks):
        """"""
        if ticks is not None:
            ticksarr = ticks
        elif pitch is not None:
            ticksarr = [pitch * i for i in range(N)]
        elif self.ticks is not None:
            ticksarr = self.ticks
        else:
            ticksarr = [self.pitch * i for i in range(N)]
        return ticksarr


    def _getpitch(self, N, pitch, config=None):
        """"""
        if config is not None:
            return config.pitch
        if isinstance(pitch, list):
            return pitch
        elif pitch is not None:
            pitch = [pitch * i for i in range(N)]
        elif self.ticks is not None:
            pitch = self.ticks
        else:
            pitch = [self.pitch * i for i in range(N)]
        return pitch


    def strt_solve(
        self,
        length=None,
        pin=None,
        xs=None,
        width=None,
        N=None,
        pitch=None,
        config=None,
    ):
        """Calculate strt solution.

        Args:
            length (float):
            pin (Node):
            xs (str):
            width (float):
            N (int):

        Returns:
            dict: solution parameters
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        N = max(1, self.N if N is None else N)
        _width = self._getwidths2(N, pin, width, xs, config)
        _pitch = self._getpitch(N, pitch, config)
        N = len(_pitch)

        if length is None:
            length = self.length

        geo_rib = {}
        start = {}
        for n in range(N):
            start[n] = (0, _pitch[n], 0)
            geo_rib[n] = [{
                'call': 'strt',
                'parameters': {
                    'xs': xs,
                    'length': length,
                    'width': _width[n],
                 },
            }]
        result = {
            'geo': geo_rib,
            'N': N,
            'solution': True,
            'message': '',
            'pin1': pin,
            'origin': 'strt',
            'start': start,
        }
        return result


    def strt(
        self,
        length=None,
        width=None,
        pin=None,
        xs=None,
        edge1=None,
        edge2=None,
        edgepoints=50,
        name=None,
        arrow=True,
        gridpatch=None,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create a straight waveguide.

        Args:
            length (float): length of guide in um
            width (float): width of guide in um
            pin (Node): optional Node for modeling info
            xs (str): optionals xsection of guide
            layer (int | str): layer number or layername
            edge1 (function): optional function F(t) describing edge1 of the waveguide, t=[0, 1].
            edge2 (function): optional function G(t) describing edge2 of the waveguide, t=[0, 1].
            edgepoints (int): optional number of points for edge1 and edge2 (default=50)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            gridpatch (float): patch gridsnap jumps at grid disconnect
                of cells with chamfers of size gridpatch. Default=0 is no patch.
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.strt_solve(
                length=length,
                width=width,
                pin=pin,
                xs=xs,
                N=N,
                pitch=pitch,
                config=config,
            )
        gridpatch = self.gridpatch if gridpatch is None else gridpatch
        result['kwargs'] = {
            'edge1': edge1,
            'edge2': edge2,
            'edgepoints': edgepoints,
            'gridpatch': gridpatch,
        }
        name = 'ic_strt' if name is None else name
        return self.ribbon(
            result=result,
            name=name,
            pinin=f"a{at1}",
            pinout=f"b{at2}",
            arrow=arrow,
            tubepins=tubepins,
        )


    def bend_solve(
        self,
        radius=None,
        angle=None,
        width=None,
        width2=None,
        pin=None,
        length1=0,
        length2=0,
        xs=None,
        offset=None,
        offset2=None,
        parabolic=None,
        N=None,
        pitch=None,
        config=None,
    ):
        """Calculate strt solution.

        Args:
            radius (float): bend radius
            angle (float): bend angle
            width (float): waveguide width
            width2 (float): Optional second width for tapered bends
            pin (Node): optional input pin to connect to
            length1 (float): optional length of straight before the bend (default=0)
            length2 (float: optional length of straight after the bend (default=0)
            xs (str): xsection name
            offset (float): optional bend offset at the beginning of the curve (default=None: xs will provide offset)
            offset2 (float): optional bend offset at the end of the curve (default=None: xs will provide offset)
            parabolic (bool): Makes the taper parabolic if true and linear otherwise. Default is True.
            N (int): optional number of waveguide in the interconnect ribbon

        Returns:
            dict: geometry solution
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        radius = self._getradius(pin, radius, xs)
        N = max(1, self.N if N is None else N)
        _width = self._getwidths2(N, pin, width, xs, config)
        if width2 is None:
            _width2 = _width
        else:
            _width2 = self._getwidths2(N, pin, width2, xs)
        _pitch = self._getpitch(N, pitch, config)
        N = len(_pitch)
        if offset is None:
            offset = self.offset
        if angle is None:
            angle = self.angle

        sections = 1
        if N > 1:
            sections = int((abs(angle) - 1e-10) // self.sectionangle + 1)
            if sections > 0:
                angle /= sections

        lengthfactor = m.tan(radians(0.5 * abs(angle)))
        geo_rib = {}
        start = {}
        for n in range(N):
            start[n] = (0, _pitch[n], 0)
            if angle < 0:
                dP = _pitch[n]
            else:
                dP = _pitch[N-1] - _pitch[n]
            elm1 = {
                'call': 'strt',
                'parameters': {
                    'xs': xs,
                    'length': length1 + dP * lengthfactor,
                    'width': _width[n],
                },
            }
            elm2 = {
                'call': 'bend',
                'parameters': {
                    'xs': xs,
                    'angle': angle,
                    'radius': radius,
                    'width': _width[n],
                    'width2': _width2[n],
                    'offset': offset,
                    'offset2': offset2,
                    'parabolic': parabolic,
                },
            }
            elm3 = {
                'call': 'strt',
                'parameters': {
                    'xs': xs,
                    'length': length2 + dP * lengthfactor,
                    'width': _width2[n],
                },
            }
            geo_rib[n] = sections * [elm1, elm2, elm3]

        result = {
            'geo': geo_rib,
            'N': N,
            'parameters': {  # to calculate xya later.
                'angle': [angle] * N,
                'radius': [radius] * N,
                'lengthfactor': lengthfactor,
           },
            'solution': True,
            'message': '',
            'origin': 'bend',
            'start': start,
        }
        return result


    def bend(
        self,
        radius=None,
        angle=None,
        width=None,
        width2=None,
        pin=None,
        xs=None,
        length1=0,
        length2=0,
        name=None,
        arrow=True,
        offset=None,
        offset2=None,
        parabolic=None,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create a bent waveguide (circular arc) with optinal in/out strt sections.

        Args:
            radius (float): radius at the center line of the arc in um
            width (float): width of the arc in um
            width2 (float): second width in case a of a taper.
            angle (float): angle of arc in degree (default=90)
            pin (Node): optional Node for modeling info
            xs (str): optiinal xsection of bend
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            offset (float): optional new offset for this bend only.
            length1 (float): length of an optional straight section before the bend
            length2 (float): length of an optional straight section after the bend
            parabolic (bool): Makes the taper parabolic if true and linear otherwise. Default is True.
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.bend_solve(
                radius=radius,
                angle=angle,
                width=width,
                width2=width2,
                pin=pin,
                length1=length1,
                length2=length2,
                xs=xs,
                offset=offset,
                offset2=offset2,
                parabolic=parabolic,
                N=N,
                pitch=pitch,
                config=config,
            )
        name = 'ic_bend' if name is None else name
        return self.ribbon(
            result=result,
            name=name,
            pinin=f"a{at1}",
            pinout=f"b{at2}",
            arrow=arrow,
            tubepins=tubepins,
        )


    def _bend_p2l_solve(self):
        """To be implemented."""
        pass


    # TODO: change to ribbon
    def bend_p2l(
        self,
        radius=None,
        angle=None,
        width=None,
        pin=None,
        xs=None,
        length1=0,
        ref=None,
        name=None,
        arrow=True,
        offset=None,
        max_length=None,
        tubepins=False,
    ):
        """Create a bent waveguide (circular arc) with optional strt input and ending at a line.

        Args:
            radius (float): radius at the center line of the arc in um
            width (float): width of the arc in um
            angle (float): angle of arc in degree (default=90)
            pin (Node): optional Node for modeling info
            xs (str): optional xsection of bend
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
        with nd.Cell(name=f'{name}_{xs}', instantiate=self.instantiate, cnt=True) as ICcell:
            refrel = nd.Pin().put(*T)
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self._line(length=length1, width=width, xs=xs).put(0)
            e2 = self._arc(radius=radius, width=width, angle=angle, xs=xs,
                offset=offset).put()
            self.strt_p2l(pin=e2.pin['b0'], ref=refrel, width=width, xs=xs,
                name=None, arrow=False, max_length=max_length).put()
            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            nd.Pin('a0', io=0, width=width, radius=radius, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, radius=radius, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
        if pin is not None:
            cfg.cp = pin
        return ICcell


    def taper_solve(
        self,
        length=None,
        pin=None,
        xs=None,
        width1=None,
        width2=None,
        shift=None,
        N=None,
        tapertype=None,
    ):
        """Calculate ptaper solution.

        Args:
           length (float):
           pin (Node):
           xs (str):
           width1 (float):
           width2 (float):
           shift (float): Output vertical position shift for skewed taper
           N (int):
           tapertype (partap): "lintap" or "partap" (default)

        Returns:
            dict: geometry solution
        """
        N = max(1, self.N if N is None else N)
        if tapertype is None:
            tapertype = 'taper'
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width1, xs)
        width2 = self._getwidth(pin, width2, xs)
        if length is None:
            length = self.length

        geo_rib = {}
        for n in range(N):
            geo_rib[n] = [{
                'call': tapertype,
                'parameters': {
                    'xs': xs,
                    'length': length,
                    'width1': width1,
                    'width2': width2,
                    'shift': shift,
                 }
            }]
            #geo_rib[n] = ((tapertype, (xs, length, width1, width2)))
        result = {
            'geo': geo_rib,
            'N': N,
            'solution': True,
            'message': '',
            'origin': 'taper',
        }
        return result


    def ptaper(
        self,
        length=None,
        width1=None,
        width2=None,
        pin=None,
        xs=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
    ):
        """Create a parabolic taper.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            pin (Node): optional Node for modeling info
            xs (str): optional xsection of taper
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.taper_solve(
                length=length,
                width1=width1,
                width2=width2,
                pin=pin,
                xs=xs,
                N=N,
                tapertype="ptaper",
            )
        name = 'ic_ptaper' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f"a{at1}", pinout=f"b{at2}", arrow=arrow, tubepins=tubepins)


    def taper(
        self,
        length=None,
        width1=None,
        width2=None,
        shift=0,
        xs=None,
        pin=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
    ):
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
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.taper_solve(
                length=length,
                width1=width1,
                width2=width2,
                pin=pin,
                xs=xs,
                N=N,
                shift=shift,
                tapertype="taper",
            )
        name = 'ic_taper' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f"a{at1}", pinout=f"b{at2}", arrow=arrow, tubepins=tubepins)


    def strt_p2l_solve(
            self,
            pin=None,
            ref=None,
            width=None,
            xs=None,
            N=None,
            max_length=None,
            pitch=None,
            config=None,
        ):
        """Calculate strt_p2l solution.

        Calculates length and uses "strt_solve" for ribbon dictionary generation.

        Args:
            pin (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            ref (Node | Instance)| tuple(x, y, a)): the reference line to intersect
            width (float): width of the interconnect in um
            xs (str): optional xsection of the strt
            N (int): optional number of waveguide in the ribbon.
            max_length (float): maximum length of the guide.

        Returns:
            dict: solution parameters"""
        pin = self._getpinout(pin)
        if ref is None:
            cfg.interconnect_errcnt += 1
            msg = "IC-{} strt_p2l needs a reference line via ref= keyword.".format(cfg.interconnect_errcnt)
            interconnect_logger(msg, 'error')
            ic_exception(msg)
            xs = 'error'
     #   width = self._getwidth2(N, pin, width, xs, config)
     #   pitch = self._getpitch(N, pitch, config)
     #   if xs != 'error':
     #       xs = self._getxs(pin, xs)
     #   N = max(1, self.N if N is None else N)
        if max_length is None:
            max_length = self.max_length

        pinb, T = nd.parse_pin(pin)
        pin = pinb.move(*T)
        refb, T = nd.parse_pin(ref)
        ref = refb.move(*T)

        x, y, a = nd.diff(pin, ref)
        L = x + y * m.tan(m.radians(a - 90))

        # rot = 0
        if abs(L) > max_length:
            L = sign(L) * max_length
            msg = "Solution for strt_p2l too large: >{}.".format(max_length)
            interconnect_logger(msg, 'warning')
            ic_exception(msg)
        if L < 0:
            # rot = 180
            pin = pin.rot(180)
            L = -L
            xs = 'error'
            cfg.interconnect_errcnt += 1
            msg = "IC-{} negative length for strt_p2l.".format(cfg.interconnect_errcnt)
            interconnect_logger(msg, 'error')
            ic_exception(msg)

        result = self.strt_solve(
            length=L,
            width=width,
            pin=pin,
            xs=xs,
            N=N,
            pitch=pitch,
            config=config,
        )
        return result


    #TODO: p2l_solve
    def strt_p2l(
        self,
        pin=None,
        ref=None,
        width=None,
        xs=None,
        edge1=None,
        edge2=None,
        edgepoints=50,
        name=None,
        arrow=True,
        max_length=None,
        gridpatch=None,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create a straight guide to intersect a reference line.

        p2l: point-to-line. Note there is no solution for a reference line
        parallel to the pointer in pin. To avoid huge (near parallel) lines,
        a max-length can be specified.

        Args:
            pin (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            ref (Node | Instance)| tuple(x, y, a)): the reference line to intersect
            width (float): width of the interconnect in um
            xs (str): optional xsection of the strt
            edge1 (function): optional function F(t) describing edge1 of the waveguide, t=[0, 1].
            edge2 (function): optional function G(t) describing edge2 of the waveguide, t=[0, 1].
            edgepoints (int): optional number of points for edge1 and edge2 (default=50)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            max_length (float): maximum length of the guide
            gridpatch (float): patch gridsnap jumps at grid disconnect
                of cells with chamfers of size gridpatch. Default=0 is no patch.
            N (int): optional number of waveguide in the ribbon.
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.strt_p2l_solve(
                pin=pin,
                ref=ref,
                width=width,
                xs=xs,
                N=N,
                max_length=max_length,
                pitch=pitch,
                config=config,
            )
        gridpatch = self.gridpatch if gridpatch is None else gridpatch
        result['kwargs'] = {
            'edge1': edge1,
            'edge2': edge2,
            'edgepoints': edgepoints,
            'gridpatch': gridpatch,
        }
        name = 'ic_strt_p2l' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f"a{at1}", pinout=f"b{at2}", arrow=arrow, tubepins=tubepins)


    def strt_p2p_solve(
        self,
        pin1=None,
        pin2=None,
        width=None,
        xs=None,
        N=None,
        at1=0,
        at2=0,
        pitch=None,
        config=None,
    ):
        """Find strt point-to-point solution.

        Args:


        Returns:
            dict: geometry solution
        """
        N = max(1, self.N if N is None else N)

        at1 = max(0, min(at1, N-1))
        at2 = max(0, min(at2, N-1))

        parse = self._p2p_parse(pin1=pin1, pin2=pin2, xs=xs, width1=width)
        pin1, pin2, xs, width, _, _, _ = parse

        x, y, a = nd.diff(pin1, pin2)
        length =  m.copysign(m.hypot(x, y), x)

        dis = self.pitch * (at2 - at1)
        anga = m.degrees(m.asin(dis / length))
        angc = m.degrees(m.atan2(y, x))
        angb = angc + anga

        geo_rib = {}
        for n in range(N):
            geo_rib[n] = [{
                'call': 'strt',
                'parameters': {
                    'xs': xs,
                    'length': abs(m.hypot(length, dis)),
                    'width': width,
                 }
            }]
        result = {
            'geo': geo_rib,
            'N': N,
            'solution': {'length': length},
            'found': True,
            'message': '',
            'N': N,
            'pin1': pin1.rot(angb),
            'origin': 'strt_p2p',
        }
        return result


    def strt_p2p(
        self,
        pin1=None,
        pin2=None,
        width=None,
        xs=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create point-to-point straight interconnect.

        Avoid usage of this component in optical interconnects.

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
        if result is None:
            result = self.strt_p2p_solve(
                pin1=pin1,
                pin2=pin2,
                width=width,
                xs=xs,
                N=N,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
            )
        name = 'ic_strt_p2p' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f"a{at1}", pinout=f"b{at2}", arrow=arrow, tubepins=tubepins)


    # TODO: add taper_p2p_solve and ribbons
    def taper_p2p(
        self,
        pin1=None,
        pin2=None,
        width1=None,
        width2=None,
        xs=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
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
        with nd.Cell(name=f'{name}_{xs}', instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self._taper(
                length=length,
                width1=width1,
                width2=width2,
                shift=shift,
                xs=xs
            ).put(0, 0, pin1.xya()[2])
            p1 = nd.Pin('a0', io=0, width=width1, radius=0, xs=xs).put(e1.pin['a0'])
            p2 = nd.Pin('b0', io=1, width=width2, radius=0, xs=xs).put(e1.pin['b0'])
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


    def _rot2ref_solve(
        self,
        pin=None,
        ref=None,
        angle=0,
        cw=None,
        length1=0,
        length2=0,
        width=None,
        radius=None,
        xs=None,
        N=None,
    ):
        """Calculate and return the angle to rotate from <pin> to reference direction <ref>.
        xs=None,
        Note that only the angle part of <ref> is used in the calcuation.

        Args:
            pin (Node): starting pin (default=cp)
            ref (Node): reference pin (default=org)
            angle (float): rotation with repect to ref in [Degrees] (default=0)
            cw (bool): angle direction clockwise or counter clockwise (default is shortest bend)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            N (int): optional number of waveguide in the ribbobnn

        Returns:
            dict: geometry solution
        """
        nd.cp.push()
        # TODO: allow to handle ref as tuple as a well
        if ref is None:
            ref = cfg.cells[-1].pin['org']
        else:
            ref = self._getpinout(ref)
        parse = self._p2p_parse(
            pin1=pin,
            xs=xs,
            width1=width,
            radius1=radius,
        )
        pin1, pin2, xs, width, _, radius1, radius2 = parse
        if pin2 is None:
            raise Exception('Source pin not specified in rot2ref.')

        # TODO: RB: does this ref need a .copy for leaving it out of the netlist?
        x, y, a = nd.diff(ref.rot(angle, drc=False), pin2)
        if a >= 180:
            a -= 360
        if cw is True:
            if a < 0:
                a += 360
        elif cw is False:
            if a > 0:
                a -= 360
        angle = -a
        nd.cp.pop()
        return self.bend_solve(
            radius=radius,
            angle=angle,
            width=width,
            pin=pin,
            length1=length1,
            length2=length2,
            xs=xs,
            N=N,
            # TODO: offset?
        )


    def rot2ref(
        self,
        pin=None,
        ref=None,
        angle=0,
        length1=0,
        length2=0,
        cw=None,
        width=None,
        xs=None,
        radius=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
    ):
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
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

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
        if result is None:
            result = self._rot2ref_solve(
                pin=pin,
                ref=ref,
                angle=angle,
                cw=cw,
                length1=length1,
                length2=length2,
                width=width,
                radius=radius,
                xs=xs,
                N=N,
                #at1=0,
                #at2=0,
            )
        if pin is not None:
            cfg.cp = pin
        if name is None:
            name = 'ic_rot2ref'
        return self.bend(result=result, name=name, at1=at1, at2=at2, arrow=arrow)


    def _construct_sbend(
        self,
        xs,
        angle,
        radius,
        width,
        offset,
        length1=0,
        length2=0,
        length3=0,
        N=None,
        removecommon=True,
        pitch=None,
        config=None,
    ):
        """Construct a sbend and filter out the common strt section between bends.

        Args:

        Returns:
            dict: sbend result
        """
        result1 = self.bend_solve(
            angle=-angle,
            radius=radius,
            width=width,
            length1=length1,
            xs=xs,
            N=N,
            pitch=pitch,
        )
        result2 = self.bend_solve(
            angle=angle,
            radius=radius,
            width=width,
            length2=length3,
            xs=xs,
            N=N,
            pitch=pitch,
        )
        geo = {}

        if removecommon:
            for n in range(N):
                result1['geo'][n] = result1['geo'][n][:-1]
                result2['geo'][n] = result2['geo'][n][1:]

        if length2 > 0:
            result_strt = self.strt_solve(xs=xs, length=length2, width=width, N=N, pitch=pitch)
            for n in range(N):
               geo[n] = result1['geo'][n] + result_strt['geo'][n] + result2['geo'][n]
        else:
            for n in range(N):
               geo[n] = result1['geo'][n] + result2['geo'][n]

        xya = (
            length1 + length3 + 2 * abs(radius * m.sin(radians(angle))) +\
                length2 * abs(radius * (1 - m.cos(radians(angle)))) +\
                0.5 * (pitch[-1] - pitch[0]) * (result1['parameters']['lengthfactor'] + result2['parameters']['lengthfactor']),
            offset,
            0
        )
        result = {
            'geo': geo,
            'N': N,
            'xya': xya,
            'origin': 'sbend',
            'start': result1['start'],

        }
        return result


    def sbend_solve(
        self,
        radius=None,
        width=None,
        pin=None,
        xs=None,
        offset=20,
        Ltot=None,
        length1=None,
        length2=None,
        Amax=90.0,
        N=None,
        at1=0,
        at2=0,
        pitch=None,
        config=None,
    ):
        """Calculate sbend solution.

        Args:
            See sbend()

        Returns:
            dict: geometry solution
        """
        N = max(1, self.N if N is None else N)
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
        radius  = self._getradius(pin, radius, xs)
        width  = self._getwidths2(N, pin, width, xs, config)
        pitch = self._getpitch(N, pitch, config)
        N = len(pitch)

        Lm, La, Lb = 0, 0, 0
        Amax = m.radians(Amax)

        if 2 * radius * (1 - m.cos(Amax)) > abs(offset):
            A = m.acos(1 - abs(offset) / (2 * radius))
            if offset > 0:
                A = -A
        else:
            Lm = (abs(offset) - abs(2 * radius * (1 - m.cos(Amax)))) / m.sin(Amax)
            A = sign(-offset) * Amax

        Lx = 2 * radius * m.sin(abs(A)) + Lm * m.cos(abs(A))
        dLx = abs(Ltot) - Lx

        if Ltot == 0:
            La = length1
            Lb = length2
        else:
            if length1 is not None:
                if La > dLx:
                    raise Exception("Ltot too short for length1={length1}")
                else:
                    La = length1
                    Lb = dLx - La
                    Ltot = abs(Ltot)
            elif length2 is not None:
                if Lb > dLx:
                    raise Exception("Ltot too short for length2={length2}")
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
        return self._construct_sbend(
            xs=xs,
            angle=m.degrees(A),
            radius=radius,
            width=width,
            offset=offset,
            length1=La,
            length2=Lm,
            length3=Lb,
            N=N,
            pitch=pitch,
            config=config,
        )


    def sbend(
        self,
        radius=None,
        width=None,
        pin=None,
        xs=None,
        offset=20,
        Ltot=None,
        length1=None,
        length2=None,
        name=None,
        arrow=True,
        Amax=90.0,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create an s-bend interconnect.

        Args:
            radius (float): bend radius at the center line of the arc in um
            width (float): width of the interconnect in um
            pin (Node): optional Node for modeling info
            xs (str): xsection of sbend
            offset (float): lateral offset of the sbend in um
            Ltot (float): optional total forward length of the sbend in um.
                When positive, additional length is added at the
                start of the s-bend, when negative it is added at the end,
                all provided the net forward length of the
                s-bend itself is shorter than abs(Ltot).
            length1 (float):
            length2 (float):
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            Amax: maximum angle of bends. default is 90
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.sbend_solve(
                xs=xs,
                width=width,
                radius=radius,
                pin=pin,
                offset=offset,
                Ltot=Ltot,
                length1=length1,
                length2=length2,
                Amax=Amax,
                N=N,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
            )
        name = 'ic_sbend' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f"a{at1}", pinout=f"b{at2}", arrow=arrow, tubepins=tubepins)


    #TODO: new tube option?
    def sinebend(
        self,
        width=None,
        pin=None,
        xs=None,
        distance=200,
        offset=20,
        name=None,
        arrow=True,
    ):
        """Create a sine-bend interconnect.

        This interconnect has zero curvature at both ends and can be
        connected without offsets to straight waveguides.
        See https://openepda.org/interconnect/interconnects.html

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
        with nd.Cell(
            name=f'{name}_{xs}_{int(distance)}_{int(offset)}',
            instantiate=self.instantiate,
            cnt=True,
        ) as ICcell:
            ICcell.group = 'interconnect'
            if abs(offset) < 1e-6:
                # Straight line will do just fine.
                e1 = self.strt(length=distance).put(0)
            else:
                e1 = self._sinecurve(width=width, distance=distance, offset=offset, xs=xs).put(0)
            nd.Pin('a0', io=0, width=width, radius=0, xs=xs).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width, radius=0, xs=xs).put(distance, offset, 0)
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)

        if pin is not None:
            cfg.cp = pin

        # Arc length of sinebend. Numeric integration is more accurate than
        # using the polyline points and doesn't take more time.
        def arclen_sb(t, d=distance, s=offset):  # function to integrate
            return s / (2 * m.pi) * sqrt((d / s) ** 2 + (1 - sin(t)) ** 2)

        ICcell.length_geo = quad(arclen_sb, 0, 2 * m.pi, args=(distance, offset))[0]

        return ICcell


    def sbend_p2p_solve(
        self,
        pin1=None,
        pin2=None,
        width=None,
        radius=None,
        xs=None,
        length1=0,
        doStrFirst=1,
        ortho=True,
        ref=None,
        Amax=90.0,
        N=None,
        at1=0,
        at2=0,
        pitch=None,
        config=None,
    ):
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
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            pitch(float | list):

        Returns:
            dict: geometry solution
        """
        N = max(1, self.N if N is None else N)
        parse = self._p2p_parse2(
            N,
            pin1,
            pin2,
            xs,
            width1=width,
            radius1=radius,
            pitch=pitch,
            config=config,
        )
        pin1, pin2, xs, width, _, radius1, radius2, pitch, _, _ = parse
        N = len(pitch)
        at1 = min(int(at1), N-1) if at1 >= 0 else -min(abs(int(at1)), N)
        at2 = min(int(at2), N-1) if at2 >= 0 else -min(abs(int(at2)), N)
        Wrib = pitch[-1] - pitch[0]
        pin1b = pin1.offset(0.5 * Wrib - pitch[at1])
        pin2b = pin2.offset(-0.5 * Wrib + pitch[at2])

        radius1N = radius1 if N == 1 else radius1 + 0.5 * Wrib
        #radius1N = radius1

        message = ''
        found = True
        A = 0
        Ltap = 0
        Amax = m.radians(Amax)

        # rotate to sbend-ref
        da1 = da2 = 0
        if ortho and ref is None:
            ref = pin1b
        if ref is not None:
            # make a pin out of ref:
            if isinstance(ref, (int, float)):
                ref = (0, 0, ref)
            ref, T = nd.parse_pin(ref)
            ref =  ref.move(*T)

            if ref is not pin1b:
                dx1, dy1, da1 = nd.diff(pin1b, ref)
                da1 = zeroDeg(da1)
                pin1b = pin1b.move(
                    abs(radius1N * m.sin(m.radians(da1))),
                    sign(da1) * radius1N * (1 - m.cos(m.radians(da1))),
                    da1
                )
            if ref is not pin2b:
                dx2, dy2, da2 = nd.diff(pin2b, ref.rot(180))
                da2 = zeroDeg(da2)
                pin2b = pin2b.move(
                    abs(radius1N * m.sin(m.radians(da2))),
                    sign(da2) * radius1N * (1 - m.cos(m.radians(da2))),
                    da2
                )

        xya = nd.diff(pin1b, pin2b)
        dx, dy, da = xya

        if length1 < 0:
            length1 = abs(length1)
            doStrFirst = 0

        # get dy from start and end position
        if length1 < Ltap:
            length1 = Ltap

        H = 2 * radius1N * (1 - m.cos(Amax)) + 2 * Ltap * m.sin(Amax)
        if H > abs(dy):
            # not enough offset for Amax--> no vertical straight guide section
            if Ltap == 0:
                A = m.acos(1 - abs(dy) / (2 * radius1N))
            else:
                tel = 0
                A = 0
                da = 0.2
                damin = 1e-8
                while abs(da) > damin and tel < 100:
                    if Ltap * m.sin(A) - radius1N * m.cos(A) < (abs(dy) - 2 * radius1N) / 2.0:
                        A += da
                    else:
                        A -= da
                        da /= 2.0
                        A += da
                    tel += 1
                if A < 10 * damin:
                    A=0

            A = sign(dy) * A
            Lstr = 0
        else: # use Amax angle
            A = sign(dy) * Amax
            Lstr = (abs(dy) - abs(2 * radius1N * (1 - m.cos(Amax))) - 2 * Ltap) / abs(m.sin(Amax))

        Lfit = dx - 2 * radius1N * abs(m.sin(A)) - (Lstr + 2 * Ltap) * abs(m.cos(A)) - length1
        La, Lb = 0, 0
        if doStrFirst == 1:
            La = length1 - Ltap
            Lb = Lfit - Ltap
        else:
            Lb = length1 - Ltap
            La = Lfit - Ltap

        if La < 0 or Lb < 0:
            found = False
            message = "(interconnect): No solution for sbend_p2p."

        sbend_result = self._construct_sbend(
            xs=xs,
            angle=m.degrees(-A),
            radius=radius1,
            width=width,
            offset=0,  # TODO: should not be 0
            length1=La,
            length2=Lstr,
            length3=Lb,
            N=N,
            removecommon=False,
            pitch=pitch,
            config=config,
        )
        sbend_result['pin1'] = pin1
        return sbend_result

        # else:
        #     if A == 0:
        #         Lstr += 4*Ltap
        # geo = geo1 + [
        #     ('strt', (xs, La, width)),
        #     ('bend', (xs, m.degrees(A), radius1, width)),
        #     ('strt', (xs, Lstr, width)),
        #     ('bend', (xs, -m.degrees(A), radius1, width)),
        #     ('strt', (xs, Lb, width))] + geo2
        # params = {
        #     'angle1': da1,
        #     'angle2': -da2,
        #     'angle': A,
        #     'length1': La,
        #     'length2': Lstr,
        #     'length3': Lb,
        #     'Lfit': Lfit,  # length left to reach end point.
        #     'message': message,
        #     'Amax': Amax,
        #     'Ltap': Ltap}
        # result = {
        #     'solution': found,
        #     'message': message,
        #     'parameters': params,
        #     'geo': geo,
        #     'xya': xya}
        # return parse, result


    def sbend_p2p(
        self,
        pin1=None,
        pin2=None,
        width=None,
        radius=None,
        Amax=90,
        xs=None,
        doStrFirst=1,
        Lstart=0,
        BendEndFlag=1,
        ref=None,
        name=None,
        arrow=True,
        bsb=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
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
            bsb (bool): if True, use bend_straight_bend_p2p() as fallback (default=True)


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
        if result is None:
            result = self.sbend_p2p_solve(
                xs=xs,
                width=width,
                radius=radius,
                pin1=pin1,
                pin2=pin2,
                length1=Lstart,
                doStrFirst=1,
                ortho=True,
                ref=None,
                Amax=Amax,
                N=N,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
            )

        # if not result['solution'] or abs(result['parameters']['angle']) < 0*m.pi:
        #     interconnect_logger(result['message'], 'info')
        #     #TODO: alterative back up if La is so long that Lb is negative?
        #     if bsb:
        #         msg = "replacing sbend with bend_strt_bend."
        #         interconnect_logger(msg, 'info')
        #         return self.bend_strt_bend_p2p(
        #             pin1=pin1,
        #             pin2=pin2,
        #             radius=radius,
        #             width=width,
        #             xs=xs,
        #             name=name,
        #             arrow=arrow,
        #         )
        #     else:
        #         return self.strt_p2p(pin1, pin2, xs='error')

        name = 'ic_sbend_p2p' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f'a{at1}', pinout=f"b{at2}", arrow=arrow, tubepins=tubepins)


    def bend_strt_bend_p2p_solve(
        self,
        pin1=None,
        pin2=None,
        radius=None,
        radius1=None,
        radius2=None,
        xs=None,
        width=None,
        ictype='shortest',
        length1=0,
        length2=0,
        N=None,
        at1=0,
        at2=0,
        pitch=None,
        config=None,
    ):
        """Calculate geometry for a bend_strt_bend interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius (float): optional bend radius for radius1 and radius2 in um
            radius1 (float): optional bend radius1 in um
            radius2 (float): optional bend radius2 in um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            ictype (str): interconnection type (default='shortest')
                options: 'shortest', 'll', 'lr', 'rl', rr', 'all'
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

        Returns:
            dict: geometry solution
        """
        N = max(1, self.N if N is None else N)
        if radius is not None:
            if radius1 is None:
                radius1 = radius
            if radius2 is None:
                radius2 = radius

        parse = self._p2p_parse2(
            N,
            pin1,
            pin2,
            xs=xs,
            width1=width,
            radius1=radius1,
            radius2=radius2,
            pitch=pitch,
            config=config,
        )
        pin1, pin2, xs, width, _, radius1, radius2, pitch, ribW, ribR = parse
        N = len(pitch)
        at1 = min(int(at1), N-1) if at1 >= 0 else max(abs(int(N + at1)), 0)
        at2 = min(int(at2), N-1) if at2 >= 0 else max(abs(int(N + at2)), 0)

        Wrib = pitch[-1] - pitch[0]
        radius1N = radius1 + 0.5 * Wrib
        radius2N = radius2 + 0.5 * Wrib

        Ltap = 0
        A = pin1.move(Ltap+length1, 0.5 * Wrib - pitch[at1], 0)  # to calculate the geometry with Ltap
        B = pin2.move(Ltap+length2, - 0.5 * Wrib + pitch[at2], 180)
        xya = nd.diff(A, B)
        dx, dy, da = xya

        radius1N -= 1e-8
        radius2N -= 1e-8
        # Calculate circle centers, where pin1 is put at (0, 0 ,0)
        c1Lx, c1Ly = 0, radius1N
        c1Rx, c1Ry = 0, -radius1N
        c2Lx, c2Ly = dx - radius2N * m.sin(m.radians(da)), dy + radius2N * m.cos(m.radians(da))
        c2Rx, c2Ry = dx + radius2N * m.sin(m.radians(da)), dy - radius2N * m.cos(m.radians(da))

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
                sx, sy = c2Rx - c1Rx, c2Ry - c1Ry
                d1, d2 = 1, -1
            if shape == 'rl':
                sx, sy = c2Lx - c1Rx, c2Ly - c1Ry
                d1, d2 = 1, 1
            if shape == 'lr':
                sx, sy = c2Rx - c1Rx, c2Ry - c1Ly
                d1, d2 = -1, -1
            if shape == 'll':
                sx, sy = c2Lx - c1Lx, c2Ly - c1Ly
                d1, d2 = -1, 1

            rr = d1 * radius1N + d2 * radius2N
            s = m.sqrt(sx**2 + sy**2)  # (sx,sy): (x,y)-distance between circle centers
            if rr > s + 1e-6:
                result = {
                    'found': False,
                    'message': "Radii too large (points to close) in bend_strt_bend in cell '{}' for shape '{}'. Maximum radius={:0.3f} or radius1+radius2<{:0.3f}. (radius1={:0.3f}, radius2={:0.3f})".\
                        format(cfg.cells[-1].cell_name, shape, 0.5*s, s, radius1, radius2)}
                #if len(shapes) == 1:
                interconnect_logger(result['message'], 'warning')
                found = False
                #return parse, result #found, 0, 0, 0, 0
            else:
                gs = m.atan2(sy, sx)  # angle through the circle centres at the start
                if abs(rr/s) <= 1:
                    found = True
                    gb = m.asin(rr / s)  # angle through the circle centers after placing connecting straight horizontal
                else:
                    found = False
                    result = {
                        'found': False,
                        'message': "No solution Found in bend_strt_bend in cell '{}'.".\
                            format(cfg.cells[-1].cell_name)}
                    interconnect_logger(result['message'], 'warning')
                    #return parse, result #found, 0, 0, 0, 0

                gt =  gb - gs  # angle of rotation of axis through circle centers from to put connection between circles horizontal.
                t1 = m.radians(0) + gt  # angle of spoke that points to start-point bsb on the circle
                t2 = m.radians(da) + gt  # angle of spoke that points to end-point bsb on the circle
                if d1 == 1:
                    b = negRad(-t1)  # angle of arc1
                else:
                    b = posRad(-t1)  # angle of arc1

                if d2 == -1:
                    e = negRad(t2)  # angle of arc2
                else:
                    e = posRad(t2)  # angle of arc2

                L = s * m.cos(gb)  # length of straight self.line
                Ltot = L + radius1N * abs(b) + radius2N * abs(e)  # total connection length
                # TODO: Ltot basad on radiusN is not the actual waveguide length.

                if b == 0:
                    L += Ltap
                else:
                    L -= Ltap
                if e == 0:
                    L += Ltap
                else:
                    L -= Ltap
                if L < 0:
                    found = True  # TODO: should be False?

            solutions[shape] = {
                'solution': {
                    'found': found,
                    'Ltot': Ltot,
                    'length1': length1,
                    'length2': L,
                    'length3': length2,
                    'angle1': b,
                    'angle2': e,
                }
            }

            result1 = self.bend_solve(
                angle=m.degrees(b),
                radius=radius1,
                width=width,
                length1=length1,
                xs=xs,
                N=N,
                pitch=pitch,
                config=config,
            )
            result2 = self.strt_solve(
                width=width,
                length=L,
                xs=xs,
                N=N,
                pitch=pitch,
                config=config,
            )
            result3 = self.bend_solve(
                angle=m.degrees(e),
                radius=radius2,
                width=width,
                length2=length2,
                xs=xs,
                N=N,
                pitch=pitch,
                config=config,
            )
            geo = {}
            for n in range(N):
                geo[n] = result1['geo'][n] + result2['geo'][n] + result3['geo'][n]

            solutions[shape]['geo'] = geo

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
                    'message': message,
                    'N': N,
                    'pin1': pin1,
                    'start': result1['start'],
                    'origin': 'bend_strt_bend',
                }
            else:
                raise Exception("No shortest solution found for bend_strt_bend.")
            return result

        elif ictype == 'all':
            #result = {'type': 'all'}
            if len(solutions) == 0:
                result['found': False]
            else:
                result = []
                for variation in solutions:
                    result.append({
                        'found': True,
                        'message': message,
                        'solution': solutions[variation]['solution'],
                        'geo': solutions[variation]['geo'],
                        'xya': xya,
                        'type': variation,
                        'N': N,
                        'pin1': pin1,
                        'origin': 'bend_strt_bend',
                   })
                #result['variations'] = curves
            return result

        elif ictype in ['rr', 'rl', 'lr', 'll']:
            result = {
                 'found': True,
                 'type': ictype,
                 'solution': solutions[ictype]['solution'],
                 'geo': solutions[ictype]['geo'],
                 'xya': xya,
                 'message': message,
                 'N': N,
                 'pin1': pin1,
            }
            return result

        else:
            result = {
                'found': False,
                'message': "No solution found in bend_strt_bend in cell '{}'.".\
                    format(cfg.cells[-1].cell_name)
            }
            return result


    def bend_strt_bend_p2p(
        self,
        pin1=None,
        pin2=None,
        radius=None,
        radius1=None,
        radius2=None,
        width=None, xs=None,
        length1=0,
        length2=0,
        ictype='shortest',
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
   ):
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
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            result (dict): pre-calculated result to draw interconnect ribbon
                (default=None: solve in place)

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
        if result is None:
            result = self.bend_strt_bend_p2p_solve(
                pin1,
                pin2,
                xs=xs,
                width=width,
                radius=radius,
                radius1=radius1,
                radius2=radius2,
                length1=length1,
                length2=length2,
                ictype=ictype,
                N=N,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
             )
        name = 'ic_bsb' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f'a{at1}', pinout=f'b{at2}', arrow=arrow, tubepins=tubepins)


# TODO: is this function needed?
    def bend_strt_bend(
        self,
        pin=None,
        radius=None,
        radius1=None,
        radius2=None,
        width=None,
        xs=None,
        ictype='shortest',
        name=None,
        arrow=True,
        pitch=None,
        config=None,
    ):
        """Generate a bend-straight-bend connection starting at the current pointer.

        This is the same connection as 'bend_strt_bend_p2p' with pin1 = cp.
        """
        return self.bend_strt_bend_p2p(
            pin1=cp.here(),
            pin2=pin,
            radius=radius,
            radius1=radius1,
            radius2=radius2,
            width=width,
            xs=xs,
            ictype=ictype,
            name=name,
            arrow=arrow,
            pitch=pitch,
            config=config,
        )


    def strt_bend_strt_p2p_solve(
        self,
        pin1=None,
        pin2=None,
        radius=None,
        xs=None,
        width=None,
        N=None,
        at1=0,
        at2=0,
        pitch=None,
        config=None,
    ):
        """Solve geometry for a strt_bend_strt interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius (float): optional first bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            pitch (float | list):

        Returns:
            dict: geometry solution
        """
        N = max(1, self.N if N is None else N)
        parse = self._p2p_parse2(
            N,
            pin1,
            pin2,
            xs=xs,
            width1=width,
            radius1=radius,
            pitch=pitch,
            config=config,
        )
        pin1, pin2, xs, width, _, radius1, _, pitch, _, _= parse
        N = len(pitch)
        at1 = min(int(at1), N-1) if at1 >= 0 else max(abs(int(N+at1)), 0)
        at2 = min(int(at2), N-1) if at2 >= 0 else max(abs(int(N+at2)), 0)
        Wrib = pitch[-1] - pitch[0]
        radius1N = radius1 + 0.5 * Wrib
        Ltap = 0
        pin1b = pin1.move(Ltap, 0.5 * Wrib - pitch[at1], 0)  # to calculate the geometry with Ltap
        pin2b = pin2.move(Ltap, -0.5 * Wrib + pitch[at2], 180)

        xya  = nd.diff(pin1b, pin2b)
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

            dx1 = abs(radius1N / m.tan(m.radians(g)))

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

        result1 = self.strt_solve(
            width=width,
            length=L1,
            xs=xs,
            N=N,
            config=config,
        )
        result2 = self.bend_solve(
            angle=da,
            radius=radius1,
            width=width,
            xs=xs,
            N=N,
            config=config,
        )
        result3 = self.strt_solve(
            width=width,
            length=L2,
            xs=xs,
            N=N,
            config=config,
        )
        geo = {}
        for n in range(N):
            geo[n] = result1['geo'][n] + result2['geo'][n] + result3['geo'][n]

        solution = {
            'length1': L1,
            'length2': L2,
            'angle': da,
            'message': msg
            }
        result = {
            'found': found,
            'solution': solution,
            'geo': geo,
            'xya': xya,
            'message': msg,
            'N': N,
            'pin1': pin1,
            'at1': at1,
            'at2': at2,
            'origin': 'strt_bend_strt',
            'start': result1['start'],
        }
        return result


    def strt_bend_strt_p2p(
        self,
        pin1=None,
        pin2=None,
        radius=None,
        width=None,
        xs=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create point-to-point straight-bend-straight interconnect.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius (float): optional bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            pitch (float | list):

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
        if result is None:
            result = self.strt_bend_strt_p2p_solve(
                pin1=pin1,
                pin2=pin2,
                xs=xs,
                width=width,
                radius=radius,
                N=N,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
            )
        if result['found']:
            name = 'ic_sbs' if name is None else name
            return self.ribbon(
                result,
                name,
                pinin=f"a{result['at1']}",
                pinout=f"b{result['at2']}",
                arrow=arrow,
                tubepins=tubepins,
            )

        else:
            interconnect_logger(result['message'], 'warning')
            ic_exception(result['message'])

        # goto bend_strt_bend:
        print("sbs switched to bsb")
        return self.bend_strt_bend_p2p(
            pin1,
            pin2,
            radius1=radius,
            radius2=radius,
            width=width,
            xs=xs,
            arrow=arrow,
            name=name,
            N=N,
            pitch=pitch,
            config=config,
        )


    def ubend_p2p_solve(
        self,
        pin1,
        pin2,
        length=0,
        xs=None,
        width=None,
        radius=None,
        balance=0,
        end_angle=False,
        N=None,
        at1=0,
        at2=None,
        pitch=None,
        config=None,
    ):
        """Calculate a ubend geometry between two pins.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            pin2 (Node | Instance | tuple(x, y, a)): end pin
            radius (float): optional bend radius for radius1 and radius2 in um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            balance (float): for a ubend <2*radius sidewyas, shift the horseshoe shape (default=0)
            end_angle (bool): Take pin2 angle into account when connecting if True (default=False)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

        Returns:
            dict: geometry solution
        """
        epsilon = 1e-6  # small distance addition to avoid floating-point-2R < 2R.

        N = max(1, self.N if N is None else N)
        parse = self._p2p_parse2(
            N,
            pin1,
            pin2,
            xs,
            width1=width,
            radius1=radius,
            pitch=pitch,
            config=config,
        )
        pin1, pin2, xs, width, _, radius1, radius2, pitch, _, _ = parse
        N = len(pitch)
        if at2 is None:
            at2 = at1
        at1 = min(int(at1), N-1) if at1 >= 0 else max(abs(int(N+at1)), 0)
        at2 = min(int(at2), N-1) if at2 >= 0 else max(abs(int(N+at2)), 0)
        Wrib = pitch[-1] - pitch[0]
        pin1b = pin1.offset(0.5 * Wrib - pitch[at1])
        pin2b = pin2.offset(-0.5 * Wrib + pitch[at2])

        radius1N = radius1 + 0.5 * Wrib
        #radius2N = radius2 + 0.5 * (N-1) * self.pitch

        xya = nd.diff(pin1b, pin2b)
        dx, dy, da = xya
        if end_angle:
            if da > 180:
                da -= 360
            ddx = radius1N * abs(m.sin(m.radians(da)))
            ddy = radius1N * m.copysign(1 - m.cos(m.radians(da)), -da)
            dx += ddx
            dy -= ddy
            p2 = pin2b.move(ddx, ddy, -da)

            result_end = self.bend_solve(
                angle=da,
                radius=radius1,
                width=width,
                N=N,
                pitch=pitch,
                config=config,
            )
        else:
            p2 = pin2b.rot(-da)
            result_end = {'geo': [[]] * N}
        da = 0

        if dx < 0:
            L2 = length - dx
            L1 = length
        else:
            L1 = length + dx
            L2 = length

        # determine if ubend should first widen by d1 and d2 on its respective pins
        if abs(dy) < 2 * radius1N:
            sign = -np.sign(dy)
            sign = sign if sign != 0 else 1
            d = sign * (epsilon + radius1N - 0.5 * abs(dy))
            d1 = d * (1 + balance) # right swing
            d2 = d * (1 - balance) # left swing
        else:
            d1, d2 = 0, 0

        geo = {}
        result1 = self.strt_solve(
            xs=xs,
            length=L1,
            width=width,
            N=N,
            pitch=pitch,
            config=config,
        )

        p1 = pin1b.move(L1, -0.5 * Wrib + pitch[0])
        p2 = p2.move(L2, 0.5 * Wrib - pitch[0]) # p2 is adjusted end pin for angle_end
        result2 = {'geo': [[]] * N}
        if abs(d1) > epsilon:
            d1 = np.sign(d1) * epsilon + d1
            result2 = self.sbend_solve(
                offset=d1,
                radius=radius1,
                width=width,
                xs=xs,
                N=N,
                pitch=pitch,
                config=config,
            )
            dx, dy, da = result2['xya']
            p1 = p1.move(dx, dy, da)

        result3 = {'geo': [[]] * N}
        if abs(d2) > epsilon:
            d2 = np.sign(d2) * epsilon + d2
            result3 = self.sbend_solve(
                offset=d2,
                radius=radius1,
                width=width,
                xs=xs,
                N=N,
                pitch=pitch,
                config=config,
            )
            # reverse result3
            for n, elms in result3['geo'].items():
                result3['geo'][n] = elms[::-1]
            dxs, dys, das = result3['xya']
            p2 = p2.move(dxs, -dys, -das)

        result4 = self.bend_strt_bend_p2p_solve(
            pin1=p1,
            pin2=p2,
            #radius=radius1,
            radius1=radius1,
            radius2=radius1,
            width=width,
            xs=xs,
            N=N,
            at1=0, #0.5 * (N-1),
            at2=0, #0.5 * (N-1),
            pitch=pitch,
            config=config,
        )
        result5 = self.strt_solve(
            xs=xs,
            length=L2,
            width=width,
            N=N,
            pitch=pitch,
            config=config,
        )

        for n in range(N):
            geo[n] = \
                result1['geo'][n] +\
                result2['geo'][n] +\
                result4['geo'][n] +\
                result3['geo'][n] +\
                result5['geo'][n] +\
                result_end['geo'][n]

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
            'xya': xya,
            'N': N,
            'pin1': pin1,
            'at1': at1,
            'at2': at2,
            'origin': 'ubend_p2p',
            'start': result1['start'],
        }
        return result


    def ubend_p2p(
        self,
        pin1=None,
        pin2=None,
        radius=None,
        width=None,
        xs=None,
        length=0,
        name=None,
        arrow=True,
        balance=0,
        end_angle=False,
        N=None,
        at1=0,
        at2=None,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        """Create point-to-point u-bend interconnect.

        An extra straight length can be added to the ubend with <length>.

        If the sideways translation needed in the ubend is < 2*radius, then
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
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

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
        if result is None:
            result = self.ubend_p2p_solve(
                pin1,
                pin2,
                length=length,
                width=width,
                radius=radius,
                balance=balance,
                end_angle=end_angle,
                N=N,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
            )
        name = 'ic_ubend_p2p' if name is None else name
        return self.ribbon(
            result=result,
            name=name,
            pinin=f"a{result['at1']}",
            pinout=f"b{result['at2']}",
            arrow=arrow,
            tubepins=tubepins,
        )


    def ubend(
        self,
        pin=None,
        offset=20.0,
        radius=None,
        width=None,
        xs=None,
        length=0,
        name=None,
        arrow=True,
        balance=0,
        end_angle=False,
        N=None,
        at1=0,
        at2=None,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
        #TODO: write a ubend_solve for this function
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
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

        Returns:
            Cell: ubend element
        """
        N = max(1, self.N if N is None else N)
        if pin is None:
            pin = cp.here()
        pitch = self._getpitch(N, pitch, config)
        N = len(pitch)
        Wrib = pitch[-1] - pitch[0]
        pin2 = pin.move(0, offset, 0) #+ Wrib - pitch[at1] - pitch[at2]
        return self.ubend_p2p(
            pin1=pin,
            pin2=pin2,
            radius=radius,
            width=width,
            xs=xs,
            length=length,
            name=name,
            arrow=arrow,
            balance=balance,
            end_angle=end_angle,
            N=N,
            at1=at1,
            at2=at2,
            result=result,
            pitch=pitch,
            config=config,
        )


    print_warning = True
    def pcurve_p2p(
        self,
        pin1=None,
        pin2=None,
        width=None,
        radius1=0,
        radius2=0,
        offset1=None,
        offset2=None,
        xs=None,
        name=None,
        arrow=True,
    ):
        if self.print_warning:
            print(
"""WARNING: the function 'pcurve_p2p' is now obsolete. You can obtain the
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
        return self.cobra_p2p(
            pin1=pin1,
            pin2=pin2,
            width1=width,
            width2=width,
            radius1=radius1,
            radius2=radius2,
            offset1=offset1,
            offset2=offset2,
            xs=xs,
            name=name,
            arrow=arrow,
        )


    def cobra_p2p_solve(
        self,
        pin1=None,
        pin2=None,
        width1=None,
        width2=None,
        radius1=0,
        radius2=0,
        parabolic=True,
        xs=None,
        at1=0,
        at2=0,
        N=None,
        pitch=None,
        config=None,
    ):
        """Find strt point-to-point solution.

        Returns:
            dict: geometry solution
        """

        N = max(1, self.N if N is None else N)
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

        parse = self._p2p_parse2(
            N,
            pin1,
            pin2,
            xs=xs,
            width1=width1,
            width2=width2,
            radius1=radius1,
            radius2=radius2,
            pitch=pitch,
            config=config,
        )
        pin1, pin2, xs, width1, width2, radius1, radius2, pitch, _, _ = parse
        N = len(pitch)
        at1 = min(int(at1), N-1) if at1 >= 0 else max(abs(int(N+at1)), 0)
        at2 = min(int(at2), N-1) if at2 >= 0 else max(abs(int(N+at2)), 0)
        Wrib = pitch[-1] - pitch[0]
        pin1b = pin1.offset(0.5 * Wrib - pitch[at1])
        pin2b = pin2.offset(-0.5 * Wrib + pitch[at2])

        dx, dy, da = nd.diff(pin1b, pin2b.rotate(180))
        xya = (dx, dy, da)

        geo_rib = {}
        start = {}
        for n in range(N):
            start[n] = (0, pitch[n] - pitch[at1], 0)
            geo_rib[n] = [{
                'call': 'cobra',
                'parameters': {
                    'xs': xs,
                    'xya': xya,
                    'radius1': radius1,
                    'radius2': radius2,
                    'width1': width1[n],
                    'width2': width2[n],
                    'shift': -0.5 * Wrib + pitch[n],
                    'parabolic': parabolic,
                 }
            }]
        result = {
            'geo': geo_rib,
            'N': N,
            'solution': True,
            'message': '',
            'N': N,
            'pin1': pin1,
            'at1': at1,
            'at2': at2,
            'origin': 'cobra',
            'start': start,
        }
        return result


    def cobra_p2p(
        self,
        pin1=None,
        pin2=None,
        width1=None,
        width2=None,
        radius1=0,
        radius2=0,
        offset1=None,
        offset2=None,
        parabolic=True,
        xs=None,
        name=None,
        arrow=True,
        N=None,
        at1=0,
        at2=0,
        result=None,
        tubepins=False,
        pitch=None,
        config=None,
    ):
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
        if result is None:
            result = self.cobra_p2p_solve(
                pin1=pin1,
                pin2=pin2,
                width1=width1,
                width2=width2,
                radius1=radius1,
                radius2=radius2,
                parabolic=parabolic,
                xs=xs,
                at1=at1,
                at2=at2,
                pitch=pitch,
                config=config,
                N=N,
            )
        name = 'ic_cobra_p2p' if name is None else name
        #name = f'{name}_{xs}_{int(dx)}_{int(dy)}_{int(da)}'
            # if abs(dy) < 1e-6 and abs(da) < 1e-6 and \
            #         abs(radius1) < 1e-6 and abs(radius2) < 1e-6:
            #     # Straight line (or taper) will do just fine.
            #     e1 = self.ptaper(length=dx, width1=width1, width2=width2).put(0)
            #     ICcell.Rmin = 0
            # else:
        return self.ribbon(
            result=result,
            name=name,
            pinin=f"a{result['at1']}",
            pinout=f"b{result['at2']}",
            arrow=arrow,
            tubepins=tubepins,
        )


    def euler_calibrate(self, radius=None, angle=None, scale=None):
        """Change calibration of the Euler bend in the Interconnect object.

        The scaling of the Euler will follow the (<radius>, <angle>) combo as
        provided to this function, or it will be set directly via the  <scale> value

        Args:
            radius (float) calibration radius. Default = self.radius
            angle (float): calibration angle. Default = 90 degrees
            scale (float): set scale directly, instead of via (radius, angle) combo

        Returns:
            None
        """
        # TODO: euler calibrate needs to go to xsection default is no scale info is provided.
        if scale is not None:
            if radius is not None or angle is not None:
                nd.main_logger(
                    f"Error: setting too many scale parameters at once. Using scale={scale}",
                    "error"
                )
        else:
            if radius is None:
                radius = self.radius
            if angle is None:
                angle = 90
        self._euler_base = nd.Tp_euler(
            width=None,
            width2=None,
            radius=radius,  # calibration
            angle=angle,  # calibration
            scale=scale,  # calibration instead of (radius, angle)
            xs=self.xs,
            layer=self.layer,
        )
        return None


    def add_scaled_euler(self, name="", scale=None, eulerangle=None):
        """Create the set of Euler functions for a specific scale.

        The scaled Euler functions are available in the dicts:
            scaled_euler
            scaled_euler_arc
            scaled_euler_arc_euler

        Args:
            name (str): key to select the specific <scale>.
            scale (float): scale factor of the euler underkey <name>
            eulerangle (float): maximum allowed euler angle for euler_arc and euler_arc_euler.

        Returns:
            None
        """
        self.scaled_euler[name] = partial(self.euler, scale=scale, radius_cal=None, angle_cal=None)
        self.scaled_euler_arc[name] = partial(self.euler_arc, scale=scale)
        self.scaled_euler_arc_euler[name] = partial(self.euler_arc_euler, scale=scale)


    def euler(
        self,
        pin: Node = None,
        width: float = None,
        width2: float = None,
        radius: float = None,
        angle: float = None,
        radius_cal: float = None,
        angle_cal: float = None,
        scale: float= None,
        xs: str = None,
        name: str = None,
        arrow: bool = True,
        parabolic: bool = None,
    ) -> Cell:
        """Create an Euler bend from a straight guide to a curvature of <radius>
        *or* to <angle>.

        NOTE: to maintain the Euler scale make sure to
        *only* provide the radius OR the angle to this method.
        If both are provided the Euler will have to scale to accommodate
        both parameters. To set the scale use radius_cal and angle_cal, optionally
        provide a "normal" angle or radius too.

        Args:
            pin1 (Node | Instance | tuple(x, y, a)): start pin (default=cp)
            width (float|function): optional waveguide width in um. This can
                be a function, w(t). In that case width2 is not used.
            width2 (float): optional waveguide width in um at end
            angle (float): end angle
            radius (float): end radius
            radius_cal (float): optional Euler radius for scale calibration for this bend only, use with angle_cal.
            angle_cal (float): optional Euler angle for scale calibration for this bend only, use with radius_cal.
            scale (float):
            xs (str): optional xsection
            radius (float): radius at start of the cobra (default=0 -> inf)
            name (str): optional new cell name for the component
            arrow (bool): draw connection arrows (default=True)
            parabolic (bool): uses a parabolic width profile if width1 != width2. Default is True.

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
        if radius is not None and angle is not None:
            nd.main_logger(
                f"Only radius *or* angle are allowed in the Euler to preserve scaling, but both are give. "
                f"If you want to overrule the scale use radius_cal and angle_cal, or the scale parameter explicitly."
                "error",
            )
        #if radius is None:
        #    radius = self._getradius(pin, radius, xs)
        elif angle is None and radius is None:
            angle = self.angle
        pinflip = not nd.get_xsection(xs).symmetry

        if name is None:
            name = 'ic_euler'
        with nd.Cell(name=f'{name}_{xs}', instantiate=self.instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            e1 = self._euler_base(
                width=width,
                width2=width2,
                radius=radius,
                angle=angle,
                radius_cal=radius_cal,
                scale=scale,
                angle_cal=angle_cal,
                xs=xs,
                parabolic=parabolic,
            ).put(0)
            ICcell.Rmin = radius  # TODO: adapt
            nd.Pin('a0', io=0, width=width, xs=xs, radius=0).put(e1.pin['a0'])
            nd.Pin('b0', io=1, width=width2, xs=xs, radius=radius).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=pinflip)
            trace.trace_stop()
            #ICcell.pin2 = pin2
        #cfg.cp = pin
        ICcell.length_geo = trace.trace_length()
        return ICcell


    def euler_arc_euler_solve(
        self,
        pin=None,
        width=None,
        width2=None,
        angle=90,
        eulerangle=None,
        arcangle=None,
        parabolic: bool = None,
        xs=None,
        name=None,
        arrow=True,
        radius_cal=None,
        angle_cal=None,
        scale=None,
        virtual_radius=False,
        get_virtual_radius=False,
        N=None,
        at1=0,
        at2=0,
        tubepins=False,
        _euler_arc=False,  # internal argument for solving for euler_arc
    ) -> Cell | float:
        """Create euler-arc-euler based bends.

        Possible bends results:
        - euler-euler
        - euler-arc-euler
        - strt-euler-arc-euler-strt

        Types of euler-arc-euler bends:
        1. "Automatic": The euler section will cover angles up to the maximum angle
            as allowed by the minimum radius. An arc bend covers the rest of angle.
        2. "Fixed <eulerangle>": Activates when explicitly setting the <eulerangle> to an
            angle smaller than the max angle. This type allows for large
            variation of angles all having equal Euler sections but different arc angles.
        3. "Virtual Radius": Similar to type "Fixed Euler", but each bend is
            extended with strt guides at the ends to create in/out positions
            corresponding to a virtual bend arc. Note that the minimum virtual radius
            is ultimatealy determined by the minimum angle required.

        Examples:
        1. Uniform Euler-uniform arc where arcangle = angle - eulerangle. Example::

            euler_arc_euler(width=1)

        2. Tapered Euler-uniform arc where arcangle = angle - eulerangle. Example::

            euler_arc_euler(width1=1, width2=2)


        Args:
            pin (Node): optional pin to connect to (default is current pin).
            width (float  | callable): alias for width1
            width1 (float  | callable): input width, can be a function w(t), t in [0, 1],
               then w(0) represents the input width.
            width2 (float): output width. Ignored if width1 is a callable.
            angle (float): angle of total bend (default=90).
            eulerangle (float): If set, each euler section will be clipped
                to this angle. (default=None).
            arcangle (float): If None, only the Eulers will be tapered.
            parabolic (bool): uses a parabolic width profile if width1 != width2. Default is True.
            virtual_radius (bool | float): activate virual radius geometry if True
                (default=False), or to set an explicit virtual radius.
            get_virtual_radius (bool): returns the virtual arc radius (no cell)
                to with all euler-arc-eulers will fit.
            radius_cal (float): optional calibration radius overrule for euler
                rescaling for this element only. Needs an angle_cal arg too.
            angle_cal (float): optional calibration angle overrule for euler
                rescaling for this element only. Needs a radius_cal arg too.
            N (int): optional number of waveguide in the ribbon.
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

        Returns:
            Cell | float: euler-arc-euler element dict | virtual radius value
        """
        nd.cfg.check_basename = False
        # avoid that eulers from two different interconnect (having other function id) raises
        # a "WARNING: Reusing a basename across function IDs"

        if arcangle is not None and arcangle < 0.0:
            raise ValueError(f"arcangle must be a float > 0 but {arcangle} was provided.")

        epsilon = 1e-6
        if _euler_arc:
            eulercnt = 1.0
        else:
            eulercnt = 2.0

        N = max(1, self.N if N is None else N)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width, xs)
        angle = self.angle if angle is None else angle
        width2 = width1 if width2 is None else width2

        # Check widths
        #if width1 is None:
        #    width1 = width
        #if width is not None and width1 is not None:
        #    nd.main_logger(
        #        f"Only provide 'width or width1, not both. Continuing with width1 = {width1:0.3f}",
        #        "warning"
        #    )


        # =============================================================================
        #   Calculate Rmax and euler-angle
        # =============================================================================
        signangle = sign(angle)
        angle = abs(angle)

        _, _, maxangle = self._euler_base(
            radius_cal=radius_cal,
            angle_cal=angle_cal,
            scale=scale,
            xs=xs,
            solve=True,
        )
        if eulerangle is None and arcangle is None:
            eulerangle = maxangle
            strictEangleFlag = False
        elif eulerangle is None and arcangle is not None:
            eulerangle = (angle - arcangle) / eulercnt
            strictEangleFlag = True
        else:
            strictEangleFlag = True

        if angle < eulercnt * eulerangle:
            if strictEangleFlag:
                nd.main_logger(
                    f"Requested angle below minimum euler angle constraint: {angle:.3f} < {2*eulerangle:.3f}. "\
                    "Will reduce the Euler angle to proceed with output.",
                    "error"
                )
            eulerangle = angle / eulercnt
            alpha = 0  # Total arc bend angle

        if eulerangle > maxangle + epsilon:
            nd.main_logger(
                f"Euler angle too large: eulerangle > maxangle: {eulerangle:.3f} > {maxangle:.3f}. "\
                "Setting eulerangle = maxangle.",
                "error"
            )
            eulerangle = maxangle
            maxangleFlag = True
        else:
            maxangleFlag = False

        Eradius, _, _ = self._euler_base(
            radius_cal=radius_cal,
            angle_cal=angle_cal,
            scale=scale,
            angle=eulerangle,
            xs=xs,
            solve=True,
        )
        Rarcmax = 2 * Eradius  # the arc radius that is always equal or larger than any double euler bend.
        if get_virtual_radius:
            return Rarcmax

        if not isinstance(virtual_radius, bool):
            if Rarcmax > virtual_radius:
                nd.main_logger(
                    f"Setting virtual euler radius below save value: {virtual_radius:.3f} < {Rarcmax:.3f}.",
                    "warning"
                )
            else:
                Rarcmax = virtual_radius

        alpha = angle - eulercnt * eulerangle  # Total arc angle (tapered + uniform)

        if arcangle is None:
            tapered_angle = 0.0
            arcangle = alpha
            width_mid = width2
        else:
            if arcangle > angle - eulerangle * eulercnt:
                nd.main_logger(
                    f"Arc angle too large: arcangle > angle - eulerangle * "\
                    f"{int(eulercnt)}: {arcangle:.3f} > {angle - eulerangle * eulercnt:.3f}. "\
                    f"Setting arcangle = angle - eulerangle * {int(eulercnt)} = {angle - eulerangle * eulercnt:.3f}.",
                    "error"
                )
                arcangle = angle - eulerangle * eulercnt
                width_mid = width2
                tapered_angle = 0.0
            else:
                tapered_angle = (alpha - arcangle) / eulercnt

                if width2 == width1:
                    width_mid = width1
                else:
                    if parabolic:
                        width_mid = parabolicvarwidth(width1, width2, eulerangle / (eulerangle + tapered_angle))
                    else:
                        width_mid = linvarwidth(width1, width2, eulerangle / (eulerangle + tapered_angle))

        euler = self._euler_base(  # always use positive angle here to make eulers same (hashed) cell.
            radius_cal=radius_cal,
            angle_cal=angle_cal,
            scale=scale,
            angle=eulerangle,
            xs=xs,
            width=width1,
            width2=width_mid,
        )
        eulerflip = True if signangle < 0 else False
        x, y, a = euler.pin['b0'].xya()
        y = abs(y)

        sections = 1  # TODO: do not allow mutiple section for euler_arc = True?
        if N > 1:
            # use for drawing in rib only mode:
            width_rib = (self.pitch - 1) * N + width1
            sections = int((abs(angle) - 1e-10) // self.sectionangle + 1)
            if sections > 0:
                angle /= sections

        Lstrt = 0

        #if strictEangleFlag:
        #    # Calculate strt for strt-euler-arc-euler-strt such that it exactly fits in the Rmax arc.
        if virtual_radius:
            angle1 = radians(0.5 * alpha + eulerangle)
            hm = Rarcmax * sin(radians(0.5 * angle))
            ha = Eradius * sin(radians(0.5 * alpha))
            t1 = y * sin(angle1)
            t2 = x * cos(angle1)
            tcomp = hm - (ha + t1 + t2)
            Lstrt = tcomp / cos(angle1)

        nd.cfg.check_basename = True  # switch back on again

        # draw strt-euler-bend-euler-strt that matches a Rmax Arc bend:
        # with nd.Cell(f'euler2_{angle:.3f}', cnt=True) as C:
        geo_rib = {}
        if not _euler_arc:
            for n in range(N):
                if signangle < 0:
                    j = n
                else:
                    j = N - n - 1
                length = j * self.pitch * m.tan(radians(0.5 * angle))
                geo_rib[n] = [
                    {'call': 'strt',
                     'parameters': {
                         'length': Lstrt + length,
                         'width': width1,
                         'xs': xs,
                     }
                    },
                    # TODO: add serializer option for euler
                    {'call': None, 'cell': euler, 'flip': eulerflip, 'parabolic': parabolic,},
                    {'call': 'arc',    # arc is static
                     'parameters': {
                         'xs':xs,
                         'angle': tapered_angle * signangle,
                         'offset': 0,
                         'radius': Eradius,
                         'width': width_mid,
                         'width2': width2,
                         'parabolic': parabolic,
                     }
                    },
                    {'call': 'arc',  # arc is static
                     'parameters': {
                         'xs':xs,
                         'angle': arcangle * signangle,
                         'offset': 0,
                         'radius': Eradius,
                         'width': width2,
                     }
                    },
                    {'call': 'arc',   # arc is static
                     'parameters': {
                         'xs':xs,
                         'angle': tapered_angle * signangle,
                         'offset': 0,
                         'radius': Eradius,
                         'width': width2,
                         'width2': width_mid,
                         'parabolic': parabolic,
                     }
                    },
                    {'call': None, 'cell': euler, 'reverse': True, 'flip': eulerflip,'parabolic': parabolic,},
                    {'call': 'strt',
                     'parameters': {
                         'length': Lstrt + length,
                         'width': width1,
                         'xs': xs,
                     }
                    },
                ]

            result = {
                'geo':  geo_rib,
                'N': N,
                'origin': 'euler_arc_euler',
            }
            return result

        else:  # _euler_arc is True
            geo_rib[0] = [
                {'call': None, 'cell': euler, 'flip': eulerflip},
                {'call': 'arc',  # arc is static
                    'parameters': {
                        'xs':xs,
                        'angle': tapered_angle * signangle,
                        'offset': 0,
                        'radius': Eradius,
                        'width': width_mid,
                        'width2': width2,
                        'parabolic': parabolic,
                    }
                },
                {'call': 'arc',  # arc is static
                 'parameters': {
                     'xs':xs,
                     'angle': arcangle * signangle,
                     'offset': 0,
                     'radius': Eradius,
                     'width': width2,
                 }
                },
            ]

            result = {
                'geo':  geo_rib,
                'N': 1,
                'origin': 'euler_arc',
            }
            return result


    def euler_arc_solve(
        self,
        pin: Node = None,
        width: float = None,
        width2: float = None,
        angle: float = 90,
        eulerangle: float = None,
        arcangle: float = None,
        parabolic: bool = None,
        xs: str = None,
        name: str = None,
        arrow: bool = True,
        radius_cal: float = None,
        angle_cal: float = None,
        scale: float= None,
        virtual_radius: bool = False,
        get_virtual_radius: bool = False,
        N: int = None,
        at1: int = 0,
        at2: int = 0,
    ) -> dict:
        """Internal routine to draw Euler-arc bends.

        Args:
            see euler_arc_euler_solve()

        Returns:
            dict: Ribbon dictionary of the Euler-arc bend
        """
        return self.euler_arc_euler_solve(
            pin=pin,
            width=width,
            width2=width2,
            angle=angle,
            eulerangle=eulerangle,
            arcangle=arcangle,
            parabolic=parabolic,
            xs=xs,
            name=name,
            arrow=arrow,
            radius_cal=radius_cal,
            angle_cal=angle_cal,
            scale=scale,
            virtual_radius=virtual_radius,
            get_virtual_radius=get_virtual_radius,
            N=N,
            at1=at1,
            at2=at2,
            _euler_arc=True,
        )


    def euler_arc(
        self,
        pin: Node = None,
        width: float = None,
        width2: float = None,
        angle: float = 90,
        eulerangle: float = None,
        arcangle: float = None,
        parabolic: bool = True,
        xs: str = None,
        name: str = None,
        arrow: bool = True,
        radius_cal: float = None,
        angle_cal: float = None,
        scale: float= None,
        N: int = None,  # remove?
        at1: int = 0,  # remove?
        at2: int = 0,  # remove?
        result: dict = None,
        tubepins: bool = False,
    ) -> Cell:
        """Draw an euler-arc bend.

        Three tapering options are available:

            1. Uniform Euler-uniform arc: euler_arc(width=1), where arcangle = angle - eulerangle;

            2. Tapered Euler-uniform arc: euler_arc(width1=1, width2=2), where arcangle = angle - eulerangle;

            3. Tapered Euler-tapered arc-uniform arc: euler_arc(width1=1, width2=2, arcangle=30),
              where the tapered arc extends for angle - eulerangle - arcangle.

        Args:
            pin (Node): optional pin to connect to (default is current pin).
            width (float): alias for width1
            width1 (float): input width
            width2 (float): output width
            angle (float): angle of total bend (default=90).
            eulerangle (float): If set, each euler section will be clipped to this angle. Default=None.
            arcangle (float): Angle of the arc section with uniform width. If None, only the Eulers will be tapered.
            If None, arcangle = angle - eulerangle, so only the Euler sections will be tapered.
            xs (str): optional xsection name
            name (str): Name of the cell
            arrow (bool): draw connection arrows (default=True)
            virtual_radius (bool | float): activate virual radius geometry if True
                (default=False), or to set an explicit virtual radius.
            get_virtual_radius (bool): returns the virtual arc radius (no cell)
                to with all euler-arc-eulers will fit.
            radius_cal (float): optional calibration radius overrule for euler
                rescaling for this element only. Needs an angle_cal arg too.
            angle_cal (float): optional calibration angle overrule for euler
                rescaling for this element only. Needs a radius_cal arg too.
            N (int): optional number of waveguide in the ribbon.
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.

        Returns:
            Cell: euler-arc element
        """
        if result is None:
            result = self.euler_arc_solve(
                pin=pin,
                width=width,
                width2=width2,
                angle=angle,
                eulerangle=eulerangle,
                arcangle=arcangle,
                parabolic=parabolic,
                xs=xs,
                arrow=arrow,
                radius_cal=radius_cal,
                angle_cal=angle_cal,
                scale=scale,
                virtual_radius=False,
                get_virtual_radius=False,
                N=1,
                at1=at1,
                at2=at2,
            )
        name = 'ic_euler_arc' if name is None else name
        return self.ribbon(
            result=result,
            name=name,
            pinin=f'a{at1}',
            pinout=f'b{at2}',
            arrow=arrow,
            tubepins=tubepins,
        )


    def euler_arc_euler(
        self,
        pin: Node = None,
        width: float = None,
        width2: float = None,
        angle: float = 90,
        eulerangle: float = None,
        arcangle: float = None,
        parabolic: bool = True,
        xs: str = None,
        name: str = None,
        arrow: bool = True,
        radius_cal: float= None,
        angle_cal: float= None,
        scale: float= None,
        virtual_radius: bool = False,
        get_virtual_radius: bool = False,
        N: int =None,
        at1: int = 0,
        at2: int = 0,
        result: dict = None,
        tubepins: bool = False,
    ) -> Cell:
        """Draw a euler-arc-euler interconnect

        Args:
           See euler_arc_euler_solve().

        Returns:
            Cell: strt_bend_strt element

        Example:
            Create and place a straight-bend-straight guide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.euler_arc_euler()
                guide.put()
                nd.export_plt()
        """
        if result is None:
            result = self.euler_arc_euler_solve(
                pin=pin,
                width=width,
                width2=width2,
                angle=angle,
                eulerangle=eulerangle,
                arcangle=arcangle,
                parabolic=parabolic,
                xs=xs,
                arrow=arrow,
                radius_cal=radius_cal,
                angle_cal=angle_cal,
                scale=scale,
                virtual_radius=virtual_radius,
                get_virtual_radius=get_virtual_radius,
                N=N,
                at1=at1,
                at2=at2,
            )
        if get_virtual_radius:
            return result
        name = 'ic_euler_arc_euler' if name is None else name
        return self.ribbon(result=result, name=name, pinin=f'a{at1}', pinout=f'b{at2}', arrow=arrow, tubepins=tubepins)


    # TODO: use ribbons:
    def mamba(
        self,
        points,
        radius=None,
        width=None,
        pin=None,
        xs=None,
        N=None,
        pitch=10,
        offset=0,
        polyline=True,
        showpins=False,
        name=None,
        arrow=True,
        at1=0,
        at2=0,
        tubepins=False,
    ):
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
            polyline (bool): boolean determining if the mamba is also drawn as polyline
                (default=True)
            showpins (bool): show the points as dots in the layout (default=False)
            name (str): optional new name for the component
            arrow (bool): draw connection arrows (default=True)
            N (int): optional number of waveguide in the ribbobnn
            at1 (int): optional ribbon pin "at" which start pin1 connects to. default=0.
            at2 (int): optional ribbon pin "at" which end pin2 connects to. default=0.
            offset (float): offset=0 (default) center the (ribbon) interconnect on the mamba points.


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
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        radius = self._getradius(pin, radius, xs)
        #pinflip = not nd.get_xsection(xs).symmetry

        ring = nd.Polygon(points=geom.circle(radius=0.5 * width), layer=self.layer)
        if N is None:
            N = self.N
        offset = int(0.5 * (N - 1) + offset)

        #create points along the mamba for interconnects:
        p1, p2 = [], []
        size = len(points)
        for i in range(size - 1):
            dx = points[i+1][0] - points[i][0]
            dy = points[i+1][1] - points[i][1]
            a = np.degrees(m.atan2(dy, dx))
            p1.append((points[i][0], points[i][1], a))
            p2.append((points[i+1][0], points[i+1][1], a + 180))

        #print('p1[0]:', p1[0])
        if name is None:
            name = 'ic_mamba'
        with nd.Cell(name=f'{name}_{xs}', instantiate=False, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            ICcell.default_pins('pla0', 'plb0')
            start = nd.Pin(width=width, xs=xs).put(p1[0])
            nd.Pin('pla0', width=width, xs=xs).put(start.rot(180))

            pin1 = start
            for i in range(size - 2): # loop over points
                pin2 = nd.Pin(width=width, xs=xs).put(p2[i+1])
                result = self.strt_bend_strt_p2p_solve(
                    pin1=pin1,
                    pin2=pin2,
                    radius=radius,
                    N=N,
                    at1=offset,
                    at2=offset,
                    pitch=pitch,
                )
                solution = result['solution']
                length1 = solution['length1']
                length2 = solution['length2']
                angle = solution['angle']
                #message  = solution['message']

                if result['found']:
                    self.strt(
                        length=length1,
                        pin=pin1,
                        N=N,
                        at1=offset,
                        at2=offset,
                        arrow=False,
                        pitch=pitch
                    ).put()
                    self.bend(
                        angle=angle,
                        radius=radius,
                        N=N,
                        at1=offset,
                        at2=offset,
                        arrow=False,
                        pitch=pitch,
                    ).put()
                    if i is size - 3:
                        self.strt(
                            length=length2,
                            N=N,
                            at1=offset,
                            at2=offset,
                            arrow=False,
                            pitch=pitch,
                        ).put()

                else:
                    if i < size - 3:
                        self.bend(
                            radius=radius,
                            pin=result['pin1'],
                            angle=angle,
                            N=N,
                            at1=offset,
                            at2=offset,
                            arrow=False,
                            pitch=pitch,
                        ).put()
                        last = cp.here()
                    else:
                        self.bend_strt_bend_p2p(
                            pin1=result['pin1'],
                            pin2=pin2,
                            radius=radius,
                            N=N,
                            at1=offset,
                            at2=offset,
                            arrow=False,
                            pitch=pitch,
                        ).put()

                pin1 = cp.here()
                if showpins:
                    plast = nd.Pin().put(p2[-1])
                    ring.put(result['pin1'])

            #nd.Pin(f'a{num}', type=1).put(cp.here())
            #end = nd.Pin(width=width, xs=xs).put(p2[-1])
            #nd.Pin('b0').put(end.rot(180))
            #if arrow:
            #    self.arrow.put(ICcell.pin['pla0'])
            #    self.arrow.put(ICcell.pin['plb0'], flip=pinflip)

            if polyline is True:
                if pin is None:
                    nd.Polyline(points=points, width=2, layer=self.mambalayer).put(0)
                else:
                    nd.Polyline(points=points, width=2, layer=self.mambalayer).put(p1[0][0], p1[0][1])
        return ICcell


    @nd.bb_util.trial_cell
    def ribbon_cells_only(self, result, cell=False):
        """Create list of cells for each ribbon tube in <result> and no mask output.

        This method cleans up any cells it creates and avoids putting them
        in the first place where possible. Also it supresses any shape generation,
        hence no polygons, etc., are calculated in mask_elements calls.

        This method can be utilized to explore cell properties along the tubes
        with less overhead than path tracing down stream.

        Args:
            result (bool): tube structure.
            cell (bool): return the actual cell (always False, needs to be there for function profile)

        Returns:
            dict: {tube number: list of Cells in tube}
        """
        shapes_store = cfg.generate_shapes
        cfg.generate_shapes = False
        kwargs = result.get('kwargs', {})
        tubes = {}
        for n in range(result['N']):
            tubes[n] = self.tube(geo=result['geo'][n], kwargs=kwargs, put=False)
        cfg.generate_shapes = shapes_store
        return(tubes)


    def ribbon(
        self,
        result,
        name=None,
        pinin='a0',
        pinout='b0',
        arrow=True,
        tubepins=False,
    ):
        """Create ribbon cell or single tube cell from a <result> dict.

        Args:
            result (dict): result of an interconnect solver function.
            name (str): interconnect name
            pinin (str): default ribbon input pin (default='a0').
            pinout (str): default ribbon output pin (default='b0').
            arrow (bool): create pin arrow at ribbon io (default=True).
            raise_pins (bool): Raise tube reference pins of outer tubes if True
                (default=False).

        Example::

            result = {
                'geo': {
                    0: [{<element>}, ...]  # tube 0
                    1: [{<element>}, ...]  # tube 1
                }
                'N': 2  # number of tubes in "result"
                'solution': True,  # optional
                'message': '',  # optional
                'pin1': "pin"?
                'origin': <origin function nam that created "result">
                'start': [  # optional, start position of the tubes.
                    (x0, y0, a0),
                    (x1, y1, a1),
                ]
            }

        Returns:
            Cell: ribbon of interconnects or single interconnect.
        """
        RIB = False
        if not isinstance(result, list):
            results = [result]
        else:
            results = result
        cnt = True
        if name is None:
        #    cnt = True
            name = "noname"
        ICcells = []
        for result in results:
            if result['N'] > 1 or cfg.always_use_ribbon:
                RIB = nd.Cell(f"ribbon_{name}_N{result['N']}", instantiate=cfg.instantiate_ribbon, cnt=cnt)
                RIB.default_pins(pinin, pinout)

            kwargs = result.get('kwargs', {})
            start = result.get('start', None)
            for n in range(result['N']):
                _tubepins = False
                if tubepins is True and n in [0, result['N'] - 1]:  # only tubepins on outer guides
                    _tubepins = True
                ICcell = self.tube(
                    name=name,
                    geo=result['geo'][n],
                    kwargs=kwargs,
                    arrow=arrow,
                    tubepins=_tubepins,
                )
                if RIB:
                    if start is None:
                        ribinfo = result.get('ribbon', False)
                        if ribinfo:
                            e2 = ICcell.put(0, ribinfo[0][n])
                        else:
                            e2 = ICcell.put(0, self.pitch * n)
                    else:
                        e2 = ICcell.put(start[n])
                    nd.Pin(f'a{n}', pin=e2.pin['a0']).put()
                    nd.Pin(f'b{n}', pin=e2.pin['b0']).put()
                    nd.put_stub([f'a{n}', f'b{n}'], pinstyle=self.pinstyle, length=0)

                    if _tubepins:
                        namei, nameo = [], []
                        side = 'r' if n == 0 else 'l'
                        for pin in e2.pin:
                            if pin.startswith(cfg.tubepin_prefix):
                                namei.append(pin)
                                nameo.append(f"{cfg.tubepin_prefix}{side}{pin[len(cfg.tubepin_prefix):]}")
                        e2.raise_pins(namei, nameo)
                        RIB.put_stub(nameo, pinstyle='tube')
            pin1 = result.get('pin1', None)
            if RIB:
                RIB.put_stub(["a0", "b0"], shape=self.ribbon0)
                RIB.put_stub([f"a{result['N']-1}", f"b{result['N']-1}"], shape=self.ribbonN)
                RIB.close()
                if pin1 is not None:
                    cfg.cp = pin1
                return RIB
            if pin1 is not None:
                cfg.cp = pin1
            ICcells.append(ICcell)

        if len(ICcells) == 1:
           return ICcells[0]
        else:
           return ICcells


    def tube(
        self,
        geo,
        showpins=False,
        name=None,
        xs=None,
        arrow=True,
        kwargs=None,
        put=True,
        tubepins=False,
        n=None,
    ):
        """Draw interconnect based on tubes.

        Each tube symbol is a dictionary that describes how to create a cell to
        place in the tube, or it points to an existing cell to place in the tube.
        The tube is a list of tube symbols.

        Recognized keys in the tube are: 'call', 'parameters', 'reverse', 'cell'.
        The 'call' value should correspond to a mask_element function name or None.
        If 'call' is None the 'cell' keyword should point to an existing cell.
        If 'call' is not None, the 'parameters' value should correspond to a
        valid dict that can be passed as the arguments to the function in 'call'.
        The 'reverse' keyword (default False id not provided)
        to swap the 'a0' and 'b0' connections (reverse place the element).

        Args:
            geo (list): list of dict describing how create a cell or to provide a cell directly.
            showpins (bool):
            name (bool): cell name
            xs (str): xsection
            arrow (bool): draw start and end pin arrows (default=True)
            kwargs (dict): additional parameters to pars to mask_elements
            put (bool); put cells (default=True), or return list of cell objects in the tube.
            raise_pins (bool): if True, put anchor pins inside the tube

        Returns:
            Cell: Tube element based on the provided elements
        """
        # TODO: if xs is not None: pars['xs'] = xs  # does this make sense as the solver uses xs specific input?
        if name is None:
            name = 'ic_tube'
        if kwargs is None:
            kwargs = {}

        if cfg.instantiate_tube is None:
            instantiate = self.instantiate
        else:
            instantiate = cfg.instantiate_tube
        with nd.Cell(name=name, instantiate=instantiate, cnt=True) as ICcell:
            ICcell.group = 'interconnect'
            trace.trace_start()
            cells = []
            swapio = []  # bool list to swap 'b0' and 'a0' connections of elms
            flips = []  # bool list of flip state of elms
            skip_empty = []  # list to store which components are skipped for having too short lengths
            for elm in geo:
                if elm is None:
                    continue

                call = elm['call']
                reverse = elm.get('reverse', False)
                flip = elm.get('flip', False)

                # handle case element is a cell:
                if call is None:
                    cell = elm.get('cell', None)
                    if cell is not None:
                        cells.append(cell)
                        swapio.append(reverse)
                        flips.append(flip)
                    continue

                # handle call case
                valid = True
                pars = elm.get('parameters', {})
                if elm['call'] == 'strt':
                    if abs(pars['length']) < 1.e-5:
                        skip_empty.append(pars)
                        continue
                    cells.append(self._line(**pars, **kwargs))
                elif elm['call'] == 'bend':
                    if abs(pars['angle']) < 1.e-5:
                        skip_empty.append(pars)
                        continue
                    cells.append(self.Farc(**pars))
                elif elm['call'] == 'arc':
                    if abs(pars['angle']) < 1.e-5:
                        skip_empty.append(pars)
                        continue
                    cells.append(self.Farc_static(**pars))
                elif elm['call'] == 'ptaper':
                    if abs(pars['length']) < 1.e-5:
                        skip_empty.append(pars)
                        continue
                    cells.append(self._ptaper(**pars))
                elif elm['call'] == 'taper':
                    if abs(pars['length']) < 1.e-5:
                        skip_empty.append(pars)
                        continue
                    cells.append(self._taper(**pars))
                elif elm['call'] == "cobra":
                    cells.append(self._cobra(**pars))
                    #ICcell.Rmin = e1.cell.properties['Rmin']
                else:
                    valid = False
                    cfg.interconnect_errcnt += 1
                    msg = "IC-{} Tube element not recognized '{}' in cell '{}'".\
                        format(cfg.interconnect_errcnt, elm['call'], cfg.cells[-1].cell_name)
                    interconnect_logger(msg, 'error')
                    ic_exception(msg)
                if valid:
                    swapio.append(reverse)
                    flips.append(flip)

            if not put:
                return cells
            insts = []  # store instances for tube pins
            for i, cell in enumerate(cells[:]):
                if i == 0:  # input pin
                    inst = cells[0].put(flip=flips[i])
                    cp.push()
                    pinin = 'b0' if swapio[i] else 'a0'
                    pin = inst.pin[pinin]
                    nd.Pin('a0', io=0, pin=pin).put()
                    cp.pop()
                else:
                    if not swapio[i]:
                        inst = cell.put(flip=flips[i])
                    else:
                        inst = cell.put('b0', cp='a0', flip=(True != flips[i]))
                insts.append(inst)
            if cells:  # output pin
                pinout = 'a0' if swapio[i] else 'b0'
                pin = inst.pin[pinout]
                pinflip = not nd.get_xsection(pin.xs).symmetry
                nd.Pin('b0', io=1, pin=pin, flip=pinflip).put()
            else:  # nothing to place. TODO: Handle in cell.put(), where empty cells can be skipped, e.g. with cell.void = True
                ICcell.void = True
                #p1 = nd.Pin('a0', io=0, width=None, xs=xs).put(0, 0, 180)
                #p2 = nd.Pin('b0', io=1, width=None, xs=xs).put(0, 0, 0)
                p1 = nd.Pin('a0', io=0, width=skip_empty[0]['width'], xs=skip_empty[0]['xs']).put(0, 0, 180)
                p2 = nd.Pin('b0', io=1, width=skip_empty[-1]['width'], xs=skip_empty[-1]['xs']).put(0, 0, 0)
                nd.connect_path(p1, p2, None)

            if tubepins and len(insts) > 1:  # do not connect these tube pins the the netlist (xs=None)
                innerpin = 'a0' if swapio[0] else 'b0'
                outpin = 'b0' if swapio[0] else 'a0'
                nd.Pin(name=f"{cfg.tubepin_prefix}0", xs='tube', type='tube', io=-1, pin=insts[0].pin[outpin].rot(180)).put(drc=False)
                nd.Pin(name=f"{cfg.tubepin_prefix}1", xs='tube', type='tube', io=-1, pin=insts[0].pin[innerpin]).put(drc=False)

                outpin = 'a0' if swapio[-1] else 'b0'
                nd.Pin(f"{cfg.tubepin_prefix}{i+1}", pin=insts[-1].pin[outpin], xs='tube', type='tube', io=-1).put(drc=False)
                for i, inst in enumerate(insts[1:-1]):
                    innerpin = 'a0' if swapio[i+1] else 'b0'
                    nd.Pin(f'{cfg.tubepin_prefix}{i+2}', xs='tube', type='tube', io=-1, pin=inst.pin[innerpin]).put(drc=False)

            trace.trace_stop()
            ICcell.length_geo = trace.trace_length()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'], flip=ICcell.pin['b0'].flip)

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


if __name__ == "__main__":
    #nd.nazca_raise(2)
    #nd.cfg.instantiate_all()

    nd.cfg_groupconnect = True
    N = 15
    nd.cfg.autoconnectdistance = 0.002

    XS = nd.add_xsection('test')
    XS.minimum_radius = 10
    ic = Interconnect(xs='test', radius=10, width=1.0, pitch=5.0)
    A = 40
    Ea = 8

    ic.strt(N=4, length=70, pitch=3).put(1000, -300)

    rib = RibIt(pitch=[0, 5, 15, 30, 40], width=[1, 3, 2, 10, 10])
    ic.strt(length=100, config=rib).put(1000, -400)
    ic.bend(config=rib).put()
    ic.strt(length=100, config=rib).put()
    ic.bend(config=rib).put()
    ic.strt(length=150, config=rib).put()
    ic.bend(angle=-180, config=rib).put()

    rib = RibIt(pitch=20, width=[1, 3, 2, 10])
    ic.bend(angle=180, config=rib).put(1200, -400)

    rib2 = RibIt(pitch=[0, 5, 15, 30], width=3)
    ic.bend(angle=180, config=rib2).put(1400, -400)
    ic.sbend(offset=100, config=rib2).put()

    rib3 = rib.add(rib2, offset=10)
    ic.bend(angle=180, config=rib3).put(1600, -400)
    ic.sbend(offset=100, config=rib3).put()

    rib3 = RibIt(width=[1, 1, 1, 1, 15, 20, 10], gap=2)
    ic.bend(angle=160, config=rib3).put(1800, -400)

    rev = rib3.reverse()
    ic.bend(angle=160, config=rev).put(1900, -400)

    ic.ubend(offset=200, config=rib3).put(2000, -400)
    ic.ubend(offset=120, config=rev).put(2100, -400)
    ic.ubend(offset=120, config=rev, at1=1).put(2250, -400)
    ic.ubend(offset=120, config=rev, at1=-1).put(2350, -400)
    ic.ubend(offset=3, config=rev, at1=-1).put(2500, -400)
    ic.ubend(offset=-12, config=rev, at1=0).put(2650, -400, 45)

    ic.cobra_p2p(pin1=(2600, 0, 90), pin2=(2100, 0, 90), config=rev).put()
    ic.cobra_p2p(pin1=(2600, 0, -90), pin2=(2100, 0, -90), config=rib3, at1=-1, at2=-1).put()


    # check virtual radius
    Rvirtual = ic.euler_arc_euler_solve(eulerangle=Ea, get_virtual_radius=True)
    print('Rarc:', Rvirtual)
    for i, an in enumerate(np.linspace(-70, 70, 11)):
        ic.bend(angle=an, radius=Rvirtual, N=1).put(-1000, i*25)
        ic.euler_arc_euler(
            width=2.0, width2=1.0, eulerangle=Ea, angle=an, N=1, virtual_radius=True
        ).put(-1000, i*25)

        # euler_arc
        ic.euler_arc(width=3.0, width2=0.5, angle=3*an, radius_cal=20, angle_cal=30).put(-400, 1000 + 10*i)
        ic.euler_arc(width=3.0, width2=0.5, eulerangle=15, angle=3*an, radius_cal=20, angle_cal=30).put(-200, 1000 + 10*i)

    # set virtual radius
    ic.euler_arc_euler(width=3.0, width2=0.5, eulerangle=15, virtual_radius=300, N=1, angle=-A).put(-200, 200)

    # check angle = 0
    ic.euler_arc_euler(angle=0).put(0)

    ic.euler_arc_euler(N=5, angle=A).put(-800)

    # check euler angle + negative angle
    ic.euler_arc_euler(width2=2.0, eulerangle=Ea, N=2, angle=-A).put(-600)
    ic.euler_arc_euler(width2=2.0, eulerangle=Ea, N=2, angle=A).put(-600)

    # width taper
    ic.euler_arc_euler(width=2.0, width2=1.0, eulerangle=Ea, angle=A, N=3).put(-400)

    # recalibrate
    ic.euler_arc_euler(width=3.0, width2=0.5, eulerangle=15, radius_cal=20, angle_cal=30, N=4, angle=-A).put(-200)


    # one pin interconnects
    with nd.Cell('test_ribbon1') as C1:
        ic = Interconnect(N=N, radius=10, width=1.0, pitch=5.0)
        s1 = ic.strt().put(0)
        sb1 = ic.sbend(offset=-5, at1=0).put()
        ic.sbend(offset=100).put()

        ic.bend(angle=150).put()
        ic.ptaper(length=50, width2=2.5).put()
        ic.taper(length=50, width1=2.5).put()

        ic.strt(length=5).put()
        ic.bend(angle=-45).put()
        r = ic.bend(angle=30).put()

        ic.bend(N=2, angle=-90).put(cp='b1')
        ic.strt(N=1, length=40).put()

        b = ic.bend(N=N-2, angle=90).put(r.pin['b2'])
        ic.bend(N=3, angle=-90).put()
        ic.rot2ref(N=3, angle=90, at1=1).put()
        ic.strt(N=3, length=40).put()

        ic.sbend_p2p(N=5, pin1=b.pin['b12'], at1=4, pin2=b.pin['b4'].move(80, 40), Lstart=40).put()

        r2 = ic.bend_strt_bend_p2p(N=20, at1=2, pin2=s1.pin['a1'], at2=10).put()

        s2 = ic.strt(N=8).put(250, -250, -90)
        ic.strt_bend_strt_p2p(N=4, pin1=r2.pin['b0'], at1=0, pin2=s2.pin['a2'], at2=0).put()

        s3 = ic.strt(N=20).put(300, -250, -90)
        ic.ubend_p2p(N=5, pin1=s2.pin['b1'], at1=5, pin2=s3.pin['b5'], at2=0, length=100).put()

        ic.ubend_p2p(N=2, pin1=s3.pin['a0'], at1=0, pin2=s3.pin['a5'], at2=0, length=10).put()

        #nd.findpath(start=s1.pin['a5'], log=True)
        #nd.findpath(start=s1.pin['a0'], log=True)
        nd.Pin('start', pin=s1.pin['a5']).put()
        for i in range(N):
            paths = nd.pathfinder.Paths(startpin=s1.pin[f'a{i}'], show=True)

    # p2p interconnects
    with nd.Cell('test_ribbon2') as C2:
        N = 5
        ic = Interconnect(N=N, radius=30, width=1.0, pitch=5.0)
        s1 = ic.strt(N=1).put(0, 10, 0)
        at1 = 1
        at2 = 3
        ic.sbend(offset=100, at1=at1, at2=at2).put()
        ic.strt(N=1).put()

        ic.strt(at1=at1, at2=at2).put(s1.pin['b0'])
        ic.strt(N=1).put()

        ic.bend(at1=at1, at2=at2).put(s1.pin['b0'])
        ic.strt(N=1).put()

        ic.taper(at1=at1, at2=at2).put(s1.pin['b0'])
        ic.strt(N=1).put()

        ic.ptaper(at1=at1, at2=at2).put(s1.pin['b0'])
        ic.strt(N=1).put()

        # p2p
        s1 = ic.strt(N=1).put(120, 10, 0)
        s2 = ic.strt(N=1).put(220, 150, 0)

        ic.sbend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['a0'], at2=at2, tubepins=True).put()
        ic.sbend_p2p(N=1, pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['a0'], at2=at2).put()

        ic.bend_strt_bend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['a0'], at2=at2).put()
        ic.bend_strt_bend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['b0'], at2=at2).put()

        ic.ubend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['b0'], at2=at2, length=100).put()
        ic.ubend_p2p(pin1=s1.pin['b0'], pin2=s2.pin['b0'], length=100, N=1).put()

        s2 = ic.strt(N=1).put(220, -100, -90)
        ic.strt_bend_strt_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['a0'], at2=at2).put()

        # ubend test:
        s1 = ic.strt(N=1).put(750, 0, 0)
        s2 = ic.strt(N=1).put(800, 25, 0)
        ic.ubend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['b0'], at2=at2, length=100).put()

        s1 = ic.strt(N=1).put(850, 150, 0)
        s2 = ic.strt(N=1).put(900, 125, 0)
        ic.ubend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['b0'], at2=at2, length=100).put()

        s1 = ic.strt(N=1).put(950, 200, 0)
        s2 = ic.strt(N=1).put(1000, 450, 0)
        ic.ubend_p2p(pin1=s1.pin['b0'], at1=at1, pin2=s2.pin['b0'], at2=at2, length=100).put()

        ic.strt_p2p(pin1=(200, 600, 0), pin2=(-100, 450, -90), at1=at1, at2=at2).put()

        ic.cobra_p2p(pin1=(200, 700, 90), pin2=(-100, 550, -45), at1=at1, at2=at2).put()
        ic.strt(length=200, at1=at2).put()


    with nd.Cell('test_ribbon3') as C3:
        ic = nd.interconnects.Interconnect(radius=10, width=2.0, N=10, pitch=5.0, mambalayer=100)
        ic.sectionangle = 120
        points = [(0, 0), (100, 100), (100, -50), (200, -20), (210, -60), (500, -50)]
        ic.mamba(points=points, showpins=True, at1=2).put('org', 0)

    C1.put(0)
    C2.put(300)
    C3.put(500, 500)

    #nd.findpath(start=C.pin['start'], log=True)
    nd.export_gds(hierarchy='full')

    sb = ic.sbend_p2p(pin1=(0, 0, 0), at1=at1, pin2=(500, 600, 180), at2=at2, tubepins=True)
    nd.export_gds(topcells=sb, bb=True, filename='sb')

