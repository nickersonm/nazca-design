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

# @author: Ronald Broeke and Xaveer Leijtens (c) 2016-2017
# @email: ronald.broeke@brightphotonics.eu
#-----------------------------------------------------------------------

"""
This module defines mask elements such as straight and bend waveguides
via template functions.

Functions in mask_elements work with closures on nested functions. They
are templates functions and not specifically intended for use
in mask design directly by the designer. Instead use the interconnect
module and/or the predefined interconnects in
design kits.

However, a number of waveguide functions have been predefined in this
module for quick access and skipping Interconnect objects:

* strt      = Tp_straight()
* bend      = Tp_arc()
* ptaper    = Tp_ptaper()
* taper     = Tp_taper()
* cobra     = Tp_cobra()
* sinebend  = Tp_sinecurve()
"""

from math import sin, cos, acos, atan, atan2, degrees, radians, pi, sqrt, hypot, tan
from functools import partial
import numpy as np

from scipy.special import fresnel  # TODO: make scipy optional?

from nazca.netlist import Cell
from nazca import cfg
from nazca.logging import logger
from nazca.gds_base import gds_db_unit as gridsize
import nazca as nd
from nazca.util import polyline_length
from nazca.generic_bend import curve2polyline, gb_point, gb_coefficients, sinebend_point
from nazca.simglobal import sim


# TODO: move to cfg.py
min_length = 1e-6  # minimum straight waveguide length to draw in um
min_angle  = 1e-6   # minimume angle in rad to draw

def layeriter(xs=None, layer=None):
    """Generator yielding all layers in a xsection.

    Args:
        xs (str): xsection name
        layer (int): layer number

    Yields:
        layer, growx, growy, accuracy: iterate over all layers in <xs> and <layer>
    """
    if xs is None and layer is None:
        xs = cfg.default_xs_name
        if xs not in cfg.XSdict.keys():
            handle_missing_xs(xs)

    if layer is not None:
        grow = ((0.5, 0), (-0.5, 0), 0, 0)  # a layer not as part of a xs has no grow
        layer_name = nd.get_layer(layer)  # cfg.layer_names[layer][0:2]
        lineitem = cfg.layer_table.loc[layer_name]
        try:
            accuracy = lineitem['accuracy']
        except Exception as E:
            print("Error: accuracy not defined for layer {}.".format(layer))
            print(E)
            accuracy = 0.001
        yield layer_name, grow, accuracy, False

    if xs is not None:
        XS = cfg.XSdict.get(xs, None)
        if XS is None:
            layer = handle_missing_xs(xs)
        ML = cfg.XSdict[xs].mask_layers
        if ML.empty:
            msg = "xsection '{0}' has no layers. "\
                "Continuing by adding fallback layer '{1}' to '{0}'.\n"\
                "Recommended solutions: Use a different xsection "\
                "or add layers to this xsection:\n"\
                "add_layer2xsection('{0}', layer=<num>).".\
                format(xs, cfg.default_layers['dump'])
            if cfg.redirect_unknown_layers:
                nd.main_logger(msg, "error")
            else:
                nd.main_logger(msg, "warning")
                nd.add_layer2xsection(xsection=xs, layer=cfg.default_layers['dump'])
            yield cfg.default_layers['dump'], ((0.5, 0), (-0.5, 0), 0, 0), 0.1, False
        else:
            for layer_name, A in ML.iterrows():
                grow = (
                    (A.leftedgefactor, A.leftedgeoffset),
                    (A.rightedgefactor, A.rightedgeoffset),
                    A.growy1,
                    A.growy2
                )
                yield layer_name, grow, A.accuracy, A.polyline


def handle_missing_xs(xs):
    """Handle missing xsections or xsection missing mask_layers.

    If xsection <xs> is not defined, create a xsection with layer.
    A warning is issued, expect if the Nazca default xsection 'nazca' is created.

    If the xsection exists but has not mask_layer attribute this function
    will add a default layer to the xsection.
    The xsextions 'dump' and 'error' are special cases.

    Returns:
        layer: layer if any has been created in this method.
    """
    layer = None
    if xs not in cfg.XSdict.keys():
        if xs != cfg.default_xs_name and xs != cfg.default_xserror_name:
            logger.warning("No xsection named '{0}'. "\
                "Correct the name or add the xsection with:\n"\
                "add_xsection('{0}') to get rid of this warning. "\
                "Already available xsections are {1}.\n".\
                format(xs, list(cfg.XSdict.keys())))
        nd.add_xsection(xs)

    add_layer = False
    if not hasattr(nd.get_xsection(xs), 'mask_layers'):
        add_layer = True
    elif nd.get_xsection(xs).mask_layers is None:
        add_layer = True
    elif nd.get_xsection(xs).mask_layers.empty:
        add_layer = True
    if add_layer:
        if xs == cfg.default_xserror_name:
            layer = 'error'
        else:
            layer = 'dump'
        layer = nd.get_layer(layer)
        nd.add_layer2xsection(xs, layer=layer)
    return layer


#==============================================================================
# Waveguide element definitions
#==============================================================================
cnt = 0  # ordinal counter for unique naming
def Tp_straight(
    length=10,
    width=1.0,
    xs=None,
    layer=None,
    edge1=None,
    edge2=None,
    edgepoints=50,
    name=None,
    modes=None,
):
    """Template for creating parametrized straight waveguide function.

    Args:
        length (float): length of waveguide
        width (float): width of waveguide
        xs (str): xsection of taper
        layer (int | str): layer number or layername
        edge1 (function): optional function F(t) describing edge1 of the waveguide
        edge2 (function): optional function G(t) describing edge2 of the waveguide
        edgepoints (int): number of edge point per edge if the edge is set (default=50)
        modes (list): list of integers, containing the labels of netlist modes

    Returns:
        function: Function returning a Cell object with a straight guide
    """
    if modes is None:
        modes = sim.modes

    def cell(
        length=length,
        width=width,
        xs=xs,
        layer=layer,
        edge1=edge1,
        edge2=edge2,
        edgepoints=edgepoints,
        name=name,
        gridpatch=0,
    ):
        """Create a straight waveguide element.

        The edge option allow for defining parametrized edges as a perturbation
        of the standard edges. If the edges have to interpreted as an absolute
        value, then set width=0. Note that if edge1 function is defined and
        edge2 is not then edge2 is -edge1.

        Args:
            length (float): length of waveguide
            width (float): width of waveguide
            xs (str): xsection of waveguide
            layer (int | str): layer number or layer name
            edge1 (function): optional function F(t) describing edge1 of the waveguide
            edge2 (function): optional function G(t) describing edge2 of the waveguide
            edgepoints (int): number of edge points per edge if edge1 is set (default=50)

        Returns:
            Cell: straight element
        """
        if name is None:
            name = 'straight'
        # assert width is not None
        if width is None:
            width = 0.0

        if length < 0:
            nd.interconnect_logger(f"Negative straight waveguide length of {length}.", "error" )

        with Cell(name=name, cnt=True) as C:

            def CM(wl, pol):
                """Optical path length model."""
                try:
                    Neff = nd.get_xsection(xs).index.Neff
                    return Neff(width=width, wl=wl, pol=pol) * length
                except:
                    nd.main_logger(
                        f"No index model found in strt for xs = '{xs}' in cell '{C.cell_name}'.",
                        "error",
                    )
                return None

            C.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(name='a0', io=0, width=width, xs=xs, show=True).put(0, 0, 180)
            p2 = nd.Pin(name='b0', io=1, width=width, xs=xs, show=True).put(length, 0, 0)

            nd.connect_path(p1, p2, length)
            for i in modes:
                nd.connect_path(p1, p2, CM, sigtype=f'm{i}')

            C.updk = {
                'call': 'strt',
                'parameters': {
                    'length': {'value': length , 'unit': 'um', 'type': 'float'},
                    'width': {'value': width, 'unit': 'um', 'type': 'float'},
                },
            }

            # to be removed:
            C.length_geo = length # TODO: used by old trace method: to be removed
            C.properties['parameters'] = {
                'length': {'value': length , 'unit': 'um', 'type': 'float'},
                'width': {'value': width, 'unit': 'um', 'type': 'float'},
                'xs': xs,
                'call': 'strt',
            }

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow                    
                sign = 1.0
                if not polyline:
                    if edge1 is None:
                        if gridpatch == 0:  # reproduce polygon before gridpatch
                            outline = [
                                (0 - c1, width * a1 + b1), 
                                (length + c2, width * a1 + b1),
                                (length + c2, width * a2 + b2),
                                (0 - c1, width * a2 + b2)
                            ]
                        else:
                            outline = [
                                (0 - c1 - 2*gridpatch, width * a1 + b1 - gridpatch), 
                                (0 - c1, width * a1 + b1), 
                                (length + c2, width * a1 + b1),
                                (length + c2 + 2*gridpatch, width * a1 + b1 - gridpatch),
                                (length + c2 + 2*gridpatch, width * a2 + b2 + gridpatch), 
                                (length + c2, width * a2 + b2), 
                                (0 - c1, width * a2 + b2),
                                (0 - c1 - 2*gridpatch, width * a2 + b2 + gridpatch), 
                              ]
                    else:
                        if edge2 is None:
                            edge2 = edge1
                            sign = -1.0
                        Fp1 = []
                        Fp2 = []
                        for t in np.linspace(0, 1, edgepoints):
                            Fp1.append((length*t, width*a1 + b1 + edge1(t)))
                            Fp2.append((length*t, width*a2 + b2 + sign*edge2(t)))
                        outline = Fp1 + list(reversed(Fp2))
                    if abs(length + c1 + c2) > min_length:
                        nd.Polygon(layer=lay, points=outline).put(0)
                else:  # polyline
                    if abs(length + c1 + c2) > min_length:
                        wpoly = width*(a1-a2) + (b1-b2)
                        centre = 0.5*(width*(a1+a2) + (b1+b2))
                        centreline = [(0-c1, centre), (length+c2, centre)]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C
    return cell


def Tp_taper(
    length=100,
    width1=2,
    width2=3,
    shift=0,
    xs=None,
    layer=None,
    name=None,
    modes=None,
):
    """Template for creating a parametrized linear taper function.

    Note that zero length taper segments may seem to overlook width pin2pin DRC.

    Args:
        length (float): length of the taper
        width1 (float): width at start
        width2 (float): width at end
        shift (float): lateral shift at taper end resulting in a skew taper
        xs (str): xsection of taper
        layer (int | str): layer number or layer name

    Returns:
        function: Function returning a Cell object with a linear taper
    """
    if modes is None:
        modes = sim.modes

    def cell(
        length=length,
        width1=width1,
        width2=width2,
        shift=shift,
        xs=xs,
        layer=layer,
        name=name,
    ):
        """Create a taper element.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            shift (float): lateral shift at taper end resulting in a skew taper
            xs (str): xsection of taper
            layer (int | str): layer number or layer name

        Returns:
            Cell: taper element
        """
        if name is None:
            name = 'taper'
        with Cell(name=name, cnt=True) as C:
            # l, w for geom/optical length/width
            l, w = length, (width1 + width2) / 2
            if shift != 0:
                l = sqrt(length ** 2 + shift ** 2)
                # projected average width
                w = w * cos(atan(shift / length))
            def CM(wl, pol):
                """Optical path length model."""
                try:
                    Neff = nd.get_xsection(xs).index.Neff
                    n1 = Neff(width=width1, wl=wl, pol=pol)
                    n2 = Neff(width=width2, wl=wl, pol=pol)
                    return 0.5* (n1+n2) * length
                except:
                    nd.main_logger(
                        f"No index model found in strt for xs = '{xs}' in cell '{C.cell_name}'.",
                        "error",
                    )
                return None

            C.instantiate = cfg.instantiate_mask_element
            if abs(length) < gridsize / 2:
                return C # Empty taper
            if width1 > width2:
                swap = True
                pin = ['b0', 'a0']
                width1, width2 = width2, width1
            else:
                swap = False
                pin = ['a0', 'b0']

            p1 = nd.Pin(name=pin[0], io=0, width=width1, xs=xs, show=True).put(0, 0, 180)
            p2 = nd.Pin(name=pin[1], io=1, width=width2, xs=xs, show=True).put(
                length,
                shift,
                0
            )

            nd.connect_path(p1, p2, l)
            for i in modes:
                nd.connect_path(p1, p2, CM, sigtype=f'm{i}')
            C.updk= {
                'call': 'taper',
                'parameters': {
                    'length': {'value': length, 'unit': 'um', 'type': 'float'},
                    'width1': {'value': width1, 'unit': 'um', 'type': 'float'},
                    'width2': {'value': width2, 'unit': 'um', 'type': 'float'},
                    'shift': {'value': shift, 'unit': 'um', 'type': 'float'},
                }
            }

            # to remove:
            C.length_geo = l
            C.properties['parameters'] = {
                'length': {'value': length, 'unit': 'um', 'type': 'float'},
                'width1': {'value': width1, 'unit': 'um', 'type': 'float'},
                'width2': {'value': width2, 'unit': 'um', 'type': 'float'},
                'shift': {'value': shift, 'unit': 'um', 'type': 'float'},
                'xs': xs,
                'call': 'taper',
            }

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if swap:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # Set widths for this layer
                twidth1a = a1*width1 + b1
                twidth1b = a2*width1 + b2
                twidth2a = a1*width2 + b1 + shift
                twidth2b = a2*width2 + b2 + shift
                if not polyline:
                    outline = [
                        (0-c1, twidth1a), (length+c2, twidth2a),
                        (length+c2, twidth2b), (0-c1, twidth1b)
                    ]
                    nd.Polygon(layer=lay, points=outline).put()
                else:
                    if abs(l) > min_length:
                        wpoly = twidth1a - twidth1b
                        centre = 0.5*(twidth1a + twidth1b)
                        centreline = [(0-c1, centre), (length+c2, centre + shift)]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C
    return cell


def Tp_ptaper(
    length=100,
    width1=1.0,
    width2=3.0,
    xs=None,
    layer=None,
    name=None,
    modes=None,
):
    """Template for creating a parametrized parabolic taper  function.

    Note that zero length taper segments may seem to overlook width pin2pin DRC.

    Args:
        length (float): length of the taper
        width1 (float): width at start
        width2 (float): width at end
        xs (str): xsection of taper
        layer (int | str): layer number or layer name

    Returns:
        function: Function returning a Cell object with a ptaper
    """
    if modes is None:
        modes = sim.modes

    def cell(
        length=length,
        width1=width1,
        width2=width2,
        xs=xs,
        layer=layer,
        name=None,
    ):
        """Create a parabolic taper element.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            xs (str): xsection of taper
            layer (int | str): layer number or layer name

        Returns:
            Cell: parabolic taper element
        """
        # TODO: add growy
        if name is None:
            name = 'ptaper'
        with Cell(name=name, cnt=True) as C:
            def CM(wl, pol):
                """Optical path length model.

                Estimate by taking index of mid width.
                """
                try:
                    Neff = nd.get_xsection(xs).index.Neff
                    n1 = Neff(width=width1, wl=wl, pol=pol)
                    n2 = Neff(width=width2, wl=wl, pol=pol)
                    return 0.5* (n1+n2) * length
                except:
                    nd.main_logger(
                        f"No index model found in strt for xs = '{xs}' in cell '{C.cell_name}'.",
                        "error",
                    )
                return None

            C.instantiate = cfg.instantiate_mask_element
            if abs(length) < gridsize / 2:
                return C # Empty taper
            if width1 > width2:
                swap = True
                pin = ['b0', 'a0']
                width1, width2 = width2, width1
            else:
                swap = False
                pin = ['a0', 'b0']

            p1 = nd.Pin(name=pin[0], io=0, width=width1, xs=xs, show=True).put(0, 0, 180)
            p2 = nd.Pin(name=pin[1], io=1, width=width2, xs=xs, show=True).put(length, 0, 0)

            nd.connect_path(p1, p2, length)
            for i in modes:
                nd.connect_path(p1, p2, CM, sigtype=f'm{i}')
            C.updk = {
                'call': 'taper',
                'parameters': {
                    'length': {'value': length, 'unit': 'um', 'type': 'float'},
                    'width1': {'value': width1, 'unit': 'um', 'type': 'float'},
                    'width2': {'value': width2, 'unit': 'um', 'type': 'float'},
                }
            }

            # to remove:
            C.length_geo = length
            C.properties['parameters'] = {
                'length': {'value': length ,'unit': 'um', 'type': 'float'},
                'width1': {'value': width1, 'unit': 'um', 'type': 'float'},
                'width2': {'value': width2, 'unit': 'um', 'type': 'float'},
                'xs': xs,
                'call': 'ptaper',
            }

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if swap:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                if not polyline:
                    # Set widths for this layer
                    twidth1a = a1*width1 + b1
                    twidth1b = a2*width1 + b2
                    twidth2a = a1*width2 + b1
                    twidth2b = a2*width2 + b2
                    if abs(0.5*width1 - 0.5*width2) / (4*length) < gridsize:
                        # No shaping required due to small delta width:
                        points = [
                            (0, twidth1a), (length, twidth2a),
                            (length, twidth2b), (0, twidth1b)
                        ]
                        nd.Polygon(layer=lay, points=points).put(0)
                        continue
                    # y = a * x**2
                    a = 4 * length / (width2**2 - width1**2)
                    y2 = y0 = a * (0.5*width1)**2
                    w_tap1 = width1
                    ptop = [(0, twidth1a)]
                    pbot = [(0, twidth1b)]
                    while y2 - y0 < length:
                        w_tap2 = w_tap1 + 4*acc + 4*sqrt(acc*(w_tap1 + acc))
                        y2 = a * (0.5*w_tap2)**2
                        if y2 - y0 < length:
                            ptop.append((y2-y0, w_tap2*a1+b1))
                            pbot.append((y2-y0, w_tap2*a2+b2))
                            w_tap1 = w_tap2
                        else:
                            break
                    ptop.append((length, twidth2a))
                    pbot.append((length, twidth2b))
                    nd.Polygon(layer=lay, points=ptop + list(reversed(pbot))).put(0)
                else:
                    if abs(length) > min_length:
                        wpoly = width1*(a1-a2)+(b1-b2)
                        centre = 0.5*(width1*(a1+a2)+(b1+b2))
                        centreline = [(0-c1, centre), (length+c2, centre)]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C
    return cell


def __get_offset(xs, width, radius, offset=None):
    """Get the offset function."""
    if offset is None:
        offset = 0
        if xs is None or abs(radius) < 1e-3:
            return offset
        try:
            F = cfg.XSdict[xs].os
        except:
            handle_missing_xs(xs)
            cfg.XSdict[xs].os = 0
        else:
            try:
                offset = float(F)
            except TypeError:
                offset = F(width=width, radius=radius)
        return offset
    else:
        try:
            os = float(offset)
        except TypeError:
            os = offset(width=width, radius=radius)
        return os


def Tp_arc(
    radius=10,
    width=1.0,
    angle=90,
    xs=None,
    layer=None,
    offset=None,
    name=None,
    modes=None,
):
    """Template for creating a parametrized circular arc waveguide function.

    Args:
        radius (float): radius at the center line of the arc in um.
        width (float): width of the arc in um.
        angle (float): angle of arc in degree (default = 90).
        xs (str): xsection of taper
        layer (int | str): layer number or layer name
        offset (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        modes (list): list of integers, containing the labels of netlist modes

    Returns:
        function: Function returning a Cell object with an arc
    """
    if modes is None:
        modes = sim.modes

    # TODO: add growy
    def cell(
        radius=radius,
        width=width,
        angle=angle,
        xs=xs,
        layer=layer,
        offset=offset,
        name=name,
    ):
        """Create a circular arc element.

        A straight-bend offset is included when it has been defined in the
        xsection used.

        Args:
            radius (float): radius at the center line of the arc in um.
            width (float): width of the arc in um.
            angle (float): angle of arc in degree (default = 90).
            xs (str): xsection of taper
            layer (int | str): layer number or layer name
            offset (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float

        Returns:
            Cell: circular arc element
        """
        if name is None:
            name = 'arc'
        ang = np.radians(angle)
        if abs(ang) < min_angle:
            ang = 0
            angle = 0
        sign = np.sign(ang)
        radius = abs(radius)

        #if offset is None:
        offset = __get_offset(xs, width, radius, offset=offset)

        Nmax = cfg.maxpolygonpoints // 2 - 2  # 2xnmax + 4
        with Cell(name=name, cnt=True) as C:
            def CM(wl, pol):
                """Optical path length model."""
                try:
                    Neff = nd.get_xsection(xs).index.Neff
                    return Neff(width=width, radius=radius, wl=wl, pol=pol) * abs(radius*ang)
                except:
                    nd.main_logger(
                        f"No index model found in strt for xs = '{xs}' in cell '{C.cell_name}'.",
                        "error"
                    )
                    return None

            #C.use_hull = True
            C.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(name='a0', io=0, width=width, xs=xs, radius=radius, show=True).put(0, 0, 180)
            p2 = nd.Pin(name='b0', io=1, width=width, xs=xs, radius=radius, show=True).\
                put(radius*sin(abs(ang)), sign*radius*(1-cos(ang)), angle)

            nd.connect_path(p1, p2, abs(radius*ang))
            for i in modes:
                nd.connect_path(p1, p2, CM, sigtype=f'm{i}')
            C.updk = {
                'call': 'bend',
                'parameters': {
                    'angle':  {'value': angle ,'unit': 'deg', 'type': 'float'},
                    'width':  {'value': width, 'unit': 'um', 'type': 'float'},
                    'radius': {'value': radius, 'unit': 'um', 'type': 'float'},
                }
            }

            # to be removed:
            C.length_geo = abs(radius*ang)
            C.properties['parameters'] = {
                'angle':  {'value': angle ,'unit': 'deg', 'type': 'float'},
                'width':  {'value': width, 'unit': 'um', 'type': 'float'},
                'radius': {'value': radius, 'unit': 'um', 'type': 'float'},
                'xs':  xs,
                'call':  'bend',
            }

            if ang == 0:
                return C

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if angle >= 0:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # Start and end sections are different from other sections
                if polyline:
                    wpoly = width*(a1-a2) + (b1-b2)
                    Ro = Ri = radius +0.5*(width*(a1+a2) + (b1+b2)) - offset
                else:
                    Ro = radius + width * a1 - offset + b1 # outer radius
                    Ri = radius + width * a2 - offset + b2 # inner radius
                    #TODO: add growy
                if Ri < 0:
                    if Ri < -1e-6:
                        nd.main_logger("Side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f}, xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".\
                            format(cfg.cells[-1].cell_name, Ri, radius, xs, width, offset),
                            "warning")
                    Ri = 0
                if Ro < 0:
                    if Ro < -1e-6:
                        nd.main_logger("Side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f}, xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".\
                            format(cfg.cells[-1].cell_name, Ro, radius, xs, width, offset),
                            "warning")
                    Ro = 0

                Rmax = max(Ri, Ro)
                if not isinstance(acc, float):
                    acc = 0.001
                    logger.error(f"Accucary not defined for layer {lay}. Setting it to {acc}.")
                da = 2 * acos(Rmax / (Rmax + acc)) # step angle for accuracy
                N = int((abs(ang) / da) + 1) # +1 for rounding
                da = ang / N # step angle for N
                u = np.linspace(da / 2, ang - da / 2, N)
                Roe = Ro * (da / 2) / sin(da / 2) # effective outer radius
                Rie = Ri * (da / 2) / sin(da / 2) # effective inner radius

                p1 = [(Roe * sin(abs(a)), sign * (radius - Roe * cos(a)))
                    for a in u]
                p2 = [(Rie * sin(np.abs(a)), sign * (radius - Rie * cos(a)))
                    for a in u]

                pstart = [(0, sign * (radius - Ri)), (0, sign * (radius - Ro))]
                pend = [
                    (Ro * sin(abs(ang)),
                    sign * (radius - Ro * cos(ang))),
                    (Ri * sin(abs(ang)),
                    sign * (radius - Ri * cos(ang)))
                ]

                section = list(range(0, len(p1), Nmax))
                if not polyline:
                    for s in section[:-1]:
                        outline = pstart\
                            + p1[s:s + Nmax]\
                            + list(reversed(p2[s:s + Nmax]))
                        nd.Polygon(layer=lay, points=outline).put(0)
                        pstart = [p2[s + Nmax - 1], p1[s + Nmax - 1]]

                    outline = pstart\
                         + p1[section[-1]:section[-1] + Nmax]\
                         + pend\
                         + list(reversed(p2[section[-1]:section[-1] + Nmax]))
                    nd.Polygon(layer=lay, points=outline).put(0)
                else:
                    for s in section[:-1]:
                        centreline = [pstart[0]] + p1[s:s + Nmax]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
                        pstart = [p2[s + Nmax - 1], p1[s + Nmax - 1]]

                    centreline = [pstart[0]]\
                        + p1[section[-1]:section[-1] + Nmax]\
                        + [pend[0]]
                    nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C
    return cell


def Tp_arc2(
    radius=10,
    width=1.0,
    angle=90,
    xs=None,
    layer=None,
    offset=None,
    name=None,
):
    """Template for creating a parametrized angled arc waveguide function.

    The arc is composed of straight sections.
    Minimum angle between sections is 45 deg.

    Args:
        radius (float): radius at the center line of the arc in um.
        width (float): width of the arc in um.
        angle (float): angle of arc in degree (default = 90).
        xs (str): xsection of taper
        layer (int | str): layer number or layer name
        offset (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float

    Returns:
        function: Function returning a Cell object with an arc composed of straight sections
    """
    # TODO: add growy
    def cell(
        radius=radius,
        width=width,
        angle=angle,
        xs=xs,
        layer=layer,
        offset=offset,
        name=name,
    ):
        """Create a single circular arc element omposed of straight sections.

        A straight-bend offset is included when it has been defined in the
        xsection used.

        Args:
            radius (float): radius at the center line of the arc in um.
            width (float): width of the arc in um.
            angle (float): angle of arc in degree (default = 90).
            xs (str): xsection of taper
            layer (int | str): layer number or layer name
            offset (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float

        Returns:
            Cell: circular arc element
        """
        if name is None:
            name = 'arc'
        ang = np.radians(angle)
        if abs(ang) < min_angle:
            ang = 0
            angle = 0
        sign = np.sign(ang)
        radius = abs(radius)

        #if offset is None:
        offset = __get_offset(xs, width, radius, offset=offset)

        Nmax = cfg.maxpolygonpoints // 2 - 2  # 2xnmax + 4
        da = pi/4.0
        Ni=int((abs(ang) / da)-1e-6)+1
        ang = ang / (1.0*Ni) # step angle for N
        with Cell(name='int', cnt=True) as INT:
            #C.use_hull = True
            INT.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(name='a0', io=0, width=width, xs=xs, show=True).put(0, 0, 180)
            p2 = nd.Pin(name='b0', io=1, width=width, xs=xs, show=True).\
                put(radius*sin(abs(ang)), sign*radius*(1-cos(ang)), 180.0/pi*ang)

            nd.connect_path(p1, p2, 2.0*abs(radius*tan(ang*0.5)))
            INT.length_geo = abs(radius*ang)

            if ang == 0:
                return INT

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if angle >= 0:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # Start and end sections are different from other sections
                if polyline:
                    wpoly = width*(a1-a2) + (b1-b2)
                    Ro = Ri = radius +0.5*(width*(a1+a2) + (b1+b2)) - offset
                else:
                    Ro = radius + width * a1 - offset + b1 # outer radius
                    Ri = radius + width * a2 - offset + b2 # inner radius
                    #TODO: add growy
                if Ri < 0:
                    if Ri < -1e-6:
                        nd.main_logger("Side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f}, xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".\
                            format(cfg.cells[-1].cell_name, Ri, radius, xs, width, offset),
                            "warning")
                    Ri = 0
                if Ro < 0:
                    if Ro < -1e-6:
                        nd.main_logger("Side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f}, xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".\
                            format(cfg.cells[-1].cell_name, Ro, radius, xs, width, offset),
                            "warning")
                    Ro = 0
                Rmax = max(Ri, Ro)
                Roe = Ro / cos(abs(ang)/2.0) # effective outer radius
                Rie = Ri / cos(abs(ang)/2.0) # effective inner radius
                p1 = [
                    (Roe * sin(abs(0.5*ang)),
                    sign * (radius - Roe * cos(0.5*ang)))]
                p2=  [
                    (Rie * sin(abs(0.5*ang)),
                    sign * (radius - Rie * cos(0.5*ang)))
                ]
                pstart = [(0, sign * (radius - Ri)), (0, sign * (radius - Ro))]
                pend = [
                    (Ro * sin(abs(ang)),
                    sign * (radius - Ro * cos(ang))),
                    (Ri * sin(abs(ang)),
                    sign * (radius - Ri * cos(ang)))
                ]
                if not polyline:
                    outline = pstart\
                         + p1\
                         + pend\
                         + list(reversed(p2))
                    nd.Polygon(layer=lay, points=outline).put(0)
                else:
                    centreline = [pstart[0]]\
                        + p1\
                        + [pend[0]]
                    nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)

        inst=[]
        with Cell(name='arc', cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            for i in range(Ni):
                inst.append(INT.put())
            nd.Pin('a0').put(inst[0].pin['a0'])
            nd.Pin('b0').put(inst[-1].pin['b0'])

        C.updk = {
            'parameters': {
                'angle':  {'value': angle ,'unit': 'deg', 'type': 'float'},
                'width':  {'value': width, 'unit': 'um', 'type': 'float'},
                'radius': {'value': radius, 'unit': 'um', 'type': 'float'},
            },
            'call':  'bend',
        }

        C.properties['parameters'] = {
            'angle':  {'value': angle ,'unit': 'deg', 'type': 'float'},
            'width':  {'value': width, 'unit': 'um', 'type': 'float'},
            'radius': {'value': radius, 'unit': 'um', 'type': 'float'},
            'xs':     {'value': xs, 'unit': None, 'type': 'str'},
            'call':   {'value': 'bend', 'unit': None, 'type': 'str'},
        }
        return C
    return cell



def Tp_sinecurve(
    width=1.0,
    distance=200,
    offset=20,
    xs=None,
    layer=None,
    name=None,
 ):
    """Template for creating parametrized sine curve waveguide function.

    Args:
        width (float): width of the interconnect in um
        pin (Node): optional Node for modeling info
        xs (str): xsection of sinebend
        distance (float): total forward length of the sinebend in um
        offset (float): lateral offset of the sinebend in um
        xs (str): xsection of waveguide
        layer (int | str): layer number or layer name

    Returns:
        function: Function returning a Cell object with the sinecurve guide
    """
    # TODO: add growy
    def cell(
        width=width,
        distance=distance,
        offset=offset,
        xs=xs,
        layer=layer,
        name=name
    ):
        """Create a (raised) sine bend waveguide element.

        Args:
            width (float): width of the interconnect in um
            pin (Node): optional Node for modeling info
            xs (str): xsection of sinebend
            offset (float): lateral offset of the sinebend in um
            distance (float): total forward length of the sinebend in um
            xs (str): xsection of waveguide
            layer (int | str): layer number or layer name

        Returns:
            Cell: sinecurve element
        """

        if name is None:
            name = 'sinecurve'

        xya = (distance, offset, 0) # End point
        # sinecurve waveguide, cwg
        with Cell(name=name, cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            nd.Pin(name='a0', io=0, width=width, xs=xs, show=True).put(0, 0, 180)
            nd.Pin(name='b0', io=1, width=width, xs=xs, show=True).put(*xya)
            # TODO: length_geo
            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # TODO: define polyline2polygon with left and right width from spine
                # sampled curve
                xy = curve2polyline(sinebend_point, xya, acc, (distance, offset))
                if not polyline:
                    # polygon of proper width
                    xy = nd.util.polyline2edge(xy, width, grow=grow)
                    nd.Polygon(layer=lay, points=xy).put(0)
                else:
                    wpoly = width*(a1-a2) + (b1-b2)
                    xy = nd.util.polyline2edge(xy, width, grow=grow, line=True)
                    nd.Polyline(layer=lay, points=xy, width=abs(wpoly)).put(0)

            C.updk = {
                 'call': 'sinecurve',
                 'parameters': {
                    'distance': {'value': distance ,'unit': 'um', 'type': 'float'},
                    'width': {'value': width, 'unit': 'um', 'type': 'float'},
                    'offset': {'value': offset, 'unit': 'um', 'type': 'float'},
                },
            }

            # to remove
            C.properties['parameters'] = {
                'distance':  {'value': distance ,'unit': 'um', 'type': 'float'},
                'width':  {'value': width, 'unit': 'um', 'type': 'float'},
                'offset': {'value': offset, 'unit': 'um', 'type': 'float'},
                'xs':     {'value': xs},
                'call':   {'value': 'sinecurve'},
            }
        return C
    return cell


def Tp_cobra(
    xya=(100, 100, 10),
    width1=1.0,
    width2=1.0,
    radius1=0,
    radius2=0,
    offset1=None,
    offset2=None,
    xs=None,
    layer=None,
    name=None,
):
    """Template for creating parametrized cobra waveguide function.

    Args:
        xya (point): point to connect to from (0,0,0)
        width1 (float): width of waveguide
        width2 (float): width of waveguide at end
        radius1 (float): radius at start (0 is no curvature)
        radius2 (float): radius at end (0 is no curvature)
        offset1 (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        offset2 (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        xs (str): xsection of waveguide
        layer (str | tuple | int): layer number or layer name

    Returns:
        function: Function returning a Cell object with the cobra guide
    """
    # TODO: add growy
    def cell(
        xya=xya,
        width1=width1,
        width2=width2,
        radius1=radius1,
        radius2=radius2,
        offset1=offset1,
        offset2=offset2,
        xs=xs,
        layer=layer,
        name=name,
    ):
        """Create a parametric waveguide element.

        Args:
            xya (point): point to connect to from (0,0,0)
            width (float): width of waveguide
            xs (str): xsection of waveguide
            layer (int | str): layer number or layer name

        Returns:
            Cell: cobra element
        """
        if name is None:
            name = 'cobra'
        if offset1 is None:
            offset1 = __get_offset(xs, width1, radius1)
        if radius1 != 0:
            radius1 = radius1 - offset1
        if offset2 is None:
            offset2 = __get_offset(xs, width2, radius2)
        if radius2 != 0:
            radius2 = radius2 - offset2
        if width2 is None:
            width2 = width1

        # cobra waveguide, pwg
        with Cell(name=name, cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(name='a0', io=0, width=width1, xs=xs, radius=radius1, show=True).\
                put(0, -offset1, 180)
            xya = (xya[0], xya[1]-offset1, xya[2])
            p2 = nd.Pin(name='b0', io=1, width=width2, xs=xs, radius=radius2, show=True).put(*xya)

            xya = (xya[0] - offset2*sin(np.radians(xya[2])), xya[1] +
                offset2*cos(np.radians(xya[2])), xya[2])
            # Solve the generic bend
            A, B, L, Rmin = gb_coefficients(xya, radius1=radius1,
                radius2=radius2)
            nd.connect_path(p1, p2, None)
            C.length_geo = None
            C.properties['Rmin'] = Rmin
            saved_acc = 10

            Rdrc = 1e8
            try:
                Rdrc = cfg.XSdict[xs].minimum_radius
                if Rdrc is None:
                    Rdrc = 0
            except KeyError:
                if xs not in cfg.XSdict.keys():
                    Rdrc = 0
            if Rdrc > Rmin:
                logger.warning('DRC minimum_radius {:.3f} < {:.3f}'.
                    format(Rmin, Rdrc))

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                # sampled curve
                xy = curve2polyline(gb_point, xya, acc, (A, B, L))
                if C.length_geo is None or saved_acc > acc:
                    # Calculate length of polyline and keep the one for the
                    # most accurate layer
                    C.length_geo = polyline_length(xy)
                    saved_acc = acc
                # polygon of proper width
                # TODO: add growy
                if not polyline:
                    if callable(width1):
                        def wfun(t):
                            return width1(t)
                        xy2 = nd.util.polyline2edge(
                            xy,
                            wfun,
                            grow=grow
                        )
                    else:
                        xy2 = nd.util.polyline2edge(
                            xy,
                            width1,
                            width2,
                            grow=grow
                        )
                        nd.Polygon(layer=lay, points=xy2).put(0)
                else:
                    wpoly = width1*(a1-a2) + (b1-b2)
                    xy = nd.util.polyline2edge(xy, width1, grow=grow, line=True)
                    nd.Polyline(layer=lay, points=xy, width=abs(wpoly)).put(0)
        return C
    return cell


def Tp_euler(
    width1=1.0,
    width2=1.0,
    radius=50,
    angle=90,
    xs=None,
    layer=None,
    name='euler',
    modes=None,
):
    """Template for creating an Euler bend.

    Call the template with which <radius> to achieve at a specified <angle>.
    The Euler function will then always follow that same shape,
    also for partial angles.

    Args:
        width1 (float): begin width
        width2 (float): end width
        radius (float): end radius for calibration
        angle (float): end angle for calibration
        xs (str): xsection na,e
        layer (str | tuple | int): mask layer
        name (str): element name (default='euler')

    Returns:
        function: Function returning a Cell object with the Euler guide
    """
    if modes is None:
        modes = sim.modes

    _radius = radius  # calibration radius
    arad = radians(angle)
    ta = sqrt(abs(arad) * 2.0 / pi)
    _scale = pi * ta * _radius

    # get size of last segment based on N0 points:
    N0 = 25
    t = np.linspace(0, ta, N0)
    t = t[-3: -1]
    y, x = fresnel(t)
    _ds = hypot(x[-1] - x[-2], y[-1] - y[-2])

    def cell(
        width1=width1,
        width2=width2,
        radius=None,
        angle=angle,
        xs=xs,
        layer=layer,
        name=name,
    ):
        """Create an Euler bend element.

        Args:
        width1 (float): begin width
        width2 (float): end width
        radius (float): optional to overrule the end-radius calibraton value.
        angle (float): end angle
        xs (str): xsection name
        layer (str | tuple | int): mask layer
        name (str): optional new element name

        Returns:
            Cell: Euler bend element
        """
        with nd.Cell(name=name, cnt=True) as C:
            def CM(wl, pol):
                """Optical path length model."""
                if not C.flagCM:
                    nd.main_logger(
                        f"No compact model found for Tp_euler in cell '{C.cell_name}'. Returning 0.0.",
                        "error",
                    )
                    C.flagCM = True
                return 0.0

            C.flagCM = False
            C.instantiate = cfg.instantiate_mask_element

            # set the "scale" to go from Fresnel to waveguide.
            if radius is not None:
                scale = _scale * radius / _radius
            else:
                scale = _scale
                radius = _radius
    
            ta = sqrt(abs(radians(angle)) * 2.0 / pi)
            tradius = scale / (pi * ta)
            if scale > 0:
                for lay, grow, acc, polyline in layeriter(xs, layer):
                    (a1, b1), (a2, b2), c1, c2 = grow
                    ds_res = sqrt(8 * acc / radius) 
                    ratio = _ds / ds_res  
                    N = max(2, int(N0 * ratio))
                    #print(f"{N:4}, {acc:6.3}, {scale}, {res}, {lay}")
                    t = np.linspace(0, ta, N)
                    y, x = fresnel(t)
                    spine = list(zip(scale * x, scale * y * np.sign(angle)))
                    shape = nd.util.polyline2edge(
                        xy=spine,
                        width1=width1,
                        width2=width2,
                        grow=grow,
                        anglei=0,
                        angleo=angle,
                    )
                    nd.Polygon(points=shape, layer=lay).put(0)
                x1, y1 = spine[-1]
            else:
                x1, y1 = 0, 0

            C.length_geo = 2 * tradius * abs(arad)
            p1 = nd.Pin('a0', io=0, width=width1, xs=xs, radius=0, show=True).put(0, 0, 180)
            p2 = nd.Pin('b0', io=1, width=width2, xs=xs, radius=radius, show=True).put(x1, y1, angle)

            nd.connect_path(p1, p2, C.length_geo)
            for i in modes:
                nd.connect_path(p1, p2, CM, sigtype=f'm{i}')
        return C
    return cell


def Tp_viper(
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
    params={},
    anglei=None,
    angleo=None,
    **kwargs,
):
    """Template for a specific Viper implementation.

    Free parameter in function x, y, and w are provides via ** kwarg. Note that
    when using Tp_Viper all free parameters used in x, y, and w methodsmust be provided.

    The input and output facet of the polygon are calculated from the discretised points.
    This may differt from the actual angle at the facets.
    The angles can be set explicitly by providing angles via "anglei" and "angleo".

    Args:
        x (function): function in at least t, t in [0,1]
        y (function): function in at least t, t in [0,1]
        w (function): function in at least t, t in [0,1]
        width1 (float): begin width
        width2 (float): end width
        xs (str): xsection na,e
        layer (str | tuple | int): mask layer
        name (str): element name (default='viper')
        N (int): number of polygon points (default=200)
        epsilon (float): infinitesimal step size to calucalte the begin and end ange
        anglei (float): set explicit angle of the input facet
        angleo (float): set explicit angle of the output facet       
        **kwarg: free parameters

    Returns:
        function: Function returning a Cell object with the Viper
    """
    if params != {}:
        kwargs = params

    reserved = ['width1', 'width2', 'xs', 'layer', 'N']
    for r in reserved:
        if r in params.keys():
            raise Exception(f"keyword '{r}' is reserved. Use another name.")
    kwargs0 = kwargs
    _width1 = width1
    _width2 = width2


    def viper_cell(
        width1=_width1,
        width2=_width2,
        xs=xs,
        layer=layer,
        N=N,    
        anglei=anglei,
        angleo=angleo,
        **kwargs,
    ):
        """Specific Viper implementation.

        Args:
        width1 (float): begin width
        width2 (float): end width
        xs (str): xsection na,e
        layer (str | tuple | int): mask layer
        anglei (float): set explicit angle of the input facet
        angleo (float): set explicit angle of the output facet
        N (int): number of polygon points (default=200)        
        **kwarg: free parameters

        Returns:
            Cell: Viper element
        """
        N = N # discretization steps should be "large enough" for mask resolution
        name = 'viper_bend'
        if width1 is None:
            width1 = nd.get_xsection(xs).width
        if width2 is None:
            width2 = width1

        kw = kwargs0.copy()
        for a, b in kw.items():
            kw[a] = kwargs.get(a, b)
        # Fill in all x, y, w function parameters except t:
        X = partial(x, **kw)
        Y = partial(y, **kw)
        if width1 is None and width2 is None:
            W = partial(w, **kw)
        else:
            W = partial(w, width1=width1, width2=width2, **kw)

        xa, ya = X(0), Y(0)
        xb, yb = X(1), Y(1)
        d = epsilon
        
        if anglei is None:
            aa = degrees(atan2( Y(0)-Y(d), X(0)-X(d)))
        else:
            aa = anglei    
        
        if angleo is None:    
            ab = degrees(atan2( Y(1)-Y(1-d), X(1)-X(1-d)))
        else:
            ab = angleo

        with nd.Cell(name=name, cnt=True) as C:
            for lay, grow, acc, line in nd.layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                xy = []
                width = []
                for i in range(N):
                    t = i / (N - 1)
                    xy.append((X(t), Y(t)))
                    width.append(W(t))
                points = nd.util.polyline2edge(
                    xy,
                    width1=width,
                    anglei=180+aa,
                    angleo=ab,
                    grow=grow)
                nd.Polygon(points=points, layer=lay).put(0)
            # TODO: radius in a0 and b0 to be implemented
            nd.Pin(name='a0', io=0, width=width1, xs=xs, show=True).put(xa, ya, aa)
            nd.Pin(name='b0', io=1, width=width2, xs=xs, show=True).put(xb, yb, ab)
        return C
    return viper_cell


# ==============================================================================
# create elements
# ==============================================================================
strt = Tp_straight()
bend = Tp_arc()
ptaper = Tp_ptaper()
taper = Tp_taper()
cobra = Tp_cobra()
sinebend = Tp_sinecurve()
euler = Tp_euler()


# ==============================================================================
# create a default cell and set cp
# ==============================================================================
cfg.defaultcell = Cell(name=cfg.defaultcellname)
cfg.topcell = cfg.defaultcell
cfg.cp = cfg.defaultcell.org
