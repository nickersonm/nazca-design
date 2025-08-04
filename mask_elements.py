# -----------------------------------------------------------------------
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

# @author: Ronald Broeke and Xaveer Leijtens (c) 2016-2021
# @email: ronald.broeke@brightphotonics.eu
# -----------------------------------------------------------------------

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
from __future__ import annotations
from typing import Callable
from math import sin, cos, atan, atan2, degrees, radians, pi, sqrt, hypot, tan
from functools import partial
import numpy as np
from scipy.special import fresnel  # TODO: make scipy optional?
from nazca.netlist import Cell
from nazca import cfg
from nazca.logging import logger
from nazca.gds_base import gds_db_unit as gridsize
import nazca as nd
from nazca.util import polyline_length, polyline2edge, arc2polyline, arc2polygon
from nazca.generic_bend import curve2polyline, gb_point, gb_coefficients, sinebend_point
from nazca.simglobal import sim
from nazca.util import ProtectedPartial

# TODO: move to cfg.py
min_length = 1e-6  # minimum straight waveguide length to draw in um
min_angle = 1e-6  # minimume angle in rad to draw
RADIUS = 50.0  # default
ANGLE = 90.0  # default
WIDTH = 1.0  # default
LENGTH = 10  # default

def layeriter(xs=None, layer=None):
    """Generator yielding all layers in a xsection.

    Args:
        xs (str): xsection name
        layer (int): layer number

    Yields:
        layer, growx, growy, accuracy: iterate over all layers in <xs> and <layer>
    """
    if not cfg.generate_shapes:
        return None

    if xs is None and layer is None:
        xs = cfg.default_xs_name
        if xs not in cfg.XSdict.keys():
            handle_missing_xs(xs)

    if layer is not None:
        grow = ((0.5, 0), (-0.5, 0), 0, 0)  # a layer not as part of a xs has no grow
        layer_name = nd.get_layer(layer)  # cfg.layer_names[layer][0:2]
        lineitem = cfg.layer_table.loc[layer_name]
        try:
            accuracy = lineitem["accuracy"]
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
            if not cfg.allow_empty_xsections:
                msg = (
                    "xsection '{0}' has no layers. "
                    "Continuing by adding fallback layer '{1}' to '{0}'.\n"
                    "Recommended solutions: Use a different xsection "
                    "or add layers to this xsection:\n"
                    "add_layer2xsection('{0}', layer=<num>).".format(
                        xs, cfg.default_layers["dump"]
                    )
                )
                if cfg.redirect_unknown_layers:
                    nd.main_logger(msg, "error")
                else:
                    nd.main_logger(msg, "warning")
                    nd.add_layer2xsection(xsection=xs, layer=cfg.default_layers["dump"])
                yield cfg.default_layers["dump"], (
                    (0.5, 0),
                    (-0.5, 0),
                    0,
                    0,
                ), 0.1, False
        else:
            for layer_name, A in ML.iterrows():
                grow = (
                    (A.leftedgefactor, A.leftedgeoffset),
                    (A.rightedgefactor, A.rightedgeoffset),
                    A.growy1,
                    A.growy2,
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
        if xs not in cfg.default_xs_list:
            nd.main_logger(
                f"No xsection named '{xs}'. "
                f"Correct the name or add the xsection with: "
                f"add_xsection('{xs}') to get rid of this warning. "
                f"Already available xsections are {list(cfg.XSdict.keys())}.",
                "warning",
            )
            nd.add_xsection(xs)
        else:
            XS = nd.add_xsection(xs)
            attrs = cfg.default_xs_list[xs]
            for attr, value in attrs.items():
                setattr(XS, attr, value)

    add_layer = False
    if not hasattr(nd.get_xsection(xs), "mask_layers"):
        add_layer = True
    elif nd.get_xsection(xs).mask_layers is None:
        add_layer = True
    elif nd.get_xsection(xs).mask_layers.empty:
        add_layer = True
    if add_layer:
        if xs == cfg.default_xserror_name:
            layer = cfg.default_xs_list[cfg.default_xserror_name]["layer"]
        else:
            layer = "dump"
        layer = nd.get_layer(layer)
        nd.add_layer2xsection(xs, layer=layer)
    return layer


# TODO: define a namespace instead of the using decorators?:
class CompactModel:
    """Group compact models together."""

    @staticmethod
    def check_xs_index(xs):
        """Check if an xs object has an index, otherwise return None.

        Args:

        Returns:
        """
        if xs is None:
            return None
        return getattr(nd.get_xsection(xs), "index", None)

    @staticmethod
    def cm_strt(*, index, width, length, wl=None, pol=0, mode=0):
        """Optical path length model for a straight waveguide.

        Args:

        Returns:
        """
        if wl is None:
            wl = sim.wl
        return index.Neff(width=width, wl=wl, pol=pol, mode=mode) * length

    @staticmethod
    def cm_taper(*, index, width1, width2, length, wl=None, pol=0, mode=0):
        """Optical path length model for a linear taper.

        Args:

        Returns:
        """
        if wl is None:
            wl = sim.wl
        # TODO: this approximate
        Neff = index.Neff
        n1 = Neff(width=width1, wl=wl, pol=pol, mode=mode)
        n2 = Neff(width=width2, wl=wl, pol=pol, mode=mode)
        return 0.5 * (n1 + n2) * length

    @staticmethod
    def cm_arc(*, index, width, radius, ang, wl=None, pol=0, mode=0):
        """Optical path length model for an arc bend.

        Args:

        Returns:
        """
        if wl is None:
            wl = sim.wl
        Neff = index.Neff(width=width, radius=radius, wl=wl, pol=pol, mode=mode)
        return Neff * abs(radius * ang)


    @staticmethod
    def cm_euler(*, index, width1, width2, length, wl=None, pol=0, mode=0):
        """PLACEHOLDER Optical path length model for an Euler bend.

        *THiS NOT EXACT*, just to get the "optlen" path through as test.

        Args:

        Returns:
        """
        if wl is None:
            wl = sim.wl
        Neff = index.Neff(width=width1, wl=wl, pol=pol, mode=mode)
        return Neff * abs(length)


cnt = 0  # ordinal counter for unique naming


def Tp_straight(
    length=LENGTH,
    width=WIDTH,
    xs=None,
    layer=None,
    edge1=None,
    edge2=None,
    edgepoints=50,
    name=None,
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

    Returns:
        function: Function returning a Cell object with a straight guide
    """

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
            name = "straight"
        # assert width is not None
        if width is None:
            width = 0.0

        if length < 0:
            nd.interconnect_logger(
                f"Negative straight waveguide length of {length}.", "error"
            )

        with Cell(name=name, cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(name="a0", io=0, width=width, radius=0, xs=xs, show=True).put(
                0, 0, 180
            )
            p2 = nd.Pin(name="b0", io=1, width=width, radius=0, xs=xs, show=True).put(
                length, 0, 0
            )

            nd.connect_path(p1, p2, length)
            index = CompactModel.check_xs_index(xs)
            if index is not None:
                compact_model = ProtectedPartial(
                    CompactModel.cm_strt, index=index, length=length, width=width
                )
                nd.connect_path(p1, p2, compact_model, sigtype="optlen")

            C.updk = {
                "call": "strt",
                "parameters": {
                    "length": {"value": length, "unit": "um", "type": "float"},
                    "width": {"value": width, "unit": "um", "type": "float"},
                },
            }

            # to be removed:
            C.length_geo = length  # TODO: used by old trace method: to be removed
            C.properties["parameters"] = {
                "length": {"value": length, "unit": "um", "type": "float"},
                "width": {"value": width, "unit": "um", "type": "float"},
                "xs": xs,
                "call": "strt",
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
                                (0 - c1, width * a2 + b2),
                            ]
                        else:
                            outline = [
                                (0 - c1 - 2 * gridpatch, width * a1 + b1 - gridpatch),
                                (0 - c1, width * a1 + b1),
                                (length + c2, width * a1 + b1),
                                (
                                    length + c2 + 2 * gridpatch,
                                    width * a1 + b1 - gridpatch,
                                ),
                                (
                                    length + c2 + 2 * gridpatch,
                                    width * a2 + b2 + gridpatch,
                                ),
                                (length + c2, width * a2 + b2),
                                (0 - c1, width * a2 + b2),
                                (0 - c1 - 2 * gridpatch, width * a2 + b2 + gridpatch),
                            ]
                    else:
                        if edge2 is None:
                            edge2 = edge1
                            sign = -1.0
                        Fp1 = []
                        Fp2 = []
                        for t in np.linspace(0, 1, edgepoints):
                            Fp1.append((length * t, width * a1 + b1 + edge1(t)))
                            Fp2.append((length * t, width * a2 + b2 + sign * edge2(t)))
                        outline = Fp1 + list(reversed(Fp2))
                    if abs(length + c1 + c2) > min_length:
                        nd.Polygon(layer=lay, points=outline).put(0)
                else:  # polyline
                    if abs(length + c1 + c2) > min_length:
                        wpoly = width * (a1 - a2) + (b1 - b2)
                        centre = 0.5 * (width * (a1 + a2) + (b1 + b2))
                        centreline = [(0 - c1, centre), (length + c2, centre)]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C

    return cell


def Tp_taper(
    length=LENGTH,
    width1=WIDTH,
    width2=2*WIDTH,
    shift=0,
    xs=None,
    layer=None,
    name=None,
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
            name = "taper"

        with Cell(name=name, cnt=True) as C:
            # l, w for geom/optical length/width
            l, w = length, (width1 + width2) / 2
            if shift != 0:
                l = sqrt(length ** 2 + shift ** 2)
                # projected average width
                w = w * cos(atan(shift / length))

            C.instantiate = cfg.instantiate_mask_element
            if abs(length) < gridsize / 2:
                return C  # Empty taper
            if width1 > width2:
                swap = True
                pin = ["b0", "a0"]
                width1, width2 = width2, width1
            else:
                swap = False
                pin = ["a0", "b0"]

            p1 = nd.Pin(
                name=pin[0], io=0, width=width1, radius=0, xs=xs, show=True
            ).put(0, 0, 180)
            p2 = nd.Pin(
                name=pin[1], io=1, width=width2, radius=0, xs=xs, show=True
            ).put(length, shift, 0)

            nd.connect_path(p1, p2, length)
            index = CompactModel.check_xs_index(xs)
            if index is not None:
                compact_model = ProtectedPartial(
                    CompactModel.cm_taper,
                    index=index,
                    length=length,
                    width1=width1,
                    width2=width2,
                )
                nd.connect_path(p1, p2, compact_model, sigtype="optlen")

            C.updk = {
                "call": "taper",
                "parameters": {
                    "length": {"value": length, "unit": "um", "type": "float"},
                    "width1": {"value": width1, "unit": "um", "type": "float"},
                    "width2": {"value": width2, "unit": "um", "type": "float"},
                    "shift": {"value": shift, "unit": "um", "type": "float"},
                },
            }

            # to remove:
            C.length_geo = sqrt(length ** 2.0 + shift ** 2.0)
            C.properties["parameters"] = {
                "length": {"value": length, "unit": "um", "type": "float"},
                "width1": {"value": width1, "unit": "um", "type": "float"},
                "width2": {"value": width2, "unit": "um", "type": "float"},
                "shift": {"value": shift, "unit": "um", "type": "float"},
                "xs": xs,
                "call": "taper",
            }

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if swap:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # Set widths for this layer
                twidth1a = a1 * width1 + b1
                twidth1b = a2 * width1 + b2
                twidth2a = a1 * width2 + b1 + shift
                twidth2b = a2 * width2 + b2 + shift
                if not polyline:
                    outline = [
                        (0 - c1, twidth1a),
                        (length + c2, twidth2a),
                        (length + c2, twidth2b),
                        (0 - c1, twidth1b),
                    ]
                    nd.Polygon(layer=lay, points=outline).put()
                else:
                    if abs(length) > min_length:
                        wpoly = twidth1a - twidth1b
                        centre = 0.5 * (twidth1a + twidth1b)
                        centreline = [(0 - c1, centre), (length + c2, centre + shift)]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C

    return cell


def Tp_ptaper(
    length=LENGTH,
    width1=WIDTH,
    width2=2*WIDTH,
    xs=None,
    layer=None,
    name=None,
):
    """Template for creating a parametrized parabolic taper  function.

    Note that zero length taper segments may seem to overlook width pin2pin DRC.

    Args:
        length (float): length of the taper
        width1 (float): width at start
        width2 (float): width at end
        xs (str): xsection of taper
        layer (int | str): layer number or layer name
        name (str): Name of the cell. Default is None.

    Returns:
        function: Function returning a Cell object with a ptaper
    """

    def cell(
        length=length,
        width1=width1,
        width2=width2,
        xs=xs,
        layer=layer,
        name=name,
        **kwargs
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
            name = "ptaper"

        with Cell(name=name, cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            if abs(length) < gridsize / 2.0:
                return C  # Empty taper
            if width1 > width2:
                swap = True
                pin = ["b0", "a0"]
                width1, width2 = width2, width1
            else:
                swap = False
                pin = ["a0", "b0"]

            p1 = nd.Pin(
                name=pin[0], io=0, width=width1, radius=0, xs=xs, show=True
            ).put(0, 0, 180)
            p2 = nd.Pin(
                name=pin[1], io=1, width=width2, radius=0, xs=xs, show=True
            ).put(length, 0, 0)

            nd.connect_path(p1, p2, length)
            index = CompactModel.check_xs_index(xs)
            if index is not None:
                compact_model = ProtectedPartial(
                    CompactModel.cm_taper,
                    index=index,
                    length=length,
                    width1=width1,
                    width2=width2,
                )
                nd.connect_path(p1, p2, compact_model, sigtype="optlen")

            C.updk = {
                "call": "taper",
                "parameters": {
                    "length": {"value": length, "unit": "um", "type": "float"},
                    "width1": {"value": width1, "unit": "um", "type": "float"},
                    "width2": {"value": width2, "unit": "um", "type": "float"},
                },
            }

            # to remove:
            C.length_geo = length
            C.properties["parameters"] = {
                "length": {"value": length, "unit": "um", "type": "float"},
                "width1": {"value": width1, "unit": "um", "type": "float"},
                "width2": {"value": width2, "unit": "um", "type": "float"},
                "xs": xs,
                "call": "ptaper",
            }

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if swap:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                if not polyline:
                    # Set widths for this layer
                    twidth1a = a1 * width1 + b1
                    twidth1b = a2 * width1 + b2
                    twidth2a = a1 * width2 + b1
                    twidth2b = a2 * width2 + b2
                    if abs(0.5 * width1 - 0.5 * width2) / (4 * length) < gridsize:
                        # No shaping required due to small delta width:
                        points = [
                            (0, twidth1a),
                            (length, twidth2a),
                            (length, twidth2b),
                            (0, twidth1b),
                        ]
                        nd.Polygon(layer=lay, points=points).put(0)
                        continue
                    # y = a * x**2
                    a = 4 * length / (width2 ** 2 - width1 ** 2)
                    y2 = y0 = a * (0.5 * width1) ** 2
                    w_tap1 = width1
                    ptop = [(0, twidth1a)]
                    pbot = [(0, twidth1b)]
                    while y2 - y0 < length:
                        w_tap2 = w_tap1 + 4 * acc + 4 * sqrt(acc * (w_tap1 + acc))
                        y2 = a * (0.5 * w_tap2) ** 2
                        if y2 - y0 < length:
                            ptop.append((y2 - y0, w_tap2 * a1 + b1))
                            pbot.append((y2 - y0, w_tap2 * a2 + b2))
                            w_tap1 = w_tap2
                        else:
                            break
                    ptop.append((length, twidth2a))
                    pbot.append((length, twidth2b))
                    nd.Polygon(layer=lay, points=ptop + list(reversed(pbot))).put(0)
                else:
                    if abs(length) > min_length:
                        wpoly = width1 * (a1 - a2) + (b1 - b2)
                        centre = 0.5 * (width1 * (a1 + a2) + (b1 + b2))
                        centreline = [(0 - c1, centre), (length + c2, centre)]
                        nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)
        return C

    return cell


def __get_offset(xs, width, radius, offset=None):
    """Get the offset function.

    Args:
        xs (str): xsection of the arc bend
        width (float): Width of the arc bend
        radius (float): Bending radius of the arc
        offset (float): Overriding straight-bend offset

    Returns:
        Any: The straight-bend offset
    """
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
    radius: float = RADIUS,
    width: float = WIDTH,
    width2: float = None,
    angle: float = 90,
    xs: str = None,
    layer: int | str = None,
    offset: float | Callable = None,
    offset2: float | Callable = None,
    parabolic: bool = True,
    name: str = None,
):
    """Template for creating a parameterized circular arc waveguide function.

    Args:
        radius (float): radius at the center line of the arc in um.
        width (float): width of the arc in um. An arbitrary profile can be set by passing a function.
        width2 (float): Second width of the arc in um. If not None, the bend will be tapered
        from width1 to width2 with a profile determined by the 'parabolic' argument. Default is None.
        angle (float): angle of arc in degree (default = 90).
        xs (str): xsection of taper
        layer (int | str): layer number or layer name
        offset (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        offset2 (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        parabolic (bool): Makes the taper parabolic if true and linear otherwise. Default is True.
        name (str): Name of the cell. Default is None, leading to "arc".

    Returns:
        function: Function returning a Cell object with an arc
    """

    # TODO: add growy
    def cell(
        radius: float = radius,
        width: float = width,
        width2: float = width2,
        angle: float = angle,
        xs: str = xs,
        layer: int | str = layer,
        offset: float | Callable = offset,
        offset2: float | Callable = offset2,
        parabolic: bool = parabolic,
        name: str = name,
    ) -> Cell:
        """Create a circular arc element.

        A straight-bend offset is included when it has been defined in the
        xsection used.

        Args:
            radius (float): radius at the center line of the arc in um.
            width (float): width of the arc in um. An arbitrary profile can be set by passing a function.
            width2 (float): Second width of the arc in um. If not None, the bend will be tapered
            from width1 to width2 with a profile determined by the 'parabolic' argument. Default is None.
            angle (float): angle of arc in degree (default = 90).
            xs (str): xsection of taper
            layer (int | str): layer number or layer name
            offset (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float
            offset2 (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float
            parabolic (bool): Makes the taper parabolic if true and linear otherwise. Default is True.
            name (str): Name of the cell. Default is None, leading to "arc".

        Returns:
            Cell: circular arc element
        """
        if name is None:
            name = "arc"
        if width2 is None:
            width2 = width
        ang = np.radians(angle)
        if abs(ang) < min_angle:
            ang = 0
            angle = 0
        sign = np.sign(ang)
        radius = abs(radius)

        offset = __get_offset(xs, width, radius, offset=offset)
        offset2 = __get_offset(xs, width2, radius, offset=offset2)

        Nmax = cfg.maxpolygonpoints // 2 - 2  # 2xnmax + 4

        with Cell(name=name, cnt=True) as C:
            # C.use_hull = True
            C.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(
                name="a0", io=0, width=width, radius=radius, xs=xs, show=True
            ).put(0, 0, 180)
            p2 = nd.Pin(
                name="b0", io=1, width=width2, radius=radius, xs=xs, show=True
            ).put(radius * sin(abs(ang)), sign * radius * (1 - cos(ang)), angle)

            nd.connect_path(p1, p2, abs(radius * ang))
            index = CompactModel.check_xs_index(xs)
            if index is not None:
                compact_model = ProtectedPartial(
                    CompactModel.cm_arc,
                    index=index,
                    ang=ang,
                    radius=radius,
                    width=(width + width2) / 2,  # TODO: not strictly correct for tapered bends
                )
                nd.connect_path(p1, p2, compact_model, sigtype="optlen")

            C.updk = {
                "call": "bend",
                "parameters": {
                    "angle": {"value": angle, "unit": "deg", "type": "float"},
                    "width": {"value": width, "unit": "um", "type": "float"},
                    "width2": {"value": width2, "unit": "um", "type": "float"},
                    "radius": {"value": radius, "unit": "um", "type": "float"},
                    "parabolic": {"value": parabolic, "unit": None, "type": "bool"},
                },
            }

            # to be removed:
            C.length_geo = abs(radius * ang)
            C.properties["parameters"] = {
                "angle": {"value": angle, "unit": "deg", "type": "float"},
                "width": {"value": width, "unit": "um", "type": "float"},
                "width2": {"value": width2, "unit": "um", "type": "float"},
                "radius": {"value": radius, "unit": "um", "type": "float"},
                "parabolic": {"value": parabolic, "unit": None, "type": "bool"},
                "xs": xs,
                "call": "bend",
            }

            if ang == 0:
                return C

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if angle >= 0:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # Start and end sections are different from other sections
                wi = width * (a1 - a2) + (b1 - b2)
                wo = width2 * (a1 - a2) + (b1 - b2)

                r = radius + 0.5 * ((width + width2) / 2 * (a1 + a2) + (b1 + b2)) - offset
                # TODO: add growy
                Ri = r - wi / 2
                Ro = r + wo / 2
                # Inner radius smaller than zero is taken care of in arc2polygon.
                if Ri < 0:
                    if Ri < -1e-6:
                        nd.main_logger(
                            "Inner side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f}, xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".format(
                                cfg.cells[-1].cell_name, Ri, radius, xs, width, offset
                            ),
                            "warning",
                        )
                if Ro < 0:
                    if Ro < -1e-6:
                        nd.main_logger(
                            "Inner side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f}, xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".format(
                                cfg.cells[-1].cell_name, Ro, radius, xs, width, offset
                            ),
                            "warning",
                        )
                if not isinstance(acc, float):
                    acc = 0.001
                    logger.error(
                        f"Accuracy not defined for layer {lay}. Setting it to {acc}."
                    )
                if not polyline:
                    outline = arc2polygon(radius=r, angle=angle, width1=wi, width2=wo, accuracy=acc, parabolic=parabolic)
                    nd.Polygon(layer=lay, points=outline).put(0, sign * (radius - r))
                else:
                    midline = arc2polyline(radius=r, angle=angle, accuracy=acc)
                    nd.Polyline(layer=lay, points=midline, width=(wi + wo) / 2).put(
                        0, sign * (radius - r)
                    )
        return C

    return cell


def Tp_arc2(
    radius: float = RADIUS,
    width: float = WIDTH,
    width2: float = None,
    angle: float = 90,
    xs: str = None,
    layer: int | str = None,
    offset: float | Callable = None,
    offset2: float | Callable = None,
    parabolic: bool = False,
    name: str = None,
) -> Callable:
    """Template for creating a parametrized angled arc waveguide function.

    The arc is composed of straight sections.
    Minimum angle between sections is 45 deg.

    Args:
        radius (float): radius at the center line of the arc in um.
        width (float): width of the arc in um. An arbitrary profile can be set by passing a function.
        width2 (float): Second width of the arc in um. If not None, the bend will be tapered
        from width1 to width2 with a profile determined by the 'parabolic' argument. Default is None.
        angle (float): angle of arc in degree (default = 90).
        xs (str): xsection of taper
        layer (int | str): layer number or layer name
        offset (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        offset2 (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float
        parabolic (bool): Makes the taper parabolic if true and linear otherwise. Default is True.
        name (str): Name of the cell. Default is None, leading to "arc".

    Returns:
        function: Function returning a Cell object with an arc composed of straight sections
    """

    # TODO: add growy
    def cell(
        radius: float = radius,
        width: float = width,
        width2: float = width2,
        angle: float = angle,
        xs: str = xs,
        layer: int | str = layer,
        offset: float | Callable = offset,
        offset2: float | Callable = offset2,
        parabolic: bool = parabolic,
        name: str = name,
    ) -> Cell:
        """Create a single circular arc element composed of straight sections.

        A straight-bend offset is included when it has been defined in the
        xsection used.

        Args:
            radius (float): radius at the center line of the arc in um.
            width (float): width of the arc in um. An arbitrary profile can be set by passing a function.
            width2 (float): Second width of the arc in um. If not None, the bend will be tapered
            from width1 to width2 with a profile determined by the 'parabolic' argument. Default is None.
            angle (float): angle of arc in degree (default = 90).
            xs (str): xsection of taper
            layer (int | str): layer number or layer name
            offset (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float
            offset2 (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float
            parabolic (bool): Makes the taper parabolic if true and linear otherwise. Default is True.
            name (str): Name of the cell. Default is None, leading to "arc".

        Returns:
            Cell: circular arc element
        """
        if name is None:
            name = "arc"
        ang = np.radians(angle)
        if abs(ang) < min_angle:
            ang = 0
            angle = 0
        sign = np.sign(ang)
        radius = abs(radius)

        # if offset is None:
        offset = __get_offset(xs, width, radius, offset=offset)

        Nmax = cfg.maxpolygonpoints // 2 - 2  # 2xnmax + 4
        da = pi / 4.0
        Ni = int((abs(ang) / da) - 1e-6) + 1
        ang = ang / (1.0 * Ni)  # step angle for N
        with Cell(name="int", cnt=True) as INT:
            # C.use_hull = True
            INT.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(name="a0", io=0, width=width, xs=xs, show=True).put(0, 0, 180)
            p2 = nd.Pin(name="b0", io=1, width=width, xs=xs, show=True).put(
                radius * sin(abs(ang)), sign * radius * (1 - cos(ang)), 180.0 / pi * ang
            )

            nd.connect_path(p1, p2, 2.0 * abs(radius * tan(ang * 0.5)))
            INT.length_geo = abs(radius * ang)

            if ang == 0:
                return INT

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                if angle >= 0:
                    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # Start and end sections are different from other sections
                if polyline:
                    wpoly = width * (a1 - a2) + (b1 - b2)
                    Ro = Ri = radius + 0.5 * (width * (a1 + a2) + (b1 + b2)) - offset
                else:
                    Ro = radius + width * a1 - offset + b1  # outer radius
                    Ri = radius + width * a2 - offset + b2  # inner radius
                    # TODO: add growy
                if Ri < 0:
                    if Ri < -1e-6:
                        nd.main_logger(
                            "Side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f},"
                            " xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".format(
                                cfg.cells[-1].cell_name, Ri, radius, xs, width, offset
                            ),
                            "warning",
                        )
                    Ri = 0
                if Ro < 0:
                    if Ro < -1e-6:
                        nd.main_logger(
                            "Side wall arc radius too small in cell '{}': side_radius={:.3f}, radius={:.3f},"
                            " xs='{}', width={:.3f}, offset={:.3f}. Setting it to 0.".format(
                                cfg.cells[-1].cell_name, Ro, radius, xs, width, offset
                            ),
                            "warning",
                        )
                    Ro = 0
                Rmax = max(Ri, Ro)
                Roe = Ro / cos(abs(ang) / 2.0)  # effective outer radius
                Rie = Ri / cos(abs(ang) / 2.0)  # effective inner radius
                p1 = [
                    (Roe * sin(abs(0.5 * ang)), sign * (radius - Roe * cos(0.5 * ang)))
                ]
                p2 = [
                    (Rie * sin(abs(0.5 * ang)), sign * (radius - Rie * cos(0.5 * ang)))
                ]
                pstart = [(0, sign * (radius - Ri)), (0, sign * (radius - Ro))]
                pend = [
                    (Ro * sin(abs(ang)), sign * (radius - Ro * cos(ang))),
                    (Ri * sin(abs(ang)), sign * (radius - Ri * cos(ang))),
                ]
                if not polyline:
                    outline = pstart + p1 + pend + list(reversed(p2))
                    nd.Polygon(layer=lay, points=outline).put(0)
                else:
                    centreline = [pstart[0]] + p1 + [pend[0]]
                    nd.Polyline(layer=lay, points=centreline, width=wpoly).put(0)

        inst = []
        with Cell(name="arc", cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            for i in range(Ni):
                inst.append(INT.put())
            nd.Pin("a0", pin=inst[0].pin["a0"]).put()
            nd.Pin("b0", pin=inst[-1].pin["b0"]).put()

        C.updk = {
            "parameters": {
                "angle": {"value": angle, "unit": "deg", "type": "float"},
                "width": {"value": width, "unit": "um", "type": "float"},
                "width2": {"value": width2, "unit": "um", "type": "float"},
                "radius": {"value": radius, "unit": "um", "type": "float"},
                "parabolic": {"value": parabolic, "unit": None, "type": "bool"},
            },
            "call": "bend",
        }

        C.properties["parameters"] = {
            "angle": {"value": angle, "unit": "deg", "type": "float"},
            "width": {"value": width, "unit": "um", "type": "float"},
            "width2": {"value": width2, "unit": "um", "type": "float"},
            "radius": {"value": radius, "unit": "um", "type": "float"},
            "parabolic": {"value": parabolic, "unit": None, "type": "bool"},
            "xs": {"value": xs, "unit": None, "type": "str"},
            "call": {"value": "bend", "unit": None, "type": "str"},
        }
        return C

    return cell


def Tp_sinecurve(
    width=WIDTH,
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
        width=width, distance=distance, offset=offset, xs=xs, layer=layer, name=name
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
            name = "sinecurve"

        xya = (distance, offset, 0)  # End point
        # sinecurve waveguide, cwg
        with Cell(name=name, cnt=True) as C:
            C.instantiate = cfg.instantiate_mask_element
            nd.Pin(name="a0", io=0, width=width, radius=0, xs=xs, show=True).put(
                0, 0, 180
            )
            nd.Pin(name="b0", io=1, width=width, radius=0, xs=xs, show=True).put(*xya)
            # TODO: length_geo
            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)
                # TODO: define polyline2polygon with left and right width from spine
                # sampled curve
                xy = curve2polyline(sinebend_point, xya, acc, (distance, offset))
                if not polyline:
                    # polygon of proper width
                    xy = polyline2edge(xy, width, grow=grow)
                    nd.Polygon(layer=lay, points=xy).put(0)
                else:
                    wpoly = width * (a1 - a2) + (b1 - b2)
                    xy = polyline2edge(xy, width, grow=grow, line=True)
                    nd.Polyline(layer=lay, points=xy, width=abs(wpoly)).put(0)

            C.updk = {
                "call": "sinecurve",
                "parameters": {
                    "distance": {"value": distance, "unit": "um", "type": "float"},
                    "width": {"value": width, "unit": "um", "type": "float"},
                    "offset": {"value": offset, "unit": "um", "type": "float"},
                },
            }

            # to remove
            C.properties["parameters"] = {
                "distance": {"value": distance, "unit": "um", "type": "float"},
                "width": {"value": width, "unit": "um", "type": "float"},
                "offset": {"value": offset, "unit": "um", "type": "float"},
                "xs": {"value": xs},
                "call": {"value": "sinecurve"},
            }
        return C

    return cell


def Tp_cobra(
    xya=(100, 100, 10),
    width1=WIDTH,
    width2=WIDTH,
    radius1=0,
    radius2=0,
    offset1=None,
    offset2=None,
    parabolic=True,
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
        parabolic=True,
        xs=xs,
        layer=layer,
        name=name,
        shift=0,
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
            name = "cobra"
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

            C.flagCM = False
            C.instantiate = cfg.instantiate_mask_element
            p1 = nd.Pin(
                name="a0", io=0, width=width1, xs=xs, radius=radius1, show=True
            ).put(0, -offset1 + shift, 180)

            p2 = nd.Pin(
                name="b0", io=1, width=width2, xs=xs, radius=radius2, show=True
            ).put(
                xya[0] - shift * sin(np.radians(xya[2])),
                xya[1] - offset1 + shift * cos(np.radians(xya[2])),
                xya[2],
            )

            xya = (
                xya[0],  # - offset2 * sin(np.radians(xya[2])),
                xya[1],  # - offset1 + offset2 * cos(np.radians(xya[2])),
                xya[2],
            )
            # Solve the generic bend

            A, B, L, Rmin = gb_coefficients(xya, radius1=radius1, radius2=radius2)

            C.length_geo = None
            C.properties["Rmin"] = Rmin
            saved_acc = 10

            Rdrc = 1e8  # initial very large min DRC radius
            try:
                Rdrc = cfg.XSdict[xs].minimum_radius
                if Rdrc is None:
                    Rdrc = 0
            except KeyError:
                if xs not in cfg.XSdict.keys():
                    Rdrc = 0
            if Rdrc > Rmin:
                nd.main_logger(
                    f"DRC minimum_radius in Cobra {Rmin:.3f} < {Rdrc:.3f}", "warning"
                )

            for lay, grow, acc, polyline in layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                # sampled curve
                spine = curve2polyline(gb_point, xya, acc, (A, B, L))
                if C.length_geo is None or saved_acc > acc:
                    # Calculate length of polyline and keep the one for the
                    # most accurate layer
                    C.length_geo = polyline_length(spine)
                    saved_acc = acc
                # polygon of proper width
                # TODO: add growy
                if not polyline:
                    if callable(width1):

                        def wfunc(t):
                            return width1(t)

                        edge = polyline2edge(
                            spine,
                            wfunc,
                            grow=grow,
                            shift=-shift,
                        )
                        if polyline:
                            nd.Polyline(points=edge, layer=lay, width=wfunc(0)).put(0)
                        else:
                            nd.Polygon(layer=lay, points=edge).put(0)
                    else:
                        edge = polyline2edge(
                            spine,
                            width1,
                            width2,
                            parabolic=parabolic,
                            grow=grow,
                            shift=-shift,
                        )
                        if polyline:
                            nd.Polyline(points=edge, layer=lay, width=width1).put(0)
                        else:
                            nd.Polygon(layer=lay, points=edge).put(0)
                else:
                    wpoly = width1 * (a1 - a2) + (b1 - b2)
                    xy = polyline2edge(spine, width1, grow=grow, line=True)
                    nd.Polyline(layer=lay, points=xy, width=abs(wpoly)).put(0)

            nd.connect_path(p1, p2, C.length_geo)

        return C

    return cell


def euler_scale(radius, angle):
    """Calculate the scale of an Euler bend having <radius> at <angle>.

    Args:
        radius (float): radius at <angle>.
        angle (float): angle in degrees.

    Returns:
        float: scale of Euler
    """
    arad = radians(angle)
    ta = sqrt(abs(arad) * 2.0 / pi)  # parameter t at angle a
    scale = pi * ta * radius
    return scale


def euler_radius(scale, angle):
    """Calculate the radius at angle <angle> of an Euler bend having scale <scale>.

    Args:
        scale (float): Euler bend scale.
        angle (float): angle in degrees.

    Returns:
        float: scale of Euler
    """
    arad = radians(angle)
    ta = sqrt(abs(arad) * 2.0 / pi)  # parameter t at angle a
    radius = scale / (pi * ta)
    return radius


def _calibrate_euler(scale=None, radius=None, angle=None, N=25, minradius=None, xs=None):
    """Calculate the Euler properties for the layout of a specific bend for a known scale.

    Calculate an Euler for a given scale. Only the radius OR the angle should be provided;
    The scale and angle will set the radius, while the scale and radius will set the angle.
    Giving both angle and radius would change the scale, which is not allowed here.

    Args:
        radius (float): radius at <angle>.
        angle (float): angle in degrees
        N (int): number of segments in initial Euler discretization.
        minradius (float): ?
        scale (float): scale of the Euler.

    Returns:
        float, float, float: scale factor, length of last of N segments, angle at which Rmin is reached.
    """
    def get_scale(scale):
        """Get scale from xsection if it can be calculated from user explicit input."""
        if scale is None:
            scale = getattr(xs, "scale", None)
        if scale is None:
            raise Exception("Euler scale is undefined. Aborting")
        return scale

    if angle is None and radius is None:
        raise Exception(
            "Euler calibration: Need an angle or radius parameter, but None is provided. "
        )
    if angle is not None and radius is not None:
        nd.main_logger(
            "Euler calibration: got an angle and radius. Only one is allowed. Will set ignore the radius",
            "error",
        )
    if angle is not None:
        if abs(angle) < 1e-6:
            ta = 0
            radius = 0
        else:
            ta = sqrt(abs(radians(angle)) * 2.0 / pi)
            radius = get_scale(scale) / (pi * ta)  # @ta
    else:
        ta = get_scale(scale) / (pi * radius)
        angle = degrees(ta ** 2 * pi / 2.0)


    if minradius is None:
        maxangle = angle
    else:
        minradius = max(minradius, 2.5)
        tamax = scale / (pi * minradius)
        maxangle = degrees(tamax**2 * pi / 2.0)
    # print(f"{maxangle=}")

    # get length of last segment based on N points:
    N = max(3, round(N * angle / 90))  # Increase N with the angle to keep ds segments from increasing too much.
    t = np.linspace(0, ta, N)
    t = t[-3:-1]
    y, x = fresnel(t)
    ds = hypot(x[-1] - x[-2], y[-1] - y[-2])
    return radius, angle, ds, maxangle, ta


def Tp_euler(
    width: float = WIDTH,
    width2: float = WIDTH,
    radius: float = RADIUS,
    angle: float = 90,
    xs: str = None,
    layer: str | tuple | int = None,
    scale: float = None,
    name: str = "euler",
    parabolic: bool = True,
) -> Callable:
    """Template for creating an Euler bend.

    The Euler is calibrated to have <radius> at a specific <angle>.

    If later only an angle or only a radius is provided, the Euler function
    will then always trace the same shape. For angle shorter than <angle> it
    will have a radius larger than <radius> and vise versa.

    Args:
        width1 (float): begin width
        width2 (float): end width
        radius (float): end radius for calibration the derived curve (default=50).
        angle (float): end angle in degrees for calibration the derived curve (default=90).
        xs (str): xsection name
        layer (str | tuple | int): mask layer
        scale (float): euler scale
        name (str): element name (default='euler')
        parabolic (bool): uses a parabolic width profile if width1 != width2. Default is True.

    Returns:
        function: Function returning a Cell object with the Euler bend
    """
    try:
        Rdrc = cfg.XSdict[xs].minimum_radius
        if Rdrc is None:
            Rdrc = 0
    except KeyError:
        if xs not in cfg.XSdict.keys():
            Rdrc = 0

    N0 = 25
    _radius_cal = radius
    _angle_cal = angle
    if scale is not None:
        _scale = scale
    else:
        _scale = euler_scale(radius=radius, angle=angle)
   # _ds, _maxangle = _calibrate_euler(
   #     radius=_radius_cal, angle=_angle_cal, N=N0, minradius=Rdrc, scale=scale
    #)

    @nd.bb_util.hashme("euler", ["xs", "angle"])
    def cell(
        width: float = width,
        width2: float = width2,
        radius: float = None,
        angle: float = None,
        scale: float = None,
        radius_cal: float = None,
        angle_cal: float = None,
        xs: str = xs,
        layer: str | tuple | int = layer,
        name: str = name,
        solve: bool = False,
        parabolic: bool = parabolic
    ) -> Cell | tuple[float]:
        """Create an Euler bend element.

        The Euler is calibrated to reach a fixed radius at a fixed angle.
        Changing the angle will only change the length of the curve,
        not its shape. A recalibration for a specific Euler call can be accomplished

        by providing both radius_cal and angle_cal, or the Euler scale.
        Note the scale can be obtained from the function euler_scale().

        Args:
            width1 (float): begin width
            width2 (float): end width
            angle (float): end angle
            radius_cal (float): optional to overrule the end-radius calibraton
                value at the calibration angle, need angle_cal as well.
            angle_cal (float): optional to overrule the end-angle calibration,
                needs radius_cal as well.
            scale (float):  Directly give the Euler scale. Alternative to providing radius_cal and angle_cal.
            xs (str): xsection name
            layer (str | tuple | int): mask layer
            name (str): optional new element name
            solve (bool): Return solved values of (radius, angle, maxangle) if True. Default is False
            parabolic (bool): uses a parabolic width profile if width1 != width2. Default is True.
            solve (bool): Return solved values of Euler (radius, angle, maxangle) if True (do not build the Cell). Default is False

        Returns:
            Cell: Euler bend element if solve is False, otherwise solver values (float, float, float): radius, angle, maxangle
        """
        # scaling options:
        epsilon = 1e-6  # margin of error to avoid a false angle DRC.
        # avoid reassigning of closure variables:
        #ds = _ds
        #maxangle = _maxangle
        Scal = _scale
        Rcal = _radius_cal
        Acal = _angle_cal

        try:
            Rdrc = cfg.XSdict[xs].minimum_radius
            if Rdrc is None:
                Rdrc = 0
        except KeyError:
            if xs not in cfg.XSdict.keys():
                Rdrc = 0

        if scale is not None and (radius_cal is not None or angle_cal is not None):
            nd.main_logger(
                "Ambiguous Euler recalibration is not allowed. "
                "Provided radius_cal and angle_cal, *or* parameters instead "
                f"and only set the radius OR angle. Continuing with the explicitly provided scale here ({scale}).",
                "error",
            )

        if radius is not None and angle is not None:
            nd.main_logger(
                "Implicit Euler recalibration is not allowed (both radius and angle were provided) "
                "Use radius_cal and angle_cal or scale parameters for calibration instead "
                f"and only set the radius OR the angle. Continuing with the angle here ({angle}), "
                " dropping the radius value and maintaining existing scaling.",
                "error",
            )
            radius = None

        # calculate Scal if new scaling info is provided:
        if scale is not None:
            Scal = scale
        else:
            if radius_cal is not None and angle_cal is not None:
                Rcal = radius_cal
                Acal = angle_cal
                Scal = euler_scale(radius=Rcal, angle=Acal)
            #elif radius_cal is None and angle_cal is None:
            #    nd.main_logger(
            #        "Incomplete calibration set of parameters provided. Need or scale or both radius_cal and angle_cal. "
            #        "Assuming no recalibration is needed. Ignoring complete parameters."
            #        "error",
            #    )

        if angle is None and radius is None:
            angle = ANGLE
        if angle is not None:
            signangle = np.sign(angle)
        else:
            signangle = 1


        taradius, taangle, ds, maxangle, ta = _calibrate_euler(
            radius=radius, angle=angle, scale=Scal, N=N0, minradius=Rdrc
        )

        if solve:
            return taradius, taangle, maxangle

        if Rdrc - epsilon > taradius:
            nd.main_logger(
                f"DRC on minimum_radius in Euler {taradius:.3f} < {Rdrc:.3f}",
                "error",
            )

        with nd.Cell(name=name, cnt=True) as C:
            C.updk = {
                "call": "euler",
                "parameters": {
                    "angle": {"value": taangle, "unit": "deg", "type": "float"},
                    "radius": {"value": taradius, "unit": "um", "type": "float"},
                    "width": {"value": width, "unit": "um", "type": "float"},
                    "width2": {"value": width2, "unit": "um", "type": "float"},
                    "radius_cal": {"value": Rcal, "unit": "um", "type": "float"},
                    "angle_cal": {"value": Acal, "unit": "deg", "type": "float"},
                    "scale": {"value": Scal, "unit": "none", "type": "float"},
                    "xs": {"value": xs},
                },
            }

            C.flagCM = False
            C.instantiate = cfg.instantiate_mask_element

            if Scal > 0 and taradius > 0:
                for lay, grow, acc, polyline in layeriter(xs, layer):
                    (a1, b1), (a2, b2), c1, c2 = grow
                    ds_res = sqrt((90 / abs(taangle)) * acc * taradius / Rcal ** 2)
                    ratio = ds / ds_res
                    N = max(2, int(N0 * ratio))
                    # print(f"{N:4}, {acc:6.3}, {Scal}, {res}, {lay}")
                    t = np.linspace(0, ta, N)
                    y, x = fresnel(t)
                    spine = list(zip(Scal * x, Scal * y * signangle))
                    shape = nd.util.polyline2edge(
                        xy=spine,
                        width1=width,
                        width2=width2,
                        grow=grow,
                        anglei=0,
                        angleo=taangle,
                    )
                    if polyline:
                        nd.Polyline(layer=lay, points=spine, width=width).put(0)
                    else:
                        nd.Polygon(points=shape, layer=lay).put(0)

                x1, y1 = spine[-1]
            else:
                x1, y1 = 0, 0

            C.length_geo = 2 * taradius * abs(radians(taangle))
            p1 = nd.Pin("a0", io=0, width=width, xs=xs, radius=0, show=True).put(
                0, 0, 180
            )
            p2 = nd.Pin(
                "b0", io=1, width=width2, xs=xs, radius=taradius, show=True
            ).put(x1, y1, angle)

            nd.connect_path(p1, p2, C.length_geo)

            index = CompactModel.check_xs_index(xs)
            if index is not None:
                compact_model = ProtectedPartial(
                    CompactModel.cm_euler, index=index, length=C.length_geo, width1=width, width2=width2
                )
                nd.connect_path(p1, p2, compact_model, sigtype="optlen")

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
    name="viper",
    params=None,
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
        **kwargs: free parameters

    Returns:
        function: Function returning a Cell object with the Viper
    """
    # if params != {}:
    #     kwargs = params

    reserved = ["width1", "width2", "xs", "layer", "N"]

    if params is not None:
        kwargs = params
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
        N = N  # discretization steps should be "large enough" for mask resolution
        name = "viper_bend"
        if width1 is None:
            width1 = nd.get_xsection(xs).width
        if width2 is None:
            width2 = width1

        if kwargs0 is not None:
            kw = kwargs0.copy()
            for a, b in kw.items():
                kw[a] = kwargs.get(a, b)
        else:
            kw = {}

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
            aa = degrees(atan2(Y(0) - Y(d), X(0) - X(d)))
        else:
            aa = anglei

        if angleo is None:
            ab = degrees(atan2(Y(1) - Y(1 - d), X(1) - X(1 - d)))
        else:
            ab = angleo

        with nd.Cell(name=name, cnt=True) as C:
            for lay, grow, acc, polyline in nd.layeriter(xs, layer):
                (a1, b1), (a2, b2), c1, c2 = grow
                xy = []
                width = []
                for i in range(N):
                    t = i / (N - 1)
                    xy.append((X(t), Y(t)))
                    width.append(W(t))
                points = nd.util.polyline2edge(
                    xy, width1=width, anglei=180 + aa, angleo=ab, grow=grow
                )

                if polyline:
                    nd.Polyline(points=xy, layer=lay, width=width[0]).put(0)
                else:
                    nd.Polygon(points=points, layer=lay).put(0)
            # TODO: radius in a0 and b0 to be implemented
            nd.Pin(name="a0", io=0, width=width1, xs=xs, show=True).put(xa, ya, aa)
            nd.Pin(name="b0", io=1, width=width2, xs=xs, show=True).put(xb, yb, ab)
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
