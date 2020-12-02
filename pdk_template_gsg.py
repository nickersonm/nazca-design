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
# (c) 2016-2017 Ronald Broeke, Katarzyna Lawniczuk
#==============================================================================

"""Module defining black box templates for PDK implementation."""


from math import atan, degrees
import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom
import nazca.cfg as cfg


def Tp_RFlineGSG(length=100, width_sig=10, width_gnd=10, gap=10, wext_bg=None,
            pinwidth=None, name='line_gsg', groupname='', xs=None, xs_bb=None,
            xs_bg=None, xs_bg_width=None):
    """Template for RF line GSG.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length, width_sig=width_sig, width_gnd=width_gnd, gap=gap):
        """Create a GSG RF line cell.

        Args:
            length (float): length of the RF section
            width_sig (float): width of the signal line
            width_gnd (float): width of the gnd line
            gap (float): gap between the signal and ground line

        Returns:
           Cell
        """
        with nd.Cell(name=bbu._hash_name, instantiate=False) as C:
            width = width_sig
            sig = nd.strt(length=length, width=width, xs=xs['a0']).put(0)
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='electrical').put(sig.pin['a0'])
            nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='electrical').put(sig.pin['b0'])
            bbu.put_stub(['a0', 'b0'])

            width = width_gnd
            nd.strt(length=length, width=width, xs=xs['a0']).put(0, -gap-0.5*(width_sig+width_gnd))
            nd.strt(length=length, width=width, xs=xs['a0']).put(0, +gap+0.5*(width_sig+width_gnd))

            if xs_bg is not None:
                nd.strt(length=length, width=xs_bg_width, xs=xs_bg).put(0)

            cfg.cp = C.pin['b0']
        return C
    return cell


def Tp_RFbendGSG(radius=50, angle=30, width_sig=10, width_gnd=10, gap=10, wext_bg=8,
        pinwidth=None, name='bend_gsg', groupname='', xs=None, xs_bg=None, xs_bg_width=None):
    """Template for GSG RF bend.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell(radius=radius, angle=angle, width_sig=width_sig,
            width_gnd=width_gnd, gap=gap):
        """Create a GSG RF bend cell.

        Args:
            radius (float): radius of the bend
            angle (float): bend angle in deg
            width_sig (float): width of the signal line
            width_gnd (float): width of the gnd line

        Returns:
            Cell
        """
        with nd.Cell(name=bbu._hash_name, instantiate=False) as C:
            C.groupname = groupname

            width = width_sig
            Sig = nd.bend(radius=radius, angle=angle, width=width, xs=xs['a0']).put(0)
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='electrical').put(Sig.pin['a0'])
            nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='electrical').put(Sig.pin['b0'])

            bbu.put_stub(['a0', 'b0'])

            width = width_gnd
            if (angle > 0 and radius > 0) or (angle < 0 and radius < 0):
                radius1 = radius-gap-0.5*width_gnd-0.5*width_sig
                radius2 = radius+gap+0.5*width_gnd+0.5*width_sig
            else:
                radius2 = radius-gap-0.5*width_gnd-0.5*width_sig
                radius1 = radius+gap+0.5*width_gnd+0.5*width_sig

            nd.bend(radius=radius1, angle=angle, width=width, xs=xs['a0']).\
                put(0, gap+0.5*width_gnd+0.5*width_sig)
            nd.bend(radius=radius2, angle=angle, width=width, xs=xs['a0']).\
                put(0, -gap-0.5*width_gnd-0.5*width_sig)

            if xs_bg is not None:
                nd.bend(radius=radius, angle=angle, width=xs_bg_width, xs=xs_bg).put(0)

            cfg.cp = C.pin['b0']

        return C
    return cell


def Tp_RFpadGSG(width_sig=10, width_gnd1=10, width_gnd2=None,
                gap1=10, gap2=62.5, height_tap=40,
                length_pad=75, width_pad_sig=45, width_pad_gnd=45,
                pinwidth=None, name='pad_gsg', groupname='',
                xs=None, layer=None):
    """Template for GSG RF pad cell.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length_pad')
    def cell(width_sig=width_sig, width_gnd1=width_gnd1,
             width_gnd2=width_gnd2,
             gap1=gap1, height_tap=height_tap, gap2=gap2,
             length_pad=length_pad, width_pad_sig=width_pad_sig,
             width_pad_gnd=width_pad_gnd):
        """Create and return a GSG RF cell.

        Args:
            ??

        Returns:
            Cell
        """
        with nd.Cell(name=bbu._hash_name) as C:
            C.groupname = groupname
            C.default_pins('a0', 'a0')

            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='electrical').put(0, 0, 180)
            bb_length = length_pad+height_tap
            bbu.put_stub([])
            bbu.cellname(C.cell_name, bb_length).put(-0.25*bb_length, 0, 180, 'a0')
            bbu.parameters(bbu._hash_params).put(-0.25*bb_length, 0, 180, 'a0')

            boundary = 10
            nd.Pin(name='rc', xs=None, width=None).\
                put(length_pad+height_tap+boundary, 0, 0)
            heightRectangle = length_pad+boundary
            widthRectangle  = 2*width_pad_gnd+width_sig+2*gap2+2*boundary
            widthTaper      = width_gnd1+width_gnd2+width_sig+2*gap1+2*boundary

            if width_gnd2 == None:
                width_gnd2 = width_pad_gnd

            # outline of the BB
            angle = -degrees(atan(height_tap/(0.5*(widthRectangle-widthTaper))))
            outline = geom.trapezoid(length=widthTaper, height=height_tap,
                angle1=angle, angle2=angle, position=4)
            nd.Polygon(layer='bbox', points=outline).put(0, 0, -90)

            outline = geom.box(length=heightRectangle, width=widthRectangle)
            nd.Polygon(layer='bbox', points=outline).put(height_tap, 0, 0)

            # Ground Taper
            outline = geom.tetragon(length=width_gnd1, height=height_tap,
                dx=gap2-gap1, x=width_pad_gnd, position=4)
            nd.Polygon(layer=layer, points=outline).\
                put(0, gap1+0.5*(width_sig+width_gnd1), -90)
            # Ground Pad
            outline = geom.rectangle(length=width_pad_gnd,
                    height=length_pad, position=7)
            nd.Polygon(layer=layer, points=outline).\
                put(height_tap, gap2+0.5*width_sig, -90)

            # Signal Line
            outline = geom.box(length=height_tap, width=width_sig)
            nd.Polygon(layer=layer, points=outline).put(0, 0, 0)
            # Signal Pad
            outline = geom.rectangle(length=width_pad_sig,
                    height=length_pad, position=4)
            nd.Polygon(layer=layer, points=outline).put(height_tap, 0, -90)

            # Ground Taper
            outline = geom.tetragon(length=width_gnd2, height=height_tap,
                dx=-gap2-width_pad_gnd+gap1+width_gnd2, x=width_pad_gnd,
                position=4)
            nd.Polygon(layer=layer, points=outline).\
                put(0, -gap1-0.5*(width_sig+width_gnd2), -90)
            # Ground Pad
            outline = geom.rectangle(length=width_pad_gnd, height=length_pad, position=1)
            nd.Polygon(layer=layer, points=outline).\
                put(height_tap, -gap2-0.5*width_sig, -90)
        return C
    return cell
