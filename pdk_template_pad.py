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
# (c) 2016-2019 Ronald Broeke, Katarzyna Lawniczuk
#==============================================================================

"""Module defining black box templates for PDK implementation."""


import nazca as nd
import nazca.bb_util as bbu


def Tp_DCpad_rectangle(length=100, tab_width=0,
        buf_length=10, buf_width=10,
        pinwidth=None,
        name='DCpad', groupname='', xs=None,
        metal_stub_length=2.0, dclayer=None,
        icon=None):
    """Template for a SQUARE DC pad in different xsections.

    Length and width are the actual metal size inside the BB.
    The buffers buf_x and buf_y are to provide extra BB space.

    Returns:
        function returning a Cell: pdc pad F(lengh, width)
    """
    @bbu.hashme(name, 'length', 'width')
    def cell(length=length):
        """Create and return a DCpad cell.

        Args:
            length (float): pad length in um
            width (float): pad width in um

        Returns:
            Cell: pad element
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.foundry_spt = []
            bb_length = length + 2*buf_length
            C.default_pins('c0','c0')
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='electrical').\
                put(0.5*bb_length, 0, 180)

            bbu.put_stub('c0', length=pinwidth['c0'], shape='circle')
            bbu.put_boundingbox('org', bb_length, bb_length)

            if icon:
                icon(bb_length, bb_length).put(0, 'cc')
        return C
    return cell


def Tp_DCpad_circle(diameter=100, pinwidth=None,
        buf_length=10, buf_width=10,
        name='DCpad', groupname='', xs=None,
        dclayer=None):
    """Template for DC pad in different xsections.

    Args:

    Returns:
        function: function generating a Cell
    """
    @bbu.hashme(name, 'diameter')
    def cell(diameter=diameter):
        """Create a DCpad cell.

        Args:
            diameter (float): diameter of circular debt

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            bb_width = diameter + 2*buf_width
            bb_length = diameter + 2*buf_length

            C.default_pins('c0', 'c0')
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='electrical').\
                put(0.5*bb_length, 0, 180)

            bbu.put_stub('c0', length=pinwidth['c0'], shape='circle')
            bbu.put_boundingbox('org', bb_length, bb_width)

            for lay, grow, acc, line in nd.layeriter(xs['c0']):
                pad = nd.geom.circle(radius=0.5*diameter, N=100)
                nd.Polygon(layer=lay, points=pad).\
                    put('cc')

        return C
    return cell


def Tp_DCpad_lw(length=100, width=100, pinwidth=None,
        name='DCpad_lw', xs=None, groupname='',
        metal_stub_length=2.0, dclayer=None):
    """Template for DC pad in different xsections.

    The current standard BB only provides square BBs.
    Length and width are the actual metal size.

    This is now actually a white BB, as it does not need replacement by the
    foundry. The building block has 5 micron extra space, since the design
    manual specifies 10 um minimum spacing between metal tracks/pads

    Args:
        xs (dict): {pin_name: xsection_name}
        pinwidth (dict): {pin_name: pin_width}

    Returns:
        function returning a Cell: pdc pad F(lengh, width)
    """
    @bbu.hashme(name, 'length', 'width')
    def cell(length=length, width=width):
        """Create a DCpad_lw cell.

        Args:
            length (float): length of the pad in um
            width (float): width of the pad in um

        Returns:
            Cell: dcpad element
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.default_pins('c0','c0')
            buf = 10
            bb_width = width + buf
            bb_length = length + buf
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='electrical').\
                put(0.5*bb_length, 0, 180)
            bbu.put_stub('c0', length=pinwidth['c0'], shape='circle')
            bbu.put_boundingbox('org', bb_length, bb_width)

            for lay, grow, acc, line in nd.layeriter(xs['c0']):
                pad = nd.geom.rounded_rect(
                    length=length+grow, height=width, position=5)
                nd.Polygon(layer=lay, points=pad).\
                    put(C.pin['c0'])
        return C
    return cell


def Tp_RFpad(length=100, width=100, pinwidth=None,
        name='pad_rf', groupname='', xs=None,
        icon=None):
    """Template for RF pad in different xsections.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length, width=width):
        """Create a RFpad cell.

        Returns:
            Cell
        """
        with nd.Cell(name=bbu._hash_name) as C:
            C.default_pins('c0','c0')
            C.groupname = groupname

            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='electrical').put(0, 0, 180)

            bbu.put_stub('c0')
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell


