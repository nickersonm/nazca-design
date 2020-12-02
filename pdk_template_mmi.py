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


import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom


OD_toStd = 1

def Tp_Modefilter(length=98.15, width=36, pinwidth=None,
        name='modefilter', groupname='', xs=None,
        icon=None, version=None, store_pins=False):
    """Template for Modefilter in different xsections.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell():
        """Create a Modefilter cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.version = version
            C.store_pins = store_pins
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='optical').put(length)
            nd.connect_path(p1, p2, length)

            bbu.put_stub(['a0', 'b0'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell


def Tp_MIR1p(length=40.2, width=35, pinwidth=None,
        name='mir1p', groupname='', xs=None, icon=None, version=None,
        store_pins=False):
    """Template for MMI1x2 in different xsections.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell():
        """Create and return a MIR_1p cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.default_pins('a0', 'a0')
            C.version = version
            C.store_pins = store_pins
            C.groupname = groupname
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, 0, 180)
            nd.connect_path(p1, p1, length)
            #nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(0, 0, 180)

            bbu.put_stub(['a0'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell


def Tp_MIR2p(length=70.6, width=39, pinwidth=None, offset=1.5,
        name='mir2p', groupname='', xs=None, icon=None, version=None,
        store_pins=False):
    """Template for MMI1x2 in different xsections.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell():
        """Create a MIR_2p cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.default_pins('a0', 'a1')
            C.version = version
            C.store_pins = store_pins
            C.groupname = groupname
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, +offset, 180)
            #nd.Pin(name='b1', xs=xs['b1'], width=pinwidth['b1']).put(0, +offset, 180)
            p2 = nd.Pin(name='a1', xs=xs['a1'], width=pinwidth['a1'], remark='optical').put(0, -offset, 180)
            #nd.Pin(name='b0', xs=xs['a0'], width=pinwidth['b0']).put(0, -offset, 180)
            nd.connect_path(p1, p1, length)
            nd.connect_path(p1, p2, length)
            nd.connect_path(p2, p2, length)
            nd.connect_path(p2, p1, length)

            bbu.put_stub(['a0', 'a1'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell


def Tp_MMI1x2(length, width, pinwidth=None, offset=0,
        name='mmi1x2', groupname='', xs=None, icon=None, version=None,
        store_pins=False):
    """Template for MMI1x2 in different xsections.

     Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell():
        """
        Create a MMI1x2 cell.

        Returns:
            Cell"""
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.version = version
            C.store_pins = store_pins
            p1 =  nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='optical').put(length, +offset)
            p3 = nd.Pin(name='b1', xs=xs['b1'], width=pinwidth['b1'], remark='optical').put(length, -offset)
            nd.connect_path(p1, p2, length)
            nd.connect_path(p1, p3, length)

            bbu.put_stub(['a0', 'b0', 'b1'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell


def Tp_MMI2x2(length, width, pinwidth=None, offset=0,
        name='mmi2x2', groupname='', xs=None, icon=None, version=None,
        store_pins=False):
    """Template for MMI2x2.

    Returns:
        function that generates a Cell object
    """
    if not isinstance(offset, list):
        offset = [offset, -offset, offset, -offset]
    @bbu.hashme(name)
    def cell():
        """
        Create a MMI2x2 cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.version = version
            C.store_pins = store_pins
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, offset[0], 180)
            p2 = nd.Pin(name='a1', xs=xs['a1'], width=pinwidth['a1'], remark='optical').put(0, offset[1], 180)
            p3 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='optical').put(length, offset[2])
            p4 = nd.Pin(name='b1', xs=xs['b1'], width=pinwidth['b1'], remark='optical').put(length, offset[3])
            nd.connect_path(p1, p3, length)
            nd.connect_path(p1, p4, length)
            nd.connect_path(p2, p3, length)
            nd.connect_path(p2, p4, length)

            bbu.put_stub(['a0', 'a1', 'b0', 'b1'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell



