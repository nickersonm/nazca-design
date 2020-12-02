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



from math import atan, tan, sin, cos, degrees, radians

import nazca as nd
import nazca.bb_util as bbu



def Tp_PhotoDetector(
        length=100, width=50, buf=20,
        name='BBname', groupname='',
        pinwidth=None, xs=None,
        icon=None):
    """Template for a photodetector with option to draw a metal DC pad.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length):
        """Create a PhotoDetector cell.

        Args:
            length (float): length of the diode in um
            pad (bool): flag to add a bond pad

        Returns:
            Cell: photo diode element
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.default_pins('a0', 'a0')
            C.foundry_spt = []
            bb_length = length+buf
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, 0, 180)
            p2 = nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='electrical').put(bb_length)
            #nd.connect_path(p1, p2=None, sigtyp=None, path=0)
            nd.connect_path(p1, p2, s=0, sigtype='o2e', path=0)

            bbu.put_stub(['a0', 'c0'])
            bbu.put_boundingbox('org', bb_length, width)
            if icon:
                icon(bb_length, width).put(0)
        return C
    return cell


def Tp_PhotoDetectorRF(
        length=100, width=50,
        name='BBname', groupname='',
        pinwidth=None, xs=None,
        spaceGS=10,
        icon=None):
    """Template for a photodetector with an option to draw a metal RF GSG pad.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell(length=length):
        """Create a photodetector RF cell.

        Returns:
            Cell
        """
        cshift = 0
        with nd.Cell(hashme=True) as C:
            C.default_pins('a0', 'a0')
            C.groupname = groupname
            C.foundry_spt = []
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').\
                put(0, 0, 180)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='gnd').\
                put(length, cshift-spaceGS-pinwidth['c0'])
            nd.Pin(name='c1', xs=xs['c1'], width=pinwidth['c1'], remark='signal').\
                put(length, cshift)
            nd.Pin(name='c2', xs=xs['c2'], width=pinwidth['c2'], remark='gnd').\
                put(length, cshift+spaceGS+pinwidth['c2'])
            nd.connect_path(p1, None, 0)


            bbu.put_stub(['a0'])
            bbu.put_boundingbox('org', length, width)
            if icon:
                icon(length, width).put(0)
            bbu.put_stub(['c0', 'c1', 'c2'])
        return C
    return cell


