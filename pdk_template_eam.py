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
import nazca.geometries as geom
import nazca.cfg as cfg


def Tp_EAM(length=100, width=50,
           name='BBname', groupname='',
           pinwidth=None, xs=None,
           icon=None):
    """Template for a EAM.

    Has optical input and output port optional GSG pad.

    Args:
        ??

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length):
        """Create a EAM cell.

        Args:
            length (float): length of the eam

        Returns:
            Cell: eam element
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 =nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length, 0)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(0.5*length, 0, 90)
            nd.connect_path(p1, p2, length)


            bbu.put_stub()
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell