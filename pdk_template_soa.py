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
import pandas as pd
import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom
import nazca.cfg as cfg


def Tp_SOA(
        length=250, width=50,
        padwidth = 200,
        pinwidth=None, dy_io=45,
        name='soa', groupname='',
        xs=None, xs_bb=None,
        pads=False,
        metallayer='BB', drawmetal=False,
        icon=None):
    """Template for SOA building block.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length', 'pads')
    def cell(length=length, pads=pads):
        """Create a SOA cell.

        Args:
            length (float): SOA length in um

        Returns:
            Cell; SOA element
        """
        nonlocal width
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, dy_io, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length, dy_io, 0)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(length/2, dy_io, 90)
            nd.connect_path(p1, p2, length)

            if pads:
                bb_width = padwidth
            else:
                bb_width = width

            bbu.put_stub(['a0', 'b0', 'c0'])
            bbu.put_boundingbox('org', length, bb_width, align='lb')
            if icon:
                icon(length, bb_width).put(0, 0, 180, 'lc')
            cfg.cp = C.pin['b0']
        return C
    return cell



