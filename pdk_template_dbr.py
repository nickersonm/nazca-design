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



def Tp_DBR(
        length=250, width=15,
        pinwidth=None, xs=None,
        name='dbr', groupname='',
        pitch=0.237, duty_cycle=0.5,
        icon=None):
    """Template for DBR building block.

    This DBR has an option to set the duty cycle as free parameter.

    Args:
        length (float):
        width (float):
        pitch (float):
        duty_cycle (float):
        name (str):
        xs (str):

    Returns:
        function returning a Cell: dbr(length, pitch, duty_cycle)
    """
    @bbu.hashme(name, 'length', 'pitch', 'duty_cycle')
    def cell(length=length, pitch=pitch, duty_cycle=duty_cycle):
        """Create a DBR cell.

        Args:
            length (float): length of DBR section in um
            pitch (float): pitch (full period) of dbr in um
            duty_cycle (float): duty cycle of grating (default = 0.5)

        Returns:
            Cell: dbr element
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.foundry_spt = []
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(length/2, 5.0, 90)
            nd.connect_path(p1, p2, length)


            bbu.put_stub(['a0', 'b0', 'c0'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

        return C
    return cell


def Tp_DBR2(
        length=250, width=15,
        name='dbr', groupname='',
        pinwidth=None, xs=None,
        pitch=0.237,
        metal=True,
        icon=None):
    """Extented template for DBR building block.

    Args:
        length (float):
        width (float):
        pinwidth (dict):
        pitch (float):
        name (float):
        cshifty (float):
        metal (int):
        pad (int):
        DCmetalPad (int):
        DCmetalLine (int):
        xs (dict):
        xs_bb (str):
        xs_metal (str):

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length', 'pitch')
    def cell(length=length, pitch=pitch, metal=metal):
        """Create a DBR cell.

        Args:
            length (float):
            pitch (float):

        Returns:
            Cell
        """
        duty_cycle = 0.5
        if metal:
            contact = 1
        else:
            contact = 0
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(0.5*length, 0, 90)
            nd.connect_path(p1, p2, length)

            bbu.put_stub(['a0', 'b0', 'c0'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)
        return C
    return cell

