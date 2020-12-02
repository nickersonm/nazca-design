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
#import pandas as pd
import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom
import nazca.cfg as cfg


def Tp_EOPM_DC(
        length=100, width=50,
        name='BBname', groupname='',
        pinwidth=None, xs=None,
        contacts=3, pads=True,
        icon=None):

    """Template for a EOPM_DC.

    The EOMP has 1 input, 1output, and 3 contact ports.
    Metal pads to contact ports are optional.

    Args:


    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length', 'contacts', 'pads',)
    def cell(length=length, contacts=contacts):
        """Create an electro-optical phase modulator (EOPM) cell.

        Their can be from 0 up to <contacts> contac points.
        They will be equidistantly spaced along the modulator.

        Args:
            length (float): length of the EOPM.
            contacts (int): number of contact
            pads (bool): draw pads

        Returns:
            Cell: eopm
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname

            # modulator section:
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='optical').put(length)
            nd.connect_path(p1, p2, length)

            # x-positions of pads:
            L = 0.5*(length - pinwidth['c0'])
            xpos = []
            if contacts > 1:
                for n in range(contacts):
                    xpos.append(0.5*pinwidth['c0'] + L*2*n/(contacts-1))
            elif contacts == 1:
                xpos = [0.5*length]
            elif contacts == 0:
                pass

            for n, x in enumerate(xpos):
                pinname = 'c'+str(n)
                p1 = nd.Pin(name=pinname, xs=xs['c0'], width=pinwidth['c0']).\
                    put(x, 0, 90)
                nd.put_stub(pinname)

            bbu.put_stub(['a0', 'b0'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

            cfg.cp = C.pin['b0']
        return C
    return cell


def Tp_EOPM_RF(
        length=100, width=50, rfpitch=10,
        pinwidth=None, xs=None,
        name='BBname', groupname='',
        dx_metal=0, dy_metal=0,
        pads=False, padangle1=-70, padangle2=70, RFpadType=None,
        bend_gsg_function=None, bend_gsg_params=None, icon=None):
    """Template for a EOPM_RF.

    Has optica linput, output, and 2 GSG ports.
    Option to add RF GSG pads to GSF ports is added.

    Args:
        ??

    Returns:
        function returning a cell: eopm(length, pad, padangle1, padangle2)
    """
    @bbu.hashme(name, 'length', 'pads')
    def cell(length=length, pads=pads):
        """Create a EOPM_RF cell.

        Args:
            length (float): modulator length in um
            pads (bool): draw RF pads
            padangle1 (float): angle of left RF pad in degrees w.r.t. modulator
            padangle2 (float): angle of right RF pad in degrees w.r.t. modulator

        Returns:
            Cell: eopm
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname

            wmet=pinwidth['c0']
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0'], remark='optical').put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0'], remark='optical').put(length, 0)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0'], remark='gnd').put(dx_metal, dy_metal+rfpitch+wmet, 180)
            nd.Pin(name='c1', xs=xs['c1'], width=pinwidth['c1'], remark='signal').put(dx_metal, dy_metal, 180)
            nd.Pin(name='c2', xs=xs['c2'], width=pinwidth['c2'], remark='gnd').put(dx_metal, dy_metal-rfpitch-wmet, 180)
            nd.Pin(name='d0', xs=xs['d0'], width=pinwidth['d0'], remark='gnd').put(length-dx_metal, dy_metal+rfpitch+wmet)
            nd.Pin(name='d1', xs=xs['d1'], width=pinwidth['d1'], remark='signal').put(length-dx_metal, dy_metal)
            nd.Pin(name='d2', xs=xs['d2'], width=pinwidth['d2'], remark='gnd').put(length-dx_metal, dy_metal-rfpitch-wmet)
            nd.connect_path(p1, p2, length)


            bbu.put_stub(['a0', 'b0'])
            bbu.put_boundingbox('org', length, width)

            if icon:
                icon(length, width).put(0)

            if pads:
                bend_gsg_function(angle=padangle1, **bend_gsg_params).put(C.pin['c1'])
                RFpadType.put()
                bend_gsg_function(angle=padangle2, **bend_gsg_params).put(C.pin['d1'])
                RFpadType.put()
            else:
                bbu.put_stub(['c0', 'c1', 'c2', 'd0', 'd1', 'd2'])

            cfg.cp = C.pin['b0']

        return C
    return cell
