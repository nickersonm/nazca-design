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


from math import tan, cos, radians
import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom
import nazca.cfg as cfg



def Tp_XStransition(length=33, width=30, pinwidth=None, name='xs1_to_xs2',
        groupname='', xs=None, icon=None):
    """Template for waveguide transitions between xsections.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell():
        """
        Create a XS-transition cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length)
            nd.connect_path(p1, p2, length)

            bbu.put_stub(['a0', 'b0'])
            bbu.put_boundingbox('org', length, width)
            if icon:
                icon(length, width).put(0)

            C.group = groupname
        return C

    return cell


def Tp_Isolation(length=30, width=20, pinwidth=None, name='iso',
        groupname='', xs=None):
    """Template for Isolation section in different xsection.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length):
        """Create an Isolation cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length)

            bbu.put_stub(['a0', 'b0'])
            bbu.put_boundingbox('org', length, width)

            nd.connect_path(p1, p2, length)
            C.group = groupname
        return C
    return cell


def Tp_CellBoundary(length=4600, height=4000, cleave=100,
                  coatingW='NO', coatingE='NO', name='cell_boundary',
                  groupname='', xs=None):
    """Template for a cell boundary.

    Args:

    Returns:
        function generating a Cell: die template with cleave street
    """

    @bbu.hashme(name)
    def cell(cleave=cleave, length=length, height=height,
             coatingW=coatingW, coatingE=coatingE):
        """Create a Cell_Boundary cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:

            bbu.put_boundingbox('org', length, height, outline=False, \
                move=(-0.5*cleave, 0.5*height-0.5*cleave))
            bbu.parameters('hashme').put(0, 'cc')

            #TODO: No foundry specific stuff here: remove.
            outline = [
                (0, 0), (length-cleave, 0),
                (length-cleave, height-cleave), (0, height-cleave)
                ]
            nd.Polygon(layer='Polyimide1Base', points=outline).put(0)

            for lay, grow, acc, line in nd.layeriter(xs):
                frame = geom.frame(sizew=cleave, sizel=length,
                    sizeh=height, grow=grow)
                nd.Polygon(layer=lay, points=frame).put(0)

            ##adding coating to nd.Cell boundary: west and east
            options = {'AR': 'Coating_AR', 'HR': 'Coating_HR',
                       'NO': 'Coating_NO', 'DC': 'Coating_DC'}
            coating_on_chip = nd.get_parameter('cell_coating_on_chip')
            box_coatno = [
                (0, 0), (cleave/2, 0),
                (cleave/2, height), (0, height)
                ]
            box_coatother = [
                (0, 0), (coating_on_chip+0.5*cleave, 0),
                (coating_on_chip+0.5*cleave, height), (0, height)
                ]
            pinW = nd.Pin().put(-0.5*cleave, -0.5*cleave)
            pinE = nd.Pin().put(length-0.5*cleave, height-0.5*cleave, 180)

            coatings = {'E': (coatingE, pinE), 'W': (coatingW, pinW) }

            for coat, pin in coatings.values():
                if coat in options.keys():
                    layer = options[coat]
                else:
                    print('Available coating options: AR | HR | NO | DC. Default is NO.')

                if coat == 'NO':
                    box = box_coatno
                else:
                    box = box_coatother
                nd.Polygon(layer=layer, points=box).put(pin)

            #Placing IOs
            pitch = 12.5
            amount_ios = round((height-cleave-pitch)/pitch)

            #IOs positions
            lay = 'bb_pin'
            for no in range(0, amount_ios):
                pinID = 'io{:03d}'.format(no)
                if no % 2 == 0:
                    angle = 7
                else:
                    angle = 0

                p = nd.Pin(name=pinID).put(-cleave/2, 12.5+no*pitch, angle)
                p = nd.Pin(name='io'+str(no)).put(-cleave/2, 12.5+no*pitch, angle)
                bbu.make_pincell().put(p)
                nd.text(pinID, layer=lay, height=0.15, align='rc').put(p.move(-0.1))

            for no in range(0, amount_ios):
                mo = amount_ios+no
                pinID = 'io{:03d}'.format(mo)
                if mo % 2 == 0:
                    angle = 7
                else:
                    angle = 0

                p = nd.Pin(name=pinID).put(length-cleave/2, 2*12.5+no*pitch, 180+angle)
                bbu.make_pincell().put(p)
                nd.text(pinID, layer=lay, height=0.15, align='lc').\
                    put(p.move(-0.1, 0, 180))

            C.group = groupname
        return C
    return cell


def Tp_CellBoundary2(length=5000, height=5000, cleave=100,
                     name='cell_boundary', groupname='',
                     pitch=None, xs=None):
    """Template for a cell boundary.

    Returns:
        function that generates a Cell object
    """

    @bbu.hashme(name)
    def cell(length=length, height=height, cleave=cleave, pitch=pitch):
        """Create a cell boundary.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            bbu.parameters('hashme').put(0)

            #TODO: No foundry specific stuff here: remove.
            for lay, grow, acc, line in nd.layeriter(xs):
                frame = geom.frame(sizew=cleave, sizel=length,
                    sizeh=height, grow=grow)
                nd.Polygon(layer=lay, points=frame).put(0)

            #Placing IOs
            amount_ios = round((height-cleave-pitch)/pitch)

            #IOs positions
            lay = 'AnnotationIO'
            for no in range(0, amount_ios):
                pinID = 'ioL{:03d}'.format(no)
                angle = 0

                p = nd.Pin(name=pinID).put(-cleave, pitch+no*pitch, angle)
                p = nd.Pin(name='ioL'+str(no)).put(-cleave, pitch+no*pitch, angle)
                bbu.arrow.put(p)
                nd.text(pinID, layer=lay, height=0.15, align='rc').put(p.move(-0.1))

            for no in range(0, amount_ios):
                mo = amount_ios+no
                pinID = 'ioR{:03d}'.format(no)
                angle = 0

                p = nd.Pin(name=pinID).put(length, pitch+no*pitch, 180+angle)
                p = nd.Pin(name='ioR'+str(no)).put(length, pitch+no*pitch, 180+angle)
                bbu.arrow.put(p)
                nd.text(pinID, layer=lay, height=0.15, align='lc').\
                    put(p.move(-0.1, 0, 180))

        C.group = groupname
        return C
    return cell


def Tp_CellID(width=200, length=200, ID='xx', name='Cell_ID',
        groupname='', xs=None):
    """Smart Photonics CellID.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'ID')
    def cell(ID=ID):
        """Create a Cell_ID cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            nd.Pin(name='a0', xs=None).put(0, 0, 180)

            bbu.put_boundingbox('org', length, width)
            C.group = groupname
        return C
    return cell


def Tp_IO(length=100, width=3.5, angle=-7, pinwidth=None, name='io',
        groupname='', xs=None):
    """Template trapezoidal IO.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length', 'angle')
    def cell(length=length, angle=angle):
        anglerad = radians(angle)
        """Create a trapezoidal IO cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True, instantiate=False) as C:
            """When instantiate is True, cell put in (0, 0, 0)"""
            L - 0.5*length*(1.0+1.0/cos(anglerad))
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(L)
            nd.connect_path(p1, p2, L)
            bbu.put_stub(['a0', 'b0'])
            bbu.cellname(length=0.5*length).put(C.pin['a0'].rot(180))

            temp_a = tan(anglerad)+0.5*width*tan(anglerad)
            temp_b = length/cos(anglerad)+length
            for lay, grow, acc, line in nd.layeriter(xs['a0']):
                outline = [
                    (0, 0.5*width+grow),
                    (0.5*temp_b-(width+grow)*temp_a, 0.5*width+grow),
                    (0.5*temp_b+grow*temp_a, -0.5*width-grow),
                    (0, -0.5*width-grow)
                    ]
                nd.Polygon(layer=lay, points=outline).put(0)
            cfg.cp = C.pin['b0']
        C.group = groupname
        return C
    return cell


def Tp_IOtrapezoid(length=100, width=3.5, angle=-7, pinwidth=None, name='io',
        groupname='', xs=None):
    """Template trapezoidal IO.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length', 'angle')
    def cell(length=length, angle=angle):
        anglerad = radians(angle)
        """Create a trapezoidal IO cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True, instantiate=False) as C:
            #When instantiate is True, cell put in (0,0,0)
            L = length/cos(anglerad)-0.5*width*tan(abs(anglerad))
            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(L)
            nd.connect_path(p1, p2, L)

            bbu.put_stub(['b0'])
            bbu.cellname(length=0.5*length).put(-0.1*length, 0, 180, 'a0')

            for lay, grow, acc, line in nd.layeriter(xs['a0']):
                outline = nd.geom.trapezoid(
                    length=length/cos(anglerad)+grow*tan(abs(anglerad)),
                    height=width+2*grow, angle1=90+angle, angle2=90,
                    position=2)
                nd.Polygon(layer=lay, points=outline).put(0)
            cfg.cp = C.pin['b0']
            C.group = groupname
        return C
    return cell


def Tp_IO2(length=100, width=3.5, angle=7, pinwidth=None, name='io',
        groupname='', xs=None):
    """Template trapezoidal IO, but flipped (XL).

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length', 'angle')
    def cell(length=length, angle=angle):
        """Create a trapezoidal IO cell.

        Returns:
            Cell
        """
        a = radians(angle)
        with nd.Cell(hashme=True, instantiate=False) as C:
            l = 0.5*length / cos(a) # half the length of the hypotenuse.
            """When instantiate is True, cell put in (0,0,0)"""

            p1 = nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(l)
            nd.connect_path(p1, p2, l)

            bbu.put_stub(['b0'])
            bbu.cellname(C.cell_name, -0.5*l).put(0, 'a0')

            for lay, grow, acc, line in nd.layeriter(xs['a0']):
                w = (width + 2 * grow) / 2
                x = w * tan(a)
                outline = [(l, w), (-l+x, w), (-l-x, -w), (l, -w)]
                nd.Polygon(layer=lay, points=outline).put(0)
            cfg.cp = C.pin['b0']
        C.group = groupname
        return C
    return cell


def Tp_BB1x1(length=100, width=50, pinwidth=None, name='BBname', groupname='',
             xs=None, xs_bb=None, ashift=0, bshift=0, parameters=None):
    """Template for a 1 input and 1 output Building Block.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name)
    def cell():
        """Create a BB1x1 cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.foundry_spt = []
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length, bshift)

            bbu.put_stub(['a0', 'b0'])
            bbu.put_boundingbox('org', length, width)

            cfg.cp = C.pin['b0']
            C.group = groupname
        return C
    return cell


def Tp_BB3ports(length=100, width=50, pinwidth=None, name='BBname',
            groupname='', xs=None, xs_bb=None,
            ashift=0, bshift=0, cshiftx=0, cshifty=0,
            parameters=None):
    """Template Block with 3 ports: input, output, contact.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length):
        """Create a BB3ports cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.foundry_spt = []
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length, bshift)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(length/2+cshiftx, cshifty, 90)

            bbu.put_stub(['a0', 'b0', 'c0'])
            bbu.put_boundingbox('org', length, width)
            cfg.cp = C.pin['b0']
            C.group = groupname
        return C
    return cell


def Tp_BB4ports(length=100, width=50, pinwidth=None, name='BBname',
            groupname='', xs=None, xs_bb=None,
            ashift=0, bshift=0, cshiftx=0, cshifty=0, parameters=None):
    """Template for a Block with 4 ports: input, output, 2 x contacts.

    Returns:
        function that generates a Cell object
    """
    @bbu.hashme(name, 'length')
    def cell(length=length):
        """Create a BB4ports cell.

        Returns:
            Cell
        """
        with nd.Cell(hashme=True) as C:
            C.foundry_spt = []
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            nd.Pin(name='b0', xs=xs['b0'], width=pinwidth['b0']).put(length, bshift)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(0+cshiftx, cshifty, 180)
            nd.Pin(name='c1', xs=xs['c1'], width=pinwidth['c0']).put(length-cshiftx, cshifty)

            bbu.put_stub(['a0', 'b0', 'c0', 'c1'])
            bbu.put_boundingbox('org', length, width)
            cfg.cp = C.pin['b0']
            C.group = groupname
        return C
    return cell


