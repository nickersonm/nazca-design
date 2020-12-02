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
# Test replacement of gds cells from one file with those from another.
# This will be used for replacing black-box building blocks with their real
# implementation.
#
# (c) 2016-2018 Katarzyna Lawniczuk, Ronald Broeke
#==============================================================================
"""
Module to create icons for building blocks from template functions.
"""

import nazca as nd
import nazca.geometries as geom
from nazca.logging import logger

scale_length = 0.7
scale_width = 0.7
max_aspect_ratio = 5

def calc_buf(length, width, bufx=None, bufy=None):
    """Calculate icon size and buffer between the icon and the BB edge.

    If no buffer value is provided (None) the icon will be scaled to the
    size of the BB based on module parameters 'scale_length' and 'scale_width'.

    If a positive buffer value is provided the icon will adjust to accomodate
    this buffer between the icon and BB at both sides of the icon.

    If a negative buffer value is provided it is interpreted as the absolute
    size of the icon. Hence, bufy=-10 results in an icon of size 10 in the
    y-direction.

    Args:
        length (float): bb-length (x-direction)
        width (float): bb_width (y-direction)
        bufx (float): space between the icon and the bb-edge in the x-direction. (default = None)
        bufy (float): space between the icon and the bb-edge in the y-direction. (default = None)

    Returns:
        float, float, float, float: icon-length (x), icon-width (y), side-buffer x, side-bufferf y
    """
    if bufx is None:
        bufx = 0.5*length*(1-scale_length)
    elif bufx < 0:
        bufx = 0.5*(length+bufx)

    if bufy is None:
        bufy = 0.5*width*(1-scale_width)
    elif bufy < 0:
        bufy = 0.5*(width+bufy)

    return length-2*bufx, width-2*bufy, bufx, bufy


def Tp_icon_mmi(Nin=2, Nout=2, bufx=None, bufy=None, layer=None, move=(0,0,0)):
    """Template to draw icons for a NxM MMI shape.

    Args:
        Nin (int): number of input guides
        Nout (int): number of output guides
        layer (int | tuple): layer to draw icon in

    Returns:
       function: function returning a Polygon with MMI shape
    """
    def icon_mmi(length=None, width=None, bufx=None, bufy=None):
        """Create a MMI icon.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)

        Nin0 = 0.5*(Nin-1)
        Nout0 = 0.5*(Nout-1)
        w_mmi = scale_width*width
        w_in = w_mmi/(2.0*max(Nin+0.5, Nout+0.5))
        w_pitch = 2*w_in
        l_in = w_in

        p = [] # list of points (x, y)
        p.append((bufx+l_in, -0.5*w_mmi)) #bottom left
        for i in range(Nin):
            p.append((bufx+l_in, (-Nin0+i)*w_pitch - 0.5*w_in))
            p.append((bufx, (-Nin0+i)*w_pitch - 0.5*w_in))
            p.append((bufx, (-Nin0+i)*w_pitch + 0.5*w_in))
            p.append((bufx+l_in, (-Nin0+i)*w_pitch + 0.5*w_in))
        p.append((bufx+l_in, + 0.5*w_mmi))
        p.append((bufx+length-l_in, +0.5*w_mmi))
        for i in range(Nout-1, -1, -1):
            p.append((bufx+length-l_in, (-Nout0+i)*w_pitch + 0.5*w_in))
            p.append((bufx+length, (-Nout0+i)*w_pitch + 0.5*w_in))
            p.append((bufx+length, (-Nout0+i)*w_pitch - 0.5*w_in))
            p.append((bufx+length-l_in, (-Nout0+i)*w_pitch - 0.5*w_in))
        p.append((bufx+length-l_in, -0.5*w_mmi))

        with nd.Cell('icon', instantiate=False) as icon:
            nd.Polygon(points=p, layer=layer).put(0)
            x, y, a = move
            nd.Pin('cc').put(x+0.5*length+bufx, y, a+180)
        return icon
    return icon_mmi


def Tp_icon_mir(Nin=2, bufx=None, bufy=None, layer=None, move=(0,0,0)):
    """Template to draw icons for a N-port MIR shape.

    Args:
        Nin (int): number of input guides
        Nout (int): number of output guides
        layer (int | tuple): layer to draw icon in

    Returns:
       function: function returning a Polygon with MIR shape
    """
    def icon_mir(length, width, bufx=None, bufy=None):
        """Create a N-port MIR icon.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length2, width2, bufx, bufy = calc_buf(length, width, bufx, bufy)

        w_mmi = scale_width*width2
        w_in = w_mmi/(2.0*Nin+0.5)
        w_pitch = 2*w_in
        l_in = w_in
        Nin0 = 0.5*(Nin-1)

        p = []
        p.append((bufx+l_in, -0.5*w_mmi)) #bottom left
        for i in range(Nin):
            p.append((bufx+l_in, (-Nin0+i)*w_pitch - 0.5*w_in))
            p.append((bufx, (-Nin0+i)*w_pitch - 0.5*w_in))
            p.append((bufx, (-Nin0+i)*w_pitch + 0.5*w_in))
            p.append((bufx+l_in, (-Nin0+i)*w_pitch + 0.5*w_in))
        p.append((bufx+l_in, +0.5*w_mmi))
        p.append((bufx+length2-0.5*w_mmi, +0.5*w_mmi))
        p.append((bufx+length2, 0))
        p.append((bufx+length2-0.5*w_mmi, -0.5*w_mmi))
        with nd.Cell('icon', instantiate=False) as icon:
            nd.Polygon(points=p, layer=layer).put(0)
            x, y, a = move
            nd.Pin('cc').put(x+0.5*length, y, a+180)
        return icon
    return icon_mir


def Tp_icon_strt(bufx=0, bufy=0, layer=None, move=(0,0,0)):
    """Template to draw straight waveguide icons.

    Returns:
        function
    """
    def icon_strt(length, width, bufx=bufx, bufy=bufy):
        """Create a straight waveguide icon.

        Args:
            length (float): length of the icon space
            width (float): width of the icon icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        with nd.Cell('icon', instantiate=False) as icon:
            rect = geom.box(length, width)
            nd.Polygon(points=rect, layer=layer).put(0)
            x, y, a = move
            nd.Pin('cc').put(x+0.5*length, y, a)
        return icon
    return icon_strt


def icon_strt(length, width, layer, bufx=None, bufy=None):
    length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
    with nd.Cell('icon', instantiate=False) as icon:
        rect = geom.box(length, width)
        nd.Polygon(points=rect, layer=layer).put(0)
        nd.Pin('cc').put(0.5*length)
    return icon


def Tp_xsection_transition(layer1=None, layer2=None, move=(0,0,0)):
    """Template to draw waveguide transistion icons.

    Returns:
        function
    """
    def xsection_transition(length, width):
        """Create a xsection transition icon.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space

        Returns:
            Cell: icon
        """
        with nd.Cell('icon', instantiate=False) as icon:
            rect = geom.box(0.25*length, 0.50*width)
            nd.Polygon(points=rect, layer=layer1).put(0.25*length)
            nd.Polygon(points=rect, layer=layer2).put(0.50*length)
            x, y, a = move
            nd.Pin('cc').put(x+0.5*length, y, a+180)
        return icon
    return xsection_transition


def xsection_transition(length, width, layer1, layer2, move=(0,0,0)):
    """Create an icon for an xs transition.

    Returns:
        Cell: icon
    """
    with nd.Cell('icon', instantiate=False) as icon:
        rect = geom.box(0.25*length, 0.50*width)
        nd.Polygon(points=rect, layer=layer1).put(0.25*length)
        nd.Polygon(points=rect, layer=layer2).put(0.50*length)
        x, y, a = move
        nd.Pin('cc').put(x+0.5*length, y, a+180)
    return icon


def Tp_icon_rounded_pad(bufx=None, bufy=None, layer=None, move=(0,0,0)):
    """Template to draw rounded pad icons.

    Returns:
        function
    """
    def icon_pad(length, width, bufx=bufx, bufy=bufy):
        """Create an icon of a rounded pad.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        with nd.Cell('icon', instantiate=False) as icon:
            pad = geom.rounded_rect(
                length=length, height=width, position=5)
            nd.Polygon(layer=layer, points=pad).put(0)
            x, y ,a, = move
            nd.Pin('cc').put(x, y, a)
        return icon
    return icon_pad


def Tp_icon_circ_pad(bufx=None, bufy=None, layer=None, move=(0,0,0)):
    """Template to draw circular pad icons.

    Returns:
        function
    """
    def icon_pad(length, width, bufx=bufx, bufy=bufy):
        """Create an icon of a circular pad.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        with nd.Cell('icon', instantiate=False) as icon:
            pad = geom.circ(diameter=min(length, width))
            nd.Polygon(layer=layer, points=pad).put(0)
            x, y ,a, = move
            nd.Pin('cc').put(x, y, a)
        return icon
    return icon_pad

def Tp_icon_directional_coupler(bufx=None, bufy=None, layer1=None, layer2=None,
        move=(0,0,0)):
    """Template for functions to draw directional coupler icons.

    Args:
        layer1 (): layer-1 to draw the icon in
        layer2 (): layer-2 to draw the icon in

    Returns:
        function
    """
    if layer2 is None:
        layer2 = layer1
    def coupler(length, width, bufx=bufx, bufy=bufy):
        """Create an icon with a directional coupler.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        w = 0.05*width
        R = min(0.7*width, 0.2*length)
        Lr = 2*1.41*R
        a = 45.0
        L = length - Lr

        bend = nd.Tp_arc(radius=R, width=w, layer=layer1)
        with nd.Cell('sbend', instantiate=False) as part:
            nd.strt(length=0.5*L, width=w, layer=layer2).put(0)
            bend(angle=-a).put()
            bend(angle=a).put()

        with nd.Cell('icon', instantiate=False) as icon:
            part.put(0, -w)
            part.put(0, -w, 180, flip=True)

            part.put(0, w, 0, flip=True)
            part.put(0, w, 180)

            x, y, a = move
            nd.Pin('cc').put(x, y, a)
        return icon
    return coupler


def Tp_icon_MZI(bufx=None, bufy=None, layer=None, move=(0,0,0)):
    """Template for functions to draw MZI icons.

    Args:
        layer (): layer to draw the icon in

    Returns:
        function
    """
    def mzi(length, width, bufx=bufx, bufy=bufy):
        """Create an icon with a directional coupler.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        sep = 0.4*width
        w = 0.1*width
        R = min(0.4*width, 0.2*length)
        Lr = 2*1.41*R
        a = 45.0
        L = length - Lr
        Lmmi = 0.05*L
        Wmmi = 5*w

        bend = nd.Tp_arc(radius=R, width=w, layer=layer)
        with nd.Cell('sbend', instantiate=False) as part:
            nd.strt(length=0.4*L, width=w, layer=layer).put(0)
            bend(angle=a).put()
            b = bend(angle=-a).put()
            bend(angle=-a).put(nd.cp.move(Lmmi))
            bend(angle=a).put()

        with nd.Cell('icon', instantiate=False) as icon:
            part.put(0, -sep)
            part.put(0, -sep, 180, flip=True)
            nd.strt(length=Lmmi, width=Wmmi, layer=layer).put(b.pin['b0'].x, 0)

            part.put(0, sep, 0, flip=True)
            part.put(0, sep, 180)
            nd.strt(length=Lmmi, width=Wmmi, layer=layer).put(-b.pin['b0'].x, 0, 180)

            x, y, a = move
            nd.Pin('cc').put(x, y, a)
        return icon
    return mzi


def Tp_icon_Yjunction(bufx=None, bufy=None, layer=None, move=(0,0,0)):
    """Template for functions to draw Y-junction icons.

    Args:
        layer (): layer to draw the icon in

    Returns:
        function
    """
    def Yjunction(length, width, bufx=bufx, bufy=bufy):
        """Create cell with a Y-junction icon.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        w = 0.05*width
        R = min(0.7*width, 0.2*length)
        Lr = 1.41*R
        a = 45.0
        L = length - Lr

        bend = nd.Tp_arc(radius=R, width=w, layer=layer)
        strt = nd.strt(length=0.5*L, width=w, layer=layer)
        with nd.Cell('icon', instantiate=False) as icon:
            s = strt.put(0)

            bend(angle=-a).put(s)
            bend(angle=a).put()
            strt.put()

            bend(angle=a).put(s)
            bend(angle=-a).put()
            strt.put()

            x, y, a = move
            nd.Pin('cc').put(x+0.5*(Lr+L), y, a+180)
        return icon
    return Yjunction


def Tp_icon_diode(layer=None, move=(0,0,0)):
    """Template for functions to draw photo-diode icons.

    Args:
        layer (): layer to draw the icon in

    Returns:
        function
    """
    def icon_diode(length=0, width=None, bufx=None, bufy=None):
        """Create an icon with a diode symbol.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: diode icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        ratio = 0.3
        w_line = 0.1*ratio*width
        with nd.Cell('icon', instantiate=False) as icon:
            dio1 = [(-0.5*ratio*width, 0.5*ratio*width),
                    (0.5*ratio*width, 0.5*ratio*width), (0, -0.5*ratio*width)]
            dio2 = geom.rectangle(ratio*width, w_line, position=5)
            dio3 = geom.rectangle(w_line, width, position=5)
            nd.Polygon(points=dio1, layer=layer).put(0.5*length)
            nd.Polygon(points=dio2, layer=layer).put(0.5*length, -0.5*ratio*(width-w_line))
            nd.Polygon(points=dio3, layer=layer).put(0.5*length)
            x, y, a = move
            nd.Pin('cc').put(x+0.5*length, y, a)
        return icon
    return icon_diode


def Tp_icon_diode_gsg(layer=None, move=(0,0,0)):
    """Template to draw GSG diode icons.

    Args:
        layer (): layer to draw the icon in

    Returns:
        function
    """
    def icon_gsg(length=None, width=None, angle=0, pitch=100):
        """Create an icon with a gsg diode symbol.

        The icon scales with the pitch and the height.

        Args:
            length (float): not used. Only for call compatibility reasons
            width (float): optional width of the icon. (default=None)
            pitch (float): pitch between ground and signal

        Returns:
            Cell: diode gsg icon
        """
        if int((angle/90)) % 2 == 1:
            length, width = width, length

        if width is None:
            width = 2*pitch

        height = width
        ratio = 0.3
        w_line = 0.1*ratio*height
        radius = 0.15*ratio*height
        height = height-2*radius-2*w_line

        with nd.Cell('icon_gsg', instantiate=False) as icon:
            # define shapes
            outline = [(-0.5*ratio*height, 0.5*ratio*height),
                (0.5*ratio*height, 0.5*ratio*height), (0, -0.5*ratio*height)]
            diode_triangle =  nd.Polygon(points=outline, layer=layer)

            outline = geom.rectangle(ratio*height, w_line, position=5)
            diode_base =  nd.Polygon(points=outline, layer=layer)

            outline = geom.ring(radius=radius, width=w_line, N=20)
            ring = nd.Polygon(points=outline, layer=layer)

            outline = geom.rectangle(w_line, height, position=5)
            pole = nd.Polygon(points=outline, layer=layer)

            outline = geom.rectangle(2*pitch, w_line, position=5)
            gnd = nd.Polygon(points=outline, layer=layer)

            # put shapes
            gnd.put(0, -0.5*(height-w_line))
            diode_triangle.put(0)
            diode_base.put(0, -0.5*ratio*(height-w_line))
            pole.put(0)
            pole.put(pitch-0.5*w_line)
            pole.put(-pitch+0.5*w_line)
            ringpos = 0.5*(height+radius+w_line)
            ring.put(0, ringpos)
            ring.put(pitch-0.5*w_line, ringpos)
            ring.put(-pitch+0.5*w_line, ringpos)

            # add pins
            x, y, a = move
            nd.Pin('cc').put(x, y, a-angle)
            nd.Pin('top').put(0, ringpos+radius+w_line, -90)
        return icon
    return icon_gsg


def Tp_icon_grating(layer=None, move=(0,0,0)):
    """Template to grating icons.

    Args:
        layer (): layer to draw the icon in

    Returns:
        function
    """
    def icon_grating(length=100, width=20, bufx=None, bufy=None):
        """Create an icon for a grating.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: grating icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)
        Nmin = 3
        h = width
        w = 0.4*width

        if length < 2*Nmin*w:
            w = length/6
            N = 3
        else:
            N = int(length/(2*w))

        with nd.Cell('icon', instantiate=False) as icon:
            t, u, v = w/4, w/2, h/2
            cross = [(u,t),(u,v),(-u,v),(-u,t),(-w,t),(-w,-t),(-u,-t),
                    (-u,-v),(u,-v),(u,-t),(w,-t),(w,t)]
            for i in range(N):
                nd.Polygon(points=cross, layer=layer).put(0+i*w*2)
            x, y, a = move
            nd.Pin('cc').put(x+(N-1)*w, y, a)
        return icon
    return icon_grating


def Tp_icon_ssc(layer=None, move=(0,0,0)):
    """Template to draw spot-size convertor icons.

    Args:
        layer (): layer to draw the icon in

    Returns:
        function
    """
    def icon_ssc(length=100, width=20, bufx=None, bufy=None, angle=0):
        """Create an icon for a SSC.

        Args:
            length (float): length of the icon space
            width (float): width of the icon space
            bufx (float): buffer size in x-direction on one side
            bufy (float): buffer size in y-direction on one side

        Returns:
            Cell: SSC icon
        """
        length, width, bufx, bufy = calc_buf(length, width, bufx, bufy)

        # define the SSC shape polygon:
        N = 20
        dx = length/N
        wout = width
        win = 0.1*wout
        wi = 0.5*win
        wo = 0.5*(wout-win)
        x, y1, y2 = [], [], []
        for i in range(N):
            x.append(i*dx)
            dy = wi+(wo/(N-1))*i*i/(N-1)
            y1.append(dy)
            y2.append(-dy)
        X = x + x[::-1]
        Y = y1 + y2[::-1]
        outline = list(zip(X, Y))

        with nd.Cell('icon', instantiate=False) as icon:
            nd.Polygon(points=outline, layer=layer).put(0)
            x, y, a = move
            nd.Pin('cc').put(x+0.5*length, y, a+180-angle)
            nd.Pin('a0').put(0, 0, 180)
            nd.Pin('b0').put(length, 0, 0)
        return icon
    return icon_ssc

