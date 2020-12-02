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
#
# (c) 2017 Xaveer Leijtens
#-----------------------------------------------------------------------

"""
Define extra structures for a mask: alignment markers, nonius, compass,
etc.

General elements, not designed specifically for photonics.
"""
import nazca as nd
import nazca.geometries as geom

def compass(size=200, layer=1):
    """Show compass image to identify the sample direction.

    Args:
        size (double): size (wxh) of the bounding box.
        layer (int | str | tuple): layer (or list of layers) to draw in.

    Returns:
        cell containing the compass.
    """
    s = size/182
    quart = [(0,0), (s,0), (s,s), (16*s,16*s), (s/2,70*s), (-s/2,70*s),
            (-16*s,16*s), (-15*s,15*s), (-s,65*s), (0,65*s)]
    p = 75 * s
    h = 20 * s

    with nd.Cell("compass_"+nd.md5(layer)) as C:
        for lay in nd.make_iter(layer):
            nd.Polygon(layer=lay, points=quart).put(0,0,0)
            nd.Polygon(layer=lay, points=quart).put(0,0,90)
            nd.Polygon(layer=lay, points=quart).put(0,0,180)
            nd.Polygon(layer=lay, points=quart).put(0,0,270)
            nd.text('E', layer=lay, height=h, align='lc').put(p, 0)
            nd.text('N', layer=lay, height=h, align='cb').put(0, 71*s)
            nd.text('W', layer=lay, height=h, align='rc').put(-p, 0)
            nd.text('S', layer=lay, height=h, align='ct').put(0, -p)
    return C


def north(size=100, layer=1):
    """Show north pointer to identify the sample direction.

    Args:
        size (double): size (wxh) of the bounding box.
        layer (int | str | tuple): layer (or list of layers) to draw in.

    Returns:
        cell containing the pointer.
    """
    s = size/148

    pointer = [(-69*s,-74*s), (-70*s,-73*s), (-s,74*s), (s,74*s),
            (70*s,-73*s), (69*s,-74*s), (0,-33*s), (0,57*s), (-s,57*s),
            (-52*s,-52*s), (-51*s,-53*s), (0,-23*s), (0,-33*s)]

    with nd.Cell("north_"+nd.md5(layer)) as C:
        for lay in nd.make_iter(layer):
            nd.Polygon(layer=lay, points=pointer).put(0,0,0)
    return C


def nonius(layer=1):
    """Nonius structure"""
    with nd.Cell(name='nonius_'+nd.md5(layer)) as C:
        for lay in nd.make_iter(layer):
            l = [20, 20, 30, 20, 20, 20, 20, 30, 20, 20]
            for i in range(10):
                nd.Polygon(layer=lay, points=
                    nd.geometries.rectangle(l[i], 5.9-0.2*i, position=2)).put(
                        0, -45.05+10*i)
            l = [30, 20, 20, 20, 20, 30, 20, 20, 20, 20, 30]
            for i in range(11):
                nd.Polygon(layer=lay, points=
                    nd.geometries.rectangle(l[i], 6-0.2*i, position=8)).put(
                        0, -50+10*i)
    return C


#==============================================================================
#  Alignment marks
#==============================================================================
def marker1(layera=1, layerb=None):
    """TU/e-Smart alignment marker 1.

    This is the first (marker1) of two matching markers.

    Args:
        layera (int | str | tuple): layer (or list of layers) in which
            marker pattern is written.
        layerb (int | str | tuple): layer (or list of layers) in which the
            background is defined.

    Returns:
        function: marker cell.

    Example:
        Place marker centered at (0,0) in layer 1 with background layer 12::

            import nazca as nd

            m1 = nd.marker1(layera=1, layerb=12)
            m1.put()

            nd.export_plt()
    """
    with nd.Cell(name='marker1_'+nd.md5((layera,layerb))) as C:
        poly = (((-75,55), (-75,75), (75,75), (75,55), (65,55),
            (65,65), (-65,65), (-65,55),), ((-37.5,17.5),
                (-37.5,32.5), (-52.5,32.5), (-52.5,37.5),
                (-37.5,37.5), (-37.5,52.5), (-32.5,52.5),
                (-32.5,37.5), (-17.5,37.5), (-17.5,32.5),
                (-32.5,32.5), (-32.5,17.5),), ((32.5,17.5),
                (32.5,32.5), (17.5,32.5), (17.5,37.5), (32.5,37.5),
                (32.5,52.5), (37.5,52.5), (37.5,37.5), (52.5,37.5),
                (52.5,32.5), (37.5,32.5), (37.5,17.5),),
                ((-38.5,-53.5), (-38.5,-38.5), (-53.5,-38.5),
                (-53.5,-31.5), (-38.5,-31.5), (-38.5,-16.5),
                (-31.5,-16.5), (-31.5,-31.5), (-16.5,-31.5),
                (-16.5,-38.5), (-31.5,-38.5), (-31.5,-53.5),),
                ((31.5,-53.5), (31.5,-38.5), (16.5,-38.5),
                (16.5,-31.5), (31.5,-31.5), (31.5,-16.5),
                (38.5,-16.5), (38.5,-31.5), (53.5,-31.5),
                (53.5,-38.5), (38.5,-38.5), (38.5,-53.5),),
                ((-75,-75), (-75,-55), (-65,-55), (-65,-65),
                (65,-65), (65,-55), (75,-55), (75,-75),))
        nd.Pin(name='a0', xs=None).put(0,0,180)
        nd.Pin(name='b0', xs=None).put(0,0,0)
        for p in poly:
            for lay in nd.make_iter(layera):
                nd.Polygon(layer=lay, points=p).put()
        for lay in nd.make_iter(layerb):
            nd.Polygon(layer=lay,
                points=geom.rectangle(180, 180, position=5)).put()
    return C

def marker2(layera=1, layerb=None, layerc=None):

    """TU/e-Smart alignment marker 2.

    This is the second (marker2) of two matching markers.

    Args:
        layera (int | str | tuple): layer (or list of layers) in which
            marker pattern is written. Default 1.
        layerb (int | str | tuple): layer (or list of layers) in which the
            background is defined. Default None.
        layerc (int | str | tuple): layer (or list of layers) in which an
            extra box around the marker is drawn for visibility in dark
            field masks. Default None.

    Returns:
        function: function generating a cell with this specific marker.

    Example:
        Place marker centered at (0,0) in layer 5 with darkfield box::

            import nazca as nd

            m2 = nd.marker2(layera=5, layerc=5)
            m2.put()

            nd.export_plt()
    """
    with nd.Cell(name='marker2_'+nd.md5((layera,layerb,layerc))) as C:
        poly = (((-55,-75), (-55,-55), (-75,-55), (-75,55),
            (-55,55), (-55,75), (-5,75), (-5,55), (-40,55),
            (-40,40), (-55,40), (-55,30), (-40,30), (-40,15),
            (-30,15), (-30,30), (-15,30), (-15,40), (-30,40),
            (-30,55), (-5,55), (-5,-16), (-39,-16), (-39,-31),
            (-54,-31), (-54,-39), (-39,-39), (-39,-54), (-31,-54),
            (-31,-39), (-16,-39), (-16,-31), (-31,-31), (-31,-16),
            (-5,-16), (-5,-75),), ((30,15), (30,30), (15,30),
            (15,40), (30,40), (30,55), (40,55), (40,40), (55,40),
            (55,30), (40,30), (40,15),), ((31,-54), (31,-39),
            (16,-39), (16,-31), (31,-31), (31,-16), (39,-16),
            (39,-31), (54,-31), (54,-39), (39,-39), (39,-54),),
            ((55,-75), (55,-55), (75,-55), (75,-75),), ((55,55),
            (55,75), (75,75), (75,55),))
        nd.Pin(name='a0', xs=None).put(0,0,180)
        nd.Pin(name='b0', xs=None).put(0,0,0)
        for p in poly:
            for lay in nd.make_iter(layera):
                nd.Polygon(layer=lay, points=p).put()
        for lay in nd.make_iter(layerb):
            nd.Polygon(layer=lay,
                points=geom.rectangle(180, 180, position=5)).put()
        for lay in nd.make_iter(layerc):
            nd.Polygon(layer=lay,
                points=geom.frame(60, 240,240)).put(-90,-90)
    return C


#==============================================================================
#  Fiducials for machine vision
#==============================================================================

def target(layera=1, layerb=None):

    """Fiducial marker 'target' for machine vision.

    Args:
        layera (layer): layer in which the target is written.
        layerb (layer): layer in which etch background is defined (trench).
        name (string): cell name.

    Returns:
        Cell: cell with this marker.

    Example:
        Place marker centered at (0,0) in layer 1 with background
        layer 12.

            import nazca as nd

            f1 = nd.target(layera=1, layerb=12)
            f1.put()

            nd.export_plt()
    """
    with nd.Cell(name='Fiducial_target_'+nd.md5((layera,layerb))) as C:
        nd.Polygon(layer=layera,
                points=geom.rectangle(200, 10, position=5)).put()
        nd.Polygon(layer=layera,
                points=geom.rectangle(10, 200, position=5)).put()
        nd.Polygon(layer=layera,
                points=geom.ring(radius=35, width=10, N=41)).put()
        nd.Polygon(layer=layera,
                points=geom.ring(radius=70, width=10, N=81)).put()
        nd.Polygon(layer=layerb,
                points=geom.rectangle(225, 225, position=5)).put()
        nd.Pin(name='a0', xs=None).put(0,0,180)
        nd.Pin(name='b0', xs=None).put(0,0,0)
    return C

def cornerUL(layera=1, layerb=None):

    """Fiducial marker ┌ with upper left corner for machine vision.

    This marker is not symmetric and can be used to specify the orientation
    of the chip.

    Args:
        layer1 (int | str | tuple): layer (or list of layers) in which the
            shape is written.
        layer2 (int | str | tuple): layer (or list of layers) in which the
            etch background is defined.

    Returns:
        Cell: cell with this marker.

    Example:
        Place marker centered at (0,0) in layer 1 with background layer 12.

            import nazca as nd

            f1 = nd.cornerUL(layera=1, layer2=12)
            f1.put()

            nd.export_plt()
    """
    with nd.Cell(name='Fiducial_cornerUL_'+nd.md5((layera,layerb))) as C:
        for lay in nd.make_iter(layera):
            nd.Polygon(layer=lay,
                    points=geom.rectangle(10, 125, position=5)).\
                            put(-67.5,0)
            nd.Polygon(layer=lay,
                    points=geom.rectangle(125, 10, position=5)).\
                            put(0,67.5)
        for lay in nd.make_iter(layerb):
            nd.Polygon(layer=lay,
                    points=geom.rectangle(200, 200, position=5)).put()
        nd.Pin(name='a0', xs=None).put(0,0,180)
        nd.Pin(name='b0', xs=None).put(0,0,0)
    return C
