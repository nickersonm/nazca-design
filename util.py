# -----------------------------------------------------------------------
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
# -----------------------------------------------------------------------
#
# Utility routines
#
# (c) 2016-2018  Xaveer Leijtens, Ronald Broeke
#
import os
import re
from functools import partial
from collections import OrderedDict
from math import hypot, radians, cos, sin, sqrt
import hashlib

import nazca as nd
from nazca import gds_base as gbase
from nazca.logging import logger

__all__ = [#'parameters_to_string', 'string_to_parameters',
           'get_cell_annotation', 'get_cell_polyline', 'get_cell_polygon',
           'make_iter', 'md5', 'isnotebook']


#def parameters_to_string(param):
#    """Create a string from a parameter dictionary.
#
#    param (dict): (parameter_name, value)
#
#    Format:
#
#    "Parameters:
#    <parameter> = <value>
#    <parameter> = <value>
#    ..."
#
#    Returns:
#        str: parameters as a string
#    """
#    plist = ['Parameters:']
#    for key, value in param.items():
#        plist.append("{} = {}".format(key, value))
#    return '\n'.join(plist)
#
#
#def string_to_parameters(string):
#    """Convert a string to a parameter dictionary.
#
#    The returned parameter values are represented as type str.
#
#    Expected format of <string>:
#
#    "parameters:
#    <parameter> = <value>
#    <parameter> = <value>
#    ..."
#
#    Header 'parameters:' is case incensitive and spaces will be stripped.
#
#    Args:
#        string (str): parameters
#
#    Returns:
#        OrderedDict: {<parameter_name>: <parameter_value>}
#    """
#    lines = string.split('\n')
#    p = OrderedDict()
#    if (lines[0].lower() == 'parameters:'):
#        for line in lines[1:]:
#            param = line.split('=', 1)
#            if len(param) == 2:
#                p[param[0].strip()] = param[1].strip()
#            else:
#                logger.warning("string_to_parameter: Expected one keyword "
#                    "and one value,\nbut found instead: {}\n"
#                    "Provided string: {}".
#                    format(param, string))
#    else:
#        logger.error("string_to_parameter: "
#            "Expected header 'parameters:'\n"
#            "Provided string: {}".format(string))
#    return p


def get_cell_annotation(cell, convert=False):
    """Yield the <cell>'s annotations one by one.

    If convert is False then return XY as integers, as in GDSII file (nm).
    If convert is True then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): default convert=False

    Yields:
        int, (int, int) | (float, float): annotation layer, position, text
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.GDS_record.TEXT:
            lay, pos, text = e.annotation
            pos[0] *= conv
            pos[1] *= conv
            yield lay, pos, text


def get_cell_polyline(cell, convert=False):
    """Yield the <cell>'s polylines one by one.

    If convert is False then return XY as integers, as in GDSII file (nm).
    If convert is True then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool):  convert polyline's values to float (default = False)

    Yields:
        int, (int, int) | (float, float): layer, XY
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.GDS_record.PATH:
            lay, points = e.polyline
            XY = []
            for i in range(0, len(points), 2):
                XY.append((points[i] * conv, points[i + 1] * conv))
            yield lay, XY


def get_cell_polygon(cell, convert=False):
    """Yield the <cell>'s polygons one by one.

    If convert is False then return XY as integers, as in GDSII file (nm).
    If convert is True then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): convert polygon's values to float (default = False)

    Yields:
        int, (int, int) | (float, float): layer, XY
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.GDS_record.BOUNDARY:
            lay, points = e.polygon
            XY = []
            for i in range(0, len(points), 2):
                XY.append((points[i] * conv, points[i + 1] * conv))
            yield lay, XY


def make_iter(x):
    """Return x as tuple, if x is not a string and not iterable."""
    if x is None:
        return tuple()
    elif type(x) is str or not hasattr(x, "__iter__"):
        return x,
    else:
        return x


def md5(x, N=None):
    """Return first N characters of md5 hash of argument x.

    If no N is given the fill md5 is returned.

    The function hashes the (default) string representation of the object.

    Args:
        x (str): string to hash
        N (int): number of characters in the hash

    Returns:
        str: hash
    """
    if N is None:
        return hashlib.md5('{}'.format(x).encode()).hexdigest()
    else:
        return hashlib.md5('{}'.format(x).encode()).hexdigest()[:N]


def file2md5(filename, save=True, suffix='.md5', fullpath=False):
    """Create a md5sum of file <filename> and optionally save to file.

    Args:
        filename (str): name of the file to hash using md5
        save (bool): save the md5sum to file under name  <filename><suffix>. Default=True
        suffix (str): suffix when saving the md5sum
        fullpath (bool): save the ms5sum using the file path in <filename>. Default=False

    Returns:
        str: hash of file content
    """
    with open(filename, 'rb') as F:
        data = F.read()
        md5sum = hashlib.md5(data).hexdigest()
    if save:
        if fullpath:
            name = filename
        else:
            name = os.path.basename(filename)
        with open(filename+suffix, 'w') as Fout:
            Fout.write("{}  {}".format(md5sum, name))
    return md5sum


def isnotebook():
    """Check if code is run in a Jupyter notebook.

    Returns:
        bool: True if call is made from a notebook
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':  # Jupyter notebook or qtconsole?
            return True
        elif shell == 'TerminalInteractiveShell':  # Terminal running IPython?
            return False
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


def klayout2nazca(string):
    """Convert a copy-pasted path (x,y points) from Klayout into list of points.

    Args:
        string (str): string with points to convert

    Returns:
        list of (float, float): list of points (x, y)

    Example::

        klayout2nazca("1.0, 3.5, 2.0, 5.5")

        out: [(1.0, 3.5), (2.0, 5.5)]
    """
    return [[float(x) for x in line.split('\t')]
        for line in string.strip().split('\n')]


def _trapezoid(xy0, xy1, w0, w1):
    """Calculate the four coordinates of the trapezoid that is defined by the
    line segment through two points xy0 and xy1 with widths w0 and w1.

    Args:
        xy0 (float, float): point (x,y), segment start
        xy1 (float, float): point (x,y), segment end
        w0 (float): width of line segment start
        w1 (float): width of line segment end

    Returns:
        list of (float, float): list of 4 (float, float) point coordinates (x, y)
        tracing the outline of the trapezoid.
    """
    dx = xy1[1] - xy0[1]
    dy = xy0[0] - xy1[0]
    try:
        length = hypot(dx, dy)
        px = dx / 2.0 / length
        py = dy / 2.0 / length
    except: # ZeroDivisionError:
        logger.error("zero-length segment trapezoid line at ({},{}) in cell '{}'".format(xy0[0], xy0[1], nd.cfg.cells[-1].cell_name))
        px = py = 0
    dx0 = px * w0
    dy0 = py * w0
    dx1 = px * w1
    dy1 = py * w1
    return [
        (xy0[0] + dx0, xy0[1] + dy0),
        (xy1[0] + dx1, xy1[1] + dy1),
        (xy1[0] - dx1, xy1[1] - dy1),
        (xy0[0] - dx0, xy0[1] - dy0)
    ]


def _trapezoid_asym(xy0, xy1, w0a, w0b, w1a, w1b):
    """Calculate the four coordinates of the trapezoid that is defined by the
    line segment through two points xy0 and xy1 with widths w0 and w1.

    Args:
        xy0 (float, float): point (x,y), segment start
        xy1 (float, float): point (x,y), segment end
        w0a (float): top width of line segment start
        w0b (float): bottom width of line segment start
        w1a (float): top width of line segment end
        w1b (float): bottom width of line segment end

    Returns:
        list of (float, float): list of 4 (float, float) point coordinates (x, y)
        tracing the outline of the trapezoid.
    """
    dx = xy1[1] - xy0[1]
    dy = xy0[0] - xy1[0]
    try:
        length = hypot(dx, dy)
        px = dx / length
        py = dy / length
    except ZeroDivisionError:
        logger.error("zero-length segment polyline at ({},{})", xy0[0], xy0[1])
        px = py = 0
    dx0a = px * w0a
    dx0b = px * w0b
    dy0a = py * w0a
    dy0b = py * w0b
    dx1a = px * w1a
    dx1b = px * w1b
    dy1a = py * w1a
    dy1b = py * w1b
    return [
        (xy0[0] + dx0a, xy0[1] + dy0a),
        (xy1[0] + dx1a, xy1[1] + dy1a),
        (xy1[0] + dx1b, xy1[1] + dy1b),
        (xy0[0] + dx0b, xy0[1] + dy0b)
    ]


def _intersect(xy0, xy1, xy2, xy3):
    """Helper function to intersect two lines.

    Intersection point (xi, yi)
    of two lines that go through (x0, y0), (x1, y1) and (x2, y2), (x3, y3).
    When the lines are parallel, return the point in between points xy1 and
    xy2. This makes the order of the points important.
    If the points are very close to parallel, return the point in between.
    Otherwise, check if the intersection is in between the two points or
    outside. In the first case, return the intersection point. In the latter
    case, also return the point in between.

    Args:
        xy0, xy1, xy2, xy3: list of four points (x, y)

    Returns:
        (float, float): intersection point (xi, yi)
    """
    x0, y0 = xy0
    x1, y1 = xy1
    x2, y2 = xy2
    x3, y3 = xy3
    D0 = (x3 - x2) * (y1 - y0) - (x1 - x0) * (y3 - y2)
    Dx = (x3 - x2) * (y2 - y0) - (x2 - x0) * (y3 - y2)
    if abs(D0) < 1e-12:  # Lines AB & CD (almost) parallel
        return ((x1 + x2) / 2, (y1 + y2) / 2)  # Point in between
    # Calculate intersection point
    xi = x0 + Dx / D0 * (x1 - x0)
    yi = y0 + Dx / D0 * (y1 - y0)
    # Check if intersection is in between B & C
    lsqr = (x2 - x1)**2 + (y2 - y1)**2
    s = (xi - x1) * (x2 - x1) + (yi - y1) * (y2 - y1)
    if 0 < s < lsqr:  # Intersection in between B & C
        return (xi, yi)
    return ((x1 + x2) / 2, (y1 + y2) / 2)  # Point in between


def polyline_length(xy):
    """Return the lenght of the polyline, which is the sum of the line
    segments in the polyline.

    Args:
        xy (list): list of (x,y) points that hold the polygon.

    Returns:
        length (float): the length of the polyline.
    """
    length = 0
    for i in range(1, len(xy)):
        length += hypot((xy[i][0] - xy[i - 1][0]), (xy[i][1] - xy[i - 1][1]))
    return length


def polyline2polygons(
    xy,
    width=2.0,
    width2=None,
    parabolic=True,
    miter=0.5,
    anglei=None,
    angleo=None,
    split=True,
):
    """Return a list of polygons that contain the outline points of a polyline
    with given width.

    Since we have to specify the outline of two or more segments that make
    an angle, we have to know what to do with the gap between those
    segments at the outside of the corners. In order to determine which
    points are on the outside of the corners we use the following algorithm:

    Given a line segment between P0 (x0, y0) and P1 (x1, y1), another point P
    (x,y) has the following relationship to the line segment. Compute
    (y - y0) (x1 - x0) - (x - x0) (y1 - y0). If it is less than 0 then P is
    to the right of the line segment, if greater than 0 it is to the left,
    if equal to 0 then it lies on the line segment.

    The routine fills an array from the top with the anticlockwise points
    and from the bottom with the clockwise points.

    Args:
        xy (list): list of (x,y) points that holds the polygon
        width (float | list | function): width of the polyline (default 2), or a
            parameterized function w(t) returning the width for the independent
            variable t wich runs from 0 to 1 from start of the polyline to the
            end, proportional to the length of the polyline segments.
        width2 (float): if width is a number, and if width2 is not None, they
            are interpreted as the start width and end width2 with a parabolic
            change of width vs length.
        parabolic (bool): if begin and end widths are specified as numbers, use
            parabolic tapering, or linear tapering if False (default True).
        miter (float): maximum fraction of the width before an extra point
            is added in outside corners (default=0.5).
        anglei (float): force input side angle in degrees (default=None,
            perpendicular to polyline). Useful if the last polyline segment is
            slightly short of the final desired angle.
        angleo (float): force output side angle  in degrees (default=None,
            perpendicular to polyline) Useful if the last polyline segment is
            slightly short of the final desired angle.
        split (bool): if the number of points becomes too large, split the polygon.

    Returns:
        list of polygons: the polygons are each lists of coordinates (float, float).
    """
    # Each iteration can add 3 points and 2 for the closing polygon
    nmax = nd.cfg.maxpolygonpoints - 5
    n = len(xy)
    if n < 2:
        raise ValueError(
            "Polyline2polygons: need at least 2 points for polyline.")

    # start angle correction
    # TODO: update xy for anglei and angleo at the start.
    N = n
    # 1/3 creates effectively 3 near-equal line segments inside the
    # original 2 of connecting elements: ---|--- -> --.-|-.--
    # with the exact angle in the middle segment and | denoting the element transition
    factor = 1.0/3.0
    if anglei is not None:
        l0 = hypot(xy[1][0] - xy[0][0], xy[1][1] - xy[0][1])
        xi = xy[0][0] + (l0/3.0) * cos(radians(anglei))
        yi = xy[0][1] + (l0/3.0) * sin(radians(anglei))
    if angleo is not None:
        l0 = hypot(xy[n-2][0] - xy[n-1][0], xy[n-2][1] - xy[n-1][1])
        xo = xy[n-1][0] + factor * l0 * cos(radians(angleo+180))
        yo = xy[n-1][1] + factor * l0 * sin(radians(angleo+180))

    xy2 = []
    for i, point in enumerate(xy):
        if i == N-1 and angleo is not None:
            n += 1
            xy2.append((xo, yo))
        xy2.append(point)
        #print(point)
        if i == 0 and anglei is not None:
            n += 1
            xy2.append((xi, yi))
    xy = xy2

    if isinstance(width, list):
        width_new = []
        for i, W in enumerate(width):
            if i == N-1 and angleo is not None:
                width_new.append((width[N-1] + factor * (width[N-2]-width[N-1]) ))
            width_new.append(W)
            #print(i, W)
            if i == 0 and anglei is not None:
                width_new.append((width[0] + factor * (width[1]-width[0]) ))
        width = width_new
    # end angle correction


    # TODO: check discretization of result against resolution.

    # Determine the width at each point of the polyline from the
    # caller-supplied width function (or fixed width)
    # First derive the independent variable from the length of each section.
    length = [0]
    ltot = 0
    for ndx in range(1, n):
        ltot += hypot(xy[ndx][0] - xy[ndx - 1][0], xy[ndx][1] - xy[ndx - 1][1])
        length.append(ltot)
    t = [lrun / ltot for lrun in length]  # Normalize to [0,1]
    # Note: if there is a curvature in the polyline, "length" is a lower boundary.

    def linvarwidth(w1, w2, t):
        return t * (w2 - w1) + w1

    def parabolicvarwidth(w1, w2, t):
        return sqrt(w1 * w1 * (1 - t) + w2 * w2 * t)

    if isinstance(width, (int, float)) and isinstance(width2, (int, float)):
        if parabolic:
            # parabolic width from w=width to w=width2.
            width = partial(parabolicvarwidth, width, width2)
        else:
            # linear width from w=width to w=width2.
            width = partial(linvarwidth, width, width2)
        # width(t) is now a function and rest will be caught in "if" below.

    if isinstance(width, (int, float)):
        # Constant width w=width
        w = [width for k in range(n)]
        dsqrmax = [(miter * width)**2 for k in range(n)]
    elif isinstance(width, list):
        # List of width: should have the same length as xy
        if len(width) != n:
            raise ValueError('Polyline2polygons: length of xy needs '
                'to match length of width.')
        w = width
        dsqrmax = [(miter * w)**2 for w in width]
    elif callable(width):
        # w=width(t) is a function of t.
        w = [width(x) for x in t]
        # the fraction of the width of the line segments that is used to
        # determine if a single point is sufficient to describe the outline, or
        # that two points are needed (miter limit).
        dsqrmax = [(miter * width(x))**2 for x in t]
    else:
        raise ValueError("Polyline2polygons: don't know what to do with "
             "this width parameter: {}.".format(width))

    # Start with the first two points:
    tr1 = _trapezoid(xy[0], xy[1], w[0], w[1])
    xyt = [tr1[0]]  # Top coordinates of polygon
    xyb = [tr1[3]]  # Bottom coordinates of polygon

    XY = []  # List of polygons to be returned
    # loop over the points in the polyline:
    for i in range(1, n - 1):
        tr0 = tr1  # Current and next trapezoid
        # Get corner points for next segment.
        tr1 = _trapezoid(xy[i], xy[i + 1], w[i], w[i + 1])
        # left or right turn
        lrt = (xy[i + 1][1] - xy[i - 1][1]) * (xy[i][0] - xy[i - 1][0]) -\
              (xy[i + 1][0] - xy[i - 1][0]) * (xy[i][1] - xy[i - 1][1])
        # Distance (squared) between the two points at the kink
        # (Top and bottom have the same distance)
        dsqr = (tr1[0][0] - tr0[1][0])**2 + (tr1[0][1] - tr0[1][1])**2
        if lrt < 0:  # Left turn
            # Inside corner: always intersect.
            xyt.append(_intersect(tr0[0], tr0[1], tr1[0], tr1[1]))
            # Outside corner: use two points, unless these points are close.
            if dsqr < dsqrmax[i]:
                xyb.append(_intersect(tr0[3], tr0[2], tr1[3], tr1[2]))
            else:
                xyb.append(tr0[2])
                xyb.append(tr1[3])
        elif lrt > 0:  # Right turn
            # Inside corner: always intersect.
            xyb.append(_intersect(tr0[3], tr0[2], tr1[3], tr1[2]))
            # Outside corner: use two points, unless these points are close.
            if dsqr < dsqrmax[i]:
                xyt.append(_intersect(tr0[0], tr0[1], tr1[0], tr1[1]))
            else:
                xyt.append(tr0[1])
                xyt.append(tr1[0])
        else:
            xyt.append(tr0[1])
            xyb.append(tr0[2])
        # Start a new polygon if the number of points becomes too large.
        if split and len(xyt) + len(xyb) >= nmax:
            XY.append(xyt + list(reversed(xyb)))
            xyt = [xyt[-1]]
            xyb = [xyb[-1]]

    # Last two points:
    xyt.append(tr1[1])
    xyb.append(tr1[2])
    XY.append(xyt + list(reversed(xyb)))
    return XY


def polyline2polygon(
    xy,
    width=2.0,
    width2=None,
    parabolic=True,
    miter=0.5,
    anglei=None,
    angleo=None,
):
    """Return a polygon that contains the outline points of a polyline with given width.

    Args:
        xy (list): list of (x,y) points that holds the polygon
        width (float | list | function): width of the polyline (default 2), or a
            parameterized function w(t) returning the width for the independent
            variable t wich runs from 0 to 1 from start of the polyline to the
            end, proportional to the length of the polyline segments.
        width2 (float): if width is a number, and if width2 is not None, they
            are interpreted as the start width and end width2 with a parabolic
            change of width vs length.
        parabolic (bool): if begin and end widths are specified as numbers, use
            parabolic tapering, or linear tapering if False (default True).
        miter (float): maximum fraction of the width before an extra point
            is added in outside corners (default=0.5).
        anglei (float): force input side angle in degrees (default=None,
            perpendicular to polyline). Useful if the last polyline segment is
            slightly short of the final desired angle.
        angleo (float): force output side angle  in degrees (default=None,
            perpendicular to polyline) Useful if the last polyline segment is
            slightly short of the final desired angle.

    Returns:
        list of (float, float): the polygon.
    """
    return polyline2polygons(xy=xy, width=width, width2=width2,
            parabolic=parabolic, miter=miter, anglei=anglei, angleo=angleo,
            split=False)[0]


def polyline2edge(
    xy,
    width1,
    width2=None,
    grow=None,
    parabolic=True,
    miter=0.5,
    anglei=None,
    angleo=None,
    line=False,
):
    """Return a polygon that contains the outline points of a polyline with
    given width.

    This method is based on edges a1 * width + b1 and a2 * width + b2
    and it does not assume symmetry along the spine.

    Since we have to specify the outline of two or more segments that make
    an angle, we have to know what to do with the gap between those
    segments at the outside of the corners. In order to determine which
    points are on the outside of the corners we use the following algorithm:
    Given a line segment between P0 (x0, y0) and P1 (x1, y1), another point P
    (x,y) has the following relationship to the line segment. Compute
    (y - y0) (x1 - x0) - (x - x0) (y1 - y0). If it is less than 0 then P is
    to the right of the line segment, if greater than 0 it is to the left,
    if equal to 0 then it lies on the line segment.
    The routine fills an array from the top with the anticlockwise points
    and from the bottom with the clockwise points.

    Args:
        xy (list): list of (x,y) points that hold the polygon
        width (float | list | function): width of the polyline (default 2), or a
            parameterized function w(t) returning the width for the independent
            variable t wich runs from 0 to 1 from start of the polyline to the
            end, proportional to the length of the polyline segments.
        width2 (float): if width is a number, and if width2 is not None, they
            are interpreted as the start width and end width2 with a parabolic
            change of width vs length.
        grow (tuple): leftedge and rightedge as in a*width+b as ((a1, b1), (a2, b2))
        parabolic (bool): if begin and end widths are specified as numbers, use
            parabolic tapering, or linear tapering if False (default True).
        miter (float): maximum fraction of the width before an extra point
            is added in outside corners (default 0.5).
        line (bool): Return only a spine (line) if True (default=False).
            The spine will be in the center of the edges as set by xy, width1 and widht2.

    Returns:
        list of (float, float): the polygon
    """
    # TODO: check discretization of result against resolution.
    n = len(xy)
    if n < 2:
        raise ValueError(
            "Polyline2polygon: need at least 2 points for polyline.")

    (a1, b1), (a2, b2), c1, c2 = grow
    (a1, b1), (a2, b2) = (-a1, -b1), (-a2, -b2)

    if line:
        position = 0.5*(a1+a2) * width1 + 0.5*(b1+b2)
        b1 = b2 = position
        a1 = a2 = 0
        width1 = 0

    # start angle correction
    # TODO: update xy for anglei and angleo at the start.
    N = n
    # 1/3 creates effectively 3 near-equal line segments inside the
    # original 2 of connecting elements: ---|--- -> --.-|-.--
    # with the exact angle in the middle segment and | denoting the element transition
    factor = 1.0/3.0
    if anglei is not None:
        l0 = hypot(xy[1][0] - xy[0][0], xy[1][1] - xy[0][1])
        xi = xy[0][0] + (l0/3.0) * cos(radians(anglei))
        yi = xy[0][1] + (l0/3.0) * sin(radians(anglei))
    if angleo is not None:
        l0 = hypot(xy[n-2][0] - xy[n-1][0], xy[n-2][1] - xy[n-1][1])
        xo = xy[n-1][0] + factor * l0 * cos(radians(angleo+180))
        yo = xy[n-1][1] + factor * l0 * sin(radians(angleo+180))

    xy2 = []
    for i, point in enumerate(xy):
        if i == N-1 and angleo is not None:
            n += 1
            xy2.append((xo, yo))
        xy2.append(point)
        #print(point)
        if i == 0 and anglei is not None:
            n += 1
            xy2.append((xi, yi))
    xy = xy2

    if isinstance(width1, list):
        width_new = []
        for i, W in enumerate(width1):
            if i == N-1 and angleo is not None:
                width_new.append((width1[N-1] + factor * (width1[N-2]-width1[N-1]) ))
            width_new.append(W)
            #print(i, W)
            if i == 0 and anglei is not None:
                width_new.append((width1[0] + factor * (width1[1]-width1[0]) ))
        width1 = width_new
    # end angle correction

    # Determine the width at each point of the polyline from the
    # caller-supplied width function (or fixed width)
    # First derive the independent variable from the length of each section.
    length = [0]
    ltot = 0
    for ndx in range(1, n):
        ltot += hypot(xy[ndx][0] - xy[ndx - 1][0], xy[ndx][1] - xy[ndx - 1][1])
        length.append(ltot)
    t = [lrun / ltot for lrun in length]  # Normalize to [0,1]

    def linvarwidth(w1, w2, t):
        return t * (w2 - w1) + w1

    def parabolicvarwidth(w1, w2, t):
        return sqrt(w1 * w1 * (1 - t) + w2 * w2 * t)

    if isinstance(width1, (int, float)) and isinstance(width2, (int, float)):
        if parabolic:
            # parabolic width from w=width to w=width2.
            width1 = partial(parabolicvarwidth, width1, width2)
        else:
            # linear width from w=width to w=width2.
            width1 = partial(linvarwidth, width1, width2)
        # width(t) is now a function and rest will be caught in "if" below.

    if isinstance(width1, (int, float)):
        # Constant width w=width
        wa = [width1*a1+b1 for k in range(n)]
        wb = [width1*a2+b2 for k in range(n)]
        dsqrmax = [(miter * width1)**2 for k in range(n)]
    elif isinstance(width1, list):
        # List of width: should have the same length as xy
        if len(width1) != n:
            raise ValueError('Polyline2polygon: length of xy needs '
                'to match length of width.')
        wa, wb, dsqrmax = [], [], []
        for w in width1:
            wa.append(w*a1+b1)
            wb.append(w*a2+b2)
            dsqrmax.append((miter * w)**2)
    elif callable(width1):
        # w=width(t) is a function of t.
        wa = [width1(x)*a1+b1 for x in t]
        wb = [width1(x)*a2+b2 for x in t]
        # the fraction of the width of the line segments that is used to
        # determine if a single point is sufficient to describe the outline, or
        # that two points are needed (miter limit).
        dsqrmax = [(miter * width1(x))**2 for x in t]
    else:
        raise ValueError("Polyline2polygon: don't know what to do with "
             "this width parameter: {}.".format(width1))

    # Start with the first two points
    tr1 = _trapezoid_asym(xy[0], xy[1], wa[0], wb[0], wa[1], wb[1])
    xyt = [tr1[0]]  # Top coordinates of polygon
    xyb = [tr1[3]]  # Bottom coordinates of polygon

    # loop over the points in the polyline:
    for i in range(1, n - 1):
        tr0 = tr1  # Current and next trapezoid
        # Get corner points for next segment.
        tr1 = _trapezoid_asym(xy[i], xy[i + 1], wa[i], wb[i], wa[i + 1], wb[i + 1])
        # left or right turn
        lrt = (
            (xy[i + 1][1] - xy[i - 1][1]) * (xy[i][0] - xy[i - 1][0]) -
            (xy[i + 1][0] - xy[i - 1][0]) * (xy[i][1] - xy[i - 1][1])
        )
        # Distance (squared) between the two points at the kink
        # (Top and bottom have the same distance)
        dsqr = (tr1[0][0] - tr0[1][0])**2 + (tr1[0][1] - tr0[1][1])**2

        if not line and lrt > 0:  # Left turn
            # Inside corner: always intersect.
            xyt.append(_intersect(tr0[0], tr0[1], tr1[0], tr1[1]))
            # Outside corner: use two points, unless these points are close.
            if dsqr < dsqrmax[i]:
                xyb.append(_intersect(tr0[3], tr0[2], tr1[3], tr1[2]))
            else:
                xyb.append(tr0[2])
                xyb.append(tr1[3])
        elif not line and lrt < 0:  # Right turn
            # Inside corner: always intersect.
            xyb.append(_intersect(tr0[3], tr0[2], tr1[3], tr1[2]))
            # Outside corner: use two points, unless these points are close.
            if dsqr < dsqrmax[i]:
                xyt.append(_intersect(tr0[0], tr0[1], tr1[0], tr1[1]))
            else:
                xyt.append(tr0[1])
                xyt.append(tr1[0])
        else:
            xyt.append(tr0[1])
            xyb.append(tr0[2])
    # Last two points.
    xyt.append(tr1[1])
    xyb.append(tr1[2])
    if line:
        return xyb
    else:
        return xyt + list(reversed(xyb))

def viper(x, y, w, N=200, anglei=None, angleo=None):
    """Parametric curve in t.

    t on interval [0, 1]

    Args:
        x (function): x-coordinate as function of t: x(t)
        y (function): y-coordinate as function of t: y(t)
        w (function): width as function of t: w(t)
        N (int): number of points
        anglei (float): force input side angle in degrees (default=None, perpendicular to polyline).
            Useful if the last polyline segment is slightly short of the final desired angle.
        angleo (float): force output side angle in degrees (default=None, perpendicular to polyline)
            Useful if the last polyline segment is slightly short of the final desired angle.

    Returns:
        list of (float, float): viper polygon
    """
    xy = []
    width = []
    for i in range(N):
        t = i / (N - 1)
        xy.append((x(t), y(t)))
        width.append(w(t))
    return polyline2polygon(xy, width=width, anglei=anglei, angleo=angleo)


def transform_polygon(points, dx=0.0, dy=0.0, da=0.0, scale=1.0, flipx=False,
        flipy=False, x=0.0, y=0.0):
    """Transform a polygon by translation, rotation, scaling and/or flipping.

    The transformation first applies (dx, dy) to reposition the origin.
    Subsequently, the scale, rotate and flips are applied, where order does not
    matter. Finally, a (x, y) translation is performed.

    Args:
        polygon (list of (float, float)): points (x, y)
        dx (float): x translation in um (default = 0.0)
        dy (float): y translation in um (default = 0.0)
        da (float): a translation in deg (default = 0.0)
        scale (float): scaling factor (default = 1.0)
        flipx (bool): flip x coordinate x -> -x (default = False)
        flipy (bool): flip y coordinate y -> -y (default = False)
        x (float): final x translation (after other transformations)
        y (float): final y translation (after other transformations)

    Returns:
        (list of (float, float)): transformed polygon points
    """
    fu, fv = 1, 1
    if flipx:
        fu = -1
    if flipy:
        fv = -1
    a = radians(da)
    xy = []
    for u, v in points:
        u = (u + dx) * fu
        v = (v + dy) * fv
        xy.append((x + scale * (cos(a) * u - sin(a) * v),
            y + scale * (sin(a) * u + cos(a) * v)))
    return xy


def read_and_filter_ascii(filename):
    """Read ascii layout export and delete some lines for diff: lines that
    contain a date/time ('bgnstr' and 'bgnlib') and the libname line, which can
    differ in length and in content.

    Args:
        filename (str): ascii layout file to read and filter.

    Returns:
        str: ascii layout in <filename> output minus problematic lines
    """
    with open(filename, 'r') as fref:
        ref = fref.readlines()
    lineiter = iter(ref)
    file = []
    for line in lineiter:
        if (
            line.startswith("bgnstr")
            or line.startswith("bgnlib")
            or line.startswith("libname")
        ):
            next(lineiter)
            continue
        file.append(line)
    return ''.join(file)


def instantiate_full_nazca_tree():
    """Set all internal Nazca cell instantiation options to True.

    A call to this function will instantiate the following cells:
    - pin symbols
    - stubs
    - mask elements
    These cells are normally not instantiated in GDS export to obtain a
    clean cell hierarchy.

    Returns:
        None
    """
    nd.cfg.instantiate_pin = True
    nd.cfg.instantiate_stub = True
    nd.cfg.instantiate_mask_element = True


def multisub(submapping, subject):
    """Simultaneously perform substitutions of a list of string instances in a subject string.

    Avoids sequential replacements, for example:
    subs = [['A', 'B'], ['B', 'CA']]
    subject = 'AB'
    returns 'BCA' (not 'CACA' or 'BCB' as in sequential replacements)

    Args:
        submapping (list of (str, str)): string replacement mapping to substite (old, new)
        subject (str): string to be subject to the substitution

    Returns:
        str: str with replacements applied.
    """
    if not submapping:
        return subject
    pattern = '|'.join('({:s})'.format(re.escape(old)) for old, new in submapping)
    substs = [new for old, new in submapping]
    replace = lambda m: substs[m.lastindex - 1]
    return re.sub(pattern, replace, subject)


def Tp_fan(N=2):
    """Tempate function for method fan with a closure on N
   
    This allows to define the fan() function with a preset N.

    Returns:
        function: fan for N.
    """
    def fan(y1=0, y2=1, i=0, N=N):
        """Auxilary fanout normalization function to scale functions for an N-ribbon.

        Note that there is always a zero point in the interval i in [0, N-1],
        Hence for both y1 and y2 positive the function will be a v-shape in i.

        Args:
            y1 (float): value at i=0
            y2 (float): value at i=N-1
            N (int): number of point
            i (int): counter between 0 and N (not including N)

        Returns:
            float: intermediate value for at i.
        """
        if y1 * y2 < 0:
            return y1 + ((y2 - y1) * i / (N - 1))
        else:
            return abs(y1 - abs((y1 + y2) * i / (N - 1)))
    return fan


def fan(y1, y2, i, N):
     """Auxilary fanout normalization function to scale functions for a N-ribbon.

     Note that there is always a zero point in the interval i in [0, N-1],
     Hence for both y1 and y2 positive the function will be a v-shape in i.

     Args:
         y1 (float): value at i=0
         y2 (float): value at i=N-1
         N (int): number of point
         i (int): counter between 0 and N (not including N)

     Returns:
         float: intermediate value for at i.
     """
     if y1 * y2 < 0:
         return y1 + ((y2 - y1) * i / (N - 1))
     else:
         return abs(y1 - abs((y1 + y2) * i / (N - 1)))


if __name__ == "__main__":
    fanN = Tp_fan(N=8)
    print(fanN(1, -1, 0))
    print(fanN(1, -1, 7))
    print(fan(1, -1, 7, 8))
