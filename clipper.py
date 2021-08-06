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
# Grow and merge polygons.
# Uses the pyclipper module
#
# (c) 2017 Xaveer Leijtens, Ronald Broeke 2018 (c)
#
"""
Module with boolean operations on polygons and grow/shrink.

Documentation of pyclipper is minimum and that of the library is better. That
can be found at:
http://www.angusj.com/delphi/clipper/documentation
e.g. documentation/Docs/Units/ClipperLib/Types/PolyFillType.htm has useful
information.
"""
import nazca.cfg as cfg
from nazca.logging import logger

try:
    import pyclipper as pc

    cfg.PYCLIPPER = True
    __all__ = [
        "merge_polygons",
        "diff_polygons",
        "xor_polygons",
        "clip_polygons",
        "grow_polygons",
        "polygons_OR",
        "polygons_NOT",
        "polygons_XOR",
        "polygons_AND",
    ]
except Exception:  # as e:
    cfg.PYCLIPPER = False
    # print('Warning: Could not import pyclipper:', e)
if cfg.PYCLIPPER:
    st = pc.scale_to_clipper
    sf = pc.scale_from_clipper


clipper_check = False


def _has_pyclipper():
    """Check if pyclipper module is loaded.

    Returns:
        bool: True if pyclipper is loaded
    """
    global clipper_check
    if cfg.PYCLIPPER:
        return True
    else:
        if not clipper_check:
            logger.error(
                "Could not load module 'pyclipper'."
                " Skipping requested pyclipper functionality."
            )
            clipper_check = True
        return False


def signed_area(XY):
    """Calculate and return the signed area of a polygon: negative is
    counter_clockwise. Polygon should be closed (start and end point should be
    the same).

    This is a general purpose routine. It may deviate from the area of the
    finally drawn GDS polygon, because it does not take into account the GDS
    gridding.

    Args:
        XY (list): polygon, a list of (x, y) coordinates.

    Returns: (float)
    """
    area = 0
    ox, oy = XY[0]
    for x, y in XY[1:]:
        area += x * oy - y * ox
        ox, oy = x, y
    return area / 2


def _x_intersect(xy0, xy1, ya):
    """Calculate the x-intersect value for line through xy0, xy1 at height ya.
    This is NOT a general purpose routine.
    """
    x0, y0 = xy0
    x1, y1 = xy1
    if y0 == y1:
        return max(x0, x1)
    return x0 + (ya - y0) * (x1 - x0) / (y1 - y0)


def _leftmost(XY):
    """Find index of leftmost point in polygon point list.
    This is NOT a general purpose routine.

    Args:
        XY (list): polygon, a list of (x, y) coordinates.

    Returns: (float)
    """
    min = XY[0][0]
    res = 0
    for i, (x, _) in enumerate(XY[1:]):
        if x < min:
            min = x
            res = i + 1
    return res


def _poly_xmin(XY):
    """Return the x-coordinate of the leftmost point in a polygon.
    This is NOT a general purpose routine.

    Args:
        XY (list): polygon, a list of (x, y) coordinates.

    Returns: (float)
    """
    return XY[_leftmost(XY)][0]


def _subtract_polygon(XYo, XYi):
    """Subtract the inner polygon from the outer polygon. Polygons should
    originate from a clipper operation. The resulting polygon is compatible with
    GDS: a single connected polygon.

    This is NOT a general purpose routine.

    The first polygon should be counter-clockwise oriented, the second should be
    clockwise oriented. This is how they are returned from clipper.

    Returned is the combined polygon.

    Args:
        XYo (TODO): TODO
        XYi (TODO): TODO

    Returns: list
    """
    # The polygons should be closed.
    if XYo[0] != XYo[-1]:
        XYo.append(XYo[0])
    if XYi[0] != XYi[-1]:
        XYi.append(XYi[0])

    ndxi = _leftmost(XYi)
    xp, yp = XYi[ndxi]
    xo, yo = XYo[0]
    dmin = float("inf")
    # Select all polygons edges in the outer polygon that intersect with the
    # horizontal line through the leftmost point of the inner polygon.
    for i, (x, y) in enumerate(XYo[1:]):
        # print(f"{i} -> {i+1}: {(xo, yo)} -> {x, y}")
        if yo >= yp and y <= yp and (xo < xp or x < xp):
            # Line segment intersecting with horizontal.
            # Calculate distance to intersection point
            xi = _x_intersect((xo, yo), (x, y), yp)
            d = xp - xi
            if d > 0 and d < dmin:
                dmin = d
                ndxo = i + 1
                point = (xi, yp)
        xo, yo = x, y
    # Construct polygon
    poly = XYo[0:ndxo] + [point] + XYi[ndxi:] + XYi[0:ndxi + 1] + [point] + XYo[ndxo:]
    return poly


def _clipper2GDS(clipper_result):
    """Cleanup the clipper polygons for GDS use. The clipper_result consists of
    a list of polygons with clockwise or counter_clockwise orientation. This
    routing will convert this into a polygon which may contain holes in a way
    compatible with GDS: holes will have a 'tether' to the outside of the
    polygon.
    Note that the input is assumed to be in clipper coordinates. Severe
    rounding will occur if that is not the case.

    The function returns the result: a list of GDS compatible polygons,
    possibly with holes, in clipper coordinates.

    Args:
        clipper_result (list): list of input polygons (clipper output)

    Returns: (list)
    """
    # Determine polygon orientation.
    # Outer polygons always have Orientation True and holes have Orientation False.
    orientation = [pc.Orientation(p) for p in clipper_result]

    # Create a tether for the holes that are fully enclosed. The tether
    # points to the left. We start with the leftmost inner polygon that is
    # inside the specific outer polygon.
    outer = []
    inner = []
    for p, outside in zip(clipper_result, orientation):
        if outside:
            outer.append(p)
        else:
            inner.append(p)
    # Sort inner polygons by their leftmost x-coordinate
    inner = sorted(inner, key=_poly_xmin)

    todo = set(range(len(inner)))  # Keep track of holes that need to be done.
    result = []
    for j, po in enumerate(outer):
        for i, pi in enumerate(inner):  # Remove each hole from the outer polygon
            if i in todo and pc.PointInPolygon(pi[0], po):
                # Subtract with tether.
                po = _subtract_polygon(po, pi)
                todo.remove(i)
        result.append(po)

    if todo:
        print(
            f"Not all polygons used...{len(todo)+1} remaining!\n"
            "This is an error that can be caused by an accuracy that is "
            "set to a too large value."
        )
    return result


def merge_polygons(paths, accuracy=1e-8):
    """Merge polygons and return a list with the resulting polygon(s). This
    function returns the union of the polygons (logical OR of all paths).


    Args:
        paths (list): list of polygons. A polygon is list of coordinates (x, y).
        accuracy (float): accuracy [µm] of the location of the intermediate points.
            There is no good reason to deviate from the default (default = 1e-8 µm).

    Returns:
        list: list of polygons  [ [(x1, y1), (x2, y2), ...], [...], ... ]
    """
    if not _has_pyclipper():
        return paths
    sc = 1 / accuracy
    clipper = pc.Pyclipper()
    # Ensure correct orientation (and remember so we can re-reverse, and do not
    # change the original polygon).
    rev = []
    for p in paths:
        if not pc.Orientation(p):
            p.reverse()
            rev.append(p)
    clipper.AddPaths(st(paths, sc), pc.PT_SUBJECT, True)
    mp = clipper.Execute(pc.CT_UNION, pc.PFT_NONZERO, pc.PFT_NONZERO)
    mp = _clipper2GDS(mp)
    for p in rev:
        p.reverse()
    return sf(mp, sc)


def diff_polygons(paths_A, paths_B, accuracy=1e-8):
    """Subtract polygons and return a list with the resulting polygon(s). This
    function returns the subtraction of the polygons (paths_A NOT paths_B)

    Args:
        paths_A (list): list of polygons that each have a list of (x,y) coordinates.
        paths_B (list): list of polygons that each have a list of (x,y) coordinates.
        accuracy (float): accuracy [µm] of the location of the intermediate points.
        There is no good reason to deviate from the default (default = 1e-8 µm).

    Returns:
        list: List of polygons that result from subtraction: paths_A - paths_B.
    """
    if not _has_pyclipper():
        return paths_A + paths_B
    sc = 1 / accuracy
    rev = []
    for p in paths_A + paths_B:
        if not pc.Orientation(p):
            p.reverse()
            rev.append(p)
    clipper = pc.Pyclipper()
    clipper.AddPaths(st(paths_A, sc), pc.PT_SUBJECT, True)
    clipper.AddPaths(st(paths_B, sc), pc.PT_CLIP, True)
    sp = clipper.Execute(pc.CT_DIFFERENCE, pc.PFT_NONZERO, pc.PFT_NONZERO)
    for p in rev:
        p.reverse()
    sp = _clipper2GDS(sp)
    return sf(sp, sc)


def clip_polygons(paths_A, paths_B, accuracy=1e-8):
    """Clip polygons and return a list with the resulting polygon(s). This
    function returns the intersection of the polygons (paths_A AND paths_B).

    Args:
        paths_A (list): list of polygons that each have a list of (x,y) coordinates.
        paths_B (list): list of polygons that each have a list of (x,y) coordinates.
        accuracy (float): accuracy [µm] of the location of the intermediate points.
        There is no good reason to deviate from the default (default = 1e-8 µm).

    Returns:
        list: List of polygons that result from clipping: paths_A AND paths_B.
    """
    if not _has_pyclipper():
        return paths_A + paths_B
    sc = 1 / accuracy
    clipper = pc.Pyclipper()
    rev = []
    for p in paths_A + paths_B:
        if not pc.Orientation(p):
            p.reverse()
            rev.append(p)
    clipper.AddPaths(st(paths_A, sc), pc.PT_SUBJECT, True)
    clipper.AddPaths(st(paths_B, sc), pc.PT_CLIP, True)
    sp = clipper.Execute(pc.CT_INTERSECTION, pc.PFT_NONZERO, pc.PFT_NONZERO)
    for p in rev:
        p.reverse()
    sp = _clipper2GDS(sp)
    return sf(sp, sc)


def xor_polygons(paths_A, paths_B, accuracy=1e-8):
    """This function returns a list of polygons that result form the logical
    XOR operation: paths_A XOR paths_B.


    Args:
        paths (list): list of polygons. A polygon is list of coordinates (x, y).
        accuracy (float): accuracy [µm] of the location of the intermediate points.
            There is no good reason to deviate from the default (default = 1e-8 µm).

    Returns:
        list: list of polygons  [ [(x1, y1), (x2, y2), ...], [...], ... ]
    """
    if not _has_pyclipper():
        return paths_A + paths_B
    sc = 1 / accuracy
    clp = pc.Pyclipper()
    # Ensure correct orientation (and remember so we can re-reverse, and do not
    # change the original polygon).
    rev = []
    for p in paths_A + paths_B:
        if not pc.Orientation(p):
            p.reverse()
            rev.append(p)
    clp.AddPaths(st(paths_A, sc), pc.PT_SUBJECT, True)
    clp.AddPaths(st(paths_B, sc), pc.PT_CLIP, True)
    for p in rev:
        p.reverse()
    xor = clp.Execute(pc.CT_XOR, pc.PFT_NONZERO, pc.PFT_NONZERO)
    xor = _clipper2GDS(xor)
    return sf(xor, sc)


def grow_polygons(paths, grow=5, accuracy=0.1, jointype="round"):
    """Grow or shrink polygons and return the resulting structures.

    Args:
        paths (list): list of polygons that each have a list of (x,y) coordinates.
        accuracy (float): accuracy [µm] of the location of the intermediate points.
            The accuracy determines the grid on which the grown path is drawn.
        jointype: specifies the type of growing that is used.
            The jointype is one of 'round' (default), 'square' or 'miter'.

    Returns:
        list: list of points [(x1, y1), (x2, y2), ...]
    """
    if not _has_pyclipper():
        return paths
    if grow == 0:  # 0-grow would affect polygons because of the grid.
        return paths
    sc = 1 / accuracy
    # Offset is the term for growing.
    pco = pc.PyclipperOffset()
    # Join type
    jt = {
        "round": pc.JT_ROUND,
        "square": pc.JT_SQUARE,
        "miter": pc.JT_MITER,
    }
    if jointype not in jt:
        print("jointype '{}' unknown.".format(jointype))
        print("jointype should be one of 'round', 'square', 'miter'.")
        print("Using default ('round')")
        jointype = "round"
    # Path orientation matters
    rev = []
    for p in paths:
        if not pc.Orientation(p):
            p.reverse()
            rev.append(p)
    pco.AddPaths(st(paths, sc), jt[jointype], pc.ET_CLOSEDPOLYGON)
    for p in rev:
        p.reverse()
    pco = pco.Execute(grow * sc)
    pco = _clipper2GDS(pco)
    return sf(pco, sc)


polygons_AND = clip_polygons
polygons_OR = merge_polygons
polygons_NOT = diff_polygons
polygons_XOR = xor_polygons
