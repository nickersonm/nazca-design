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
# Split polygons if they contain more than a maximum number of points.
#
# (c) 2021 Xaveer Leijtens
#
"""
Polygon splitter.

This module has a routine to limit the maximum number of polygon points, by splitting a
polygon that contains too many points into smaller polygons. The shape of the combined
polygons stays the same.

Typical use is indirect, by setting a maximum number in the cfg module, e.g.:
    nd.cfg.maxpolygonpoints = 600

According to the GDSII standard, this should be 600, but there is no direct reason for
that. Older software of equipment may not be able to handle larger polygons. However, it
is assumed to be safe (and more efficient) to set this to a larger limit.

The maxiumum value for Nazca is 8190 and the current value can be found from the cfg
module (nd.cfg.maxpolygonpoints).
"""
import nazca as nd
from nazca.util import boundingbox
from operator import itemgetter
from itertools import groupby
from numpy import ndarray, sqrt

__all__ = ["limit_polygon_points"]


class Point:
    """Simple point class to keep the implementation close to the algorithm"""

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def cross(self, p):
        """2-dimensional cross product (scalar result)"""
        return self.x * p.y - self.y * p.x

    def dot(self, p):
        """2-dimensional cross product (scalar result)"""
        return self.x * p.x + self.y * p.y

    def __str__(self):
        return f"({self.x}, {self.y})"

    def __sub__(self, p):
        return Point(self.x - p.x, self.y - p.y)

    def __add__(self, p):
        return Point(self.x + p.x, self.y + p.y)

    def __mul__(self, t):
        return Point(self.x * t, self.y * t)


# https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
def intersect(a, b, c, d, firstisline=False):
    """Intersect AB with CD

    Args:
        a (list): [x, y] of start of first line segment
        b (list): [x, y] of end of first line segment
        c (list): [x, y] of start of second line segment
        d (list): [x, y] of end of second line segment
        firstisline (boolean): treat first segment as full line (default False)

    Return:
        list of length 0, 1 or 2 that contains the intersection point(s)
    """
    # The points will eventually be on a grid of size gds_db_unit, so a factor 100
    # smaller is a safe value.
    eps = nd.gds_base.gds_db_unit / 100
    # Define vectors p, p+r and q, q+s
    p = Point(a[0], a[1])
    r = Point(b[0] - a[0], b[1] - a[1])
    q = Point(c[0], c[1])
    s = Point(d[0] - c[0], d[1] - c[1])

    r_cross_s = r.cross(s)
    if abs(r_cross_s) < eps:
        if abs((q - p).cross(r)) < eps:
            # Case 1: collinear
            t0 = (q - p).dot(r) / (r.dot(r))
            t1 = t0 + s.dot(r) / (r.dot(r))
            result = []
            kind = 0
            if firstisline or -eps < t0 < 1 + eps:
                result.append(list(c))
                kind = 1
            if firstisline or -eps < t1 < 1 + eps:
                result.append(list(d))
                kind = 2 + kind
            # kind = 0 (no), 1 (1st point), 2 (2nd point), 3 (1st and 2nd point)
            return kind, result
        else:
            # Case 2: parallel, non-intersecting
            return 0, []
    else:
        t = (q - p).cross(s) / r_cross_s
        u = (q - p).cross(r) / r_cross_s
        if (firstisline or -eps < t < 1 + eps) and -eps <= u <= 1 + eps:
            # Case 3: intersecting
            pi = p + r * t
            if abs(u) < eps:
                return 4, [[pi.x, pi.y]]  # Starts at line (segment)
            elif abs(u - 1) < eps:
                return 5, [[pi.x, pi.y]]  # Ends at line (segment)
            else:
                return 6, [[pi.x, pi.y]]  # Intersects line (segment)
        else:
            # Case 4: not parallel, not intersecting
            return 0, []


class Node(object):
    """Simple node to use in doubly-linked list."""

    def __init__(self, data):
        self.data = data
        self.next = None
        self.prev = None
        self.visited = None


def polygon_split(s0, s1, p):
    """Split polygon p along the line segment s0-s1.
    This is implementation is an adaptation of the algorithm outlined by David Geier:
    https://geidav.wordpress.com/2015/03/21/splitting-an-arbitrary-polygon-by-a-line/
    But his algorithm doesn't properly take into account coinciding points on
    the split line. This one does by removing consecutive polygon points that are on the
    split line and only leaves sections that eventually cross the split line (LOR and
    ROL). Only these are then used to split the polygon.

    Args:
        s0 (float, float): (x, y) coordinate of segment start.
        s1 (float, float): (x, y) coordinate of segment end.
        p (list): polygon, a list (or ndarray) of (x, y) coordinates.

    Returns:
        list of polygons that are each a list of (x, y) coordinates.
    """
    eps = 1e-5
    ds = sqrt((s0[0] - s1[0]) ** 2 + (s0[1] - s1[1]) ** 2)

    if p[0] == p[-1]:  # Open polygon if closed.
        p.pop()
    if nd.util.signed_area(p) > 0:  # Needs to be counter clockwise
        p.reverse()

    index_on_line = []  # index of points on intersection line.
    # Since we will modify the polygon as we go, we cannot just loop over the
    # elements, but use indices in the while loop:
    lp = len(p)
    i = 0
    while i < lp:
        # Candidate segments
        kind, points = intersect(s0, s1, p[i], p[(i + 1) % lp], True)
        if kind in {1, 3, 4}:  # Add index of first point p[i]
            index_on_line.append(i)
        elif kind == 6:  # Point is new: insert. It will be added on next iteration.
            p.insert(i + 1, points[0])
            lp += 1
        i += 1
    l = len(index_on_line)
    i = 0
    while i < l:  # index_on_line is typically a small list.
        # Remove consecutive points on the line.
        if (index_on_line[i] + 1) % lp == index_on_line[(i + 1) % l]:
            del index_on_line[i]
            l -= 1
        else:
            i += 1
    # Now we'll use a linked list to represent the polygon. This is needed to
    # later be able to split off separate parts.
    node = [Node(i) for i in range(lp)]
    # Link them
    for i in range(lp):
        node[i - 1].next = node[i]
        node[(i + 1) % lp].prev = node[i]
        node[i].prev = node[i - 1]
        node[i].next = node[(i + 1) % lp]

    # Distance from start of the dividing line segment is later used for sorting.
    def signed_distance(i):
        # This is using the 2D (scalar) cross product.
        d = (p[i][0] - s0[0]) * (s1[0] - s0[0]) + (p[i][1] - s0[1]) * (s1[1] - s0[1])
        # s0 and s1 are ds apart. Scale such that the proper digits are significant.
        return d / ds

    def lor(node):
        """Return L, O, R"""
        i = node.data
        # This is using the dot product: left, right or on the line.
        d = (p[i][0] - s0[0]) * (s1[1] - s0[1]) - (p[i][1] - s0[1]) * (s1[0] - s0[0])
        if d < -eps:
            return "L"
        if d > eps:
            return "R"
        return "O"

    on_line = []  # Stores all ROL and LOR points on intersection line.
    # We only want LOR, ROL. There can be successive points on the line. Those we will
    # skipped until we find either "L" or "R".
    for ndx in index_on_line:
        nod = node[ndx]
        lr = lor(nod.prev)
        while lr == "O":  # Previous points: repeat until different from "O"
            nod = nod.prev
            lr = lor(nod.prev)
        LOR = lr + "O"
        lr = lor(node[ndx].next)
        while lr == "O":  # Next points: repeat until different from "O"
            nod = nod.next
            lr = lor(nod.next)
        LOR = LOR + lr
        if LOR in {"LOR", "ROL"}:
            # Replace the index stored in index_on_line with a tuple:
            # - node
            # - signed distance
            # - LOR
            on_line.append([node[ndx], signed_distance(ndx), LOR])

    on_line.sort(key=itemgetter(1))
    # Only LOR and ROL in on_line and they should be in pairs, LOR-ROL, starting with
    # LOR, because of the chosen orientation.
    lorrol = ["LOR", "ROL"]
    for i in range(len(on_line) - 1):
        node, d, LOR = on_line[i]
        if LOR != lorrol[i % 2]:
            # Not in the right order, d must be the same and swap them
            # (allow a change in the last significant digit).
            if abs(d - on_line[i + 1][1]) < eps:
                on_line[i], on_line[i + 1] = on_line[i + 1], on_line[i]
            else:
                raise ValueError("Self-intersecting polygon")

    # Now guaranteed to be correct (LOR-ROL) * N.
    for a, b in zip(on_line[::2], on_line[1::2]):  # Connection-split magic:
        src_node, dst_node = a[0], b[0]
        # Bridge source and destination
        new_src = Node(src_node.data)  # New source node
        new_dst = Node(dst_node.data)  # New destination node
        new_src.next = dst_node
        new_dst.next = src_node
        new_dst.prev = dst_node.prev
        new_dst.prev.next = new_dst
        dst_node.prev = new_src
        new_src.prev = src_node.prev
        src_node.prev.next = new_src
        src_node.prev = new_dst

    pols = []  # Resulting polygons
    for node, _, _ in on_line:
        pol = []
        points = []
        while not node.visited:
            points.append(node.data)
            pol.append(p[node.data])
            node.visited = True
            node = node.next
        if pol:
            pols.append(pol)
    return pols


def limit_polygon_points(p, nmax=None):
    """Limit the maximum number in a polygon by splitting it into multiple
    polygons if needed. The returned (open) polygon(s) will have at most nmax-1
    points. The '-1' is for the closing point that will be added when written
    to GDS. If nmax is not specified, then he maximum number comes from the cfg
    module: maxpolygonpoints.

    Args:
        p (list): polygon with list of (x, y) values.
        nmax (int): maximum number of points in resulting polygon(s).

    Returns:
        list of polygons, that each are a list of (x, y) values.
    """
    if nmax is None:
        nmax = nd.cfg.maxpolygonpoints
    assert 3 < nmax < 8191
    # Most polygons don't need anything. Make that fast.
    if len(p) < nmax:
        return [p]
    # First cleanup the polygon by removing consecutive identical points
    if isinstance(p, ndarray):
        p = [v for i, v in enumerate(p) if i == 0 or (v != p[i - 1]).any()]
    else:
        p = [uniq[0] for uniq in groupby(p)]
    todo = [p]
    result = []
    while todo:
        pol = todo.pop()
        npts = len(pol)
        if npts < nmax:
            result.append(pol)
            continue
        # Find a good intersection line:
        # Split horizontal or vertical in half.
        xmin, ymin, xmax, ymax = boundingbox(pol)
        if xmax - xmin > ymax - ymin:  # Split vertical
            s0 = [xmin + (xmax - xmin) / 2, 0]
            s1 = [xmin + (xmax - xmin) / 2, 100]
        else:  # Split horizontal
            s0 = [0, ymin + (ymax - ymin) / 2]
            s1 = [100, ymin + (ymax - ymin) / 2]
        pols = polygon_split(s0, s1, pol)
        for p in pols:
            l = len(p)
            if l < nmax:
                if l > 2:
                    result.append(p)
            else:
                todo.append(p)
    return result


def limit_polyline_points(p, nmax=None):
    """Limit the maximum number in a polyline by splitting it into multiple
    polylines if needed. The returned polyline(s) will have at most nmax
    points. If nmax is not specified, then he maximum number comes from the cfg
    module: maxpolylinepoints.

    Args:
        p (list): polyline with list of (x, y) values.
        nmax (int): maximum number of points in resulting polyline(s).

    Returns:
        list of polylines, that each are a list of (x, y) values.
    """
    result = []
    if nmax is None:
        nmax = nd.cfg.maxpolylinepoints
    assert 2 < nmax < 8191
    npts = len(p)
    # Segments should overlap (have two points in common) to avoid gaps.
    Nsegment = int(npts / (nmax - 2)) + 1
    # Points per segment
    N = int(npts / Nsegment) + 1
    for i in range(Nsegment):
        # Start and end segment should overlap to avoid gaps.
        start = max(i * N - 2, 0)
        stop = min((i + 1) * N, npts)
        result.append(p[start:stop])
        if stop >= npts:
            return result
