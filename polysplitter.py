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
from operator import itemgetter

__all__ = ["limit_polygon_points"]


# https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
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


def intersect(a, b, c, d, firstisline=False):
    """Intersect AB with CD

    Args:
        a (list): (x, y) of start of first line segment
        b (list): (x, y) of end of first line segment
        c (list): (x, y) of start of second line segment
        c (list): (x, y) of end of second line segment
        firstisline (boolean): treat first segment as full line (default False)

    Return:
        list of length 0, 1 or 2 that contains the intersection point(s)
    """
    eps = 1e-5  # Unit would be um or nm, so this is safe?
    # Define vectors p, p+r and q, q+r
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
            type = 0
            if firstisline or -eps < t0 < 1 + eps:
                pi = p + r * t0
                result.append(
                    (pi.x, pi.y),
                )
                type = 1
            if firstisline or -eps < t1 < 1 + eps:
                pi = p + r * t1
                result.append((pi.x, pi.y))
                type = 2 + type
            # type = 0 (no), 1 (1st point), 2 (2nd point), 3 (1st and 2nd point)
            return type, result
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
                return 4, [(pi.x, pi.y)]  # Starts at line (segment)
            elif abs(u - 1) < eps:
                return 5, [(pi.x, pi.y)]  # Ends at line (segment)
            else:
                return 6, [(pi.x, pi.y)]  # Intersects line (segment)
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
    This is an implementation of the algorithm outlined by David Geier:
    https://geidav.wordpress.com/2015/03/21/splitting-an-arbitrary-polygon-by-a-line/
    But his algorithm doesn't properly take into account coinciding points on
    the split line. This one does.

    Args:
        s0 (float, float): (x, y) coordinate of segment start.
        s1 (float, float): (x, y) coordinate of segment end.
        p (list): polygon, a list of (x, y) coordinates.

    Returns:
        list of polygons that are each a list of (x, y) coordinates.
    """
    # To prevent modifying the original polygon, we first make a copy.
    p = p[:]
    # Ensure it is closed.
    if p[0] != p[-1]:
        p.append(p[0])
    if nd.clipper.signed_area(p) > 0:  # Clockwise
        p.reverse()

    s0x, s0y = s0
    s1x, s1y = s1

    index_on_line = []  # index of edges on intersection line.
    # Since we will modify the polygon as we go, we cannot just loop over the
    # elements, but use indices in the while loop:
    l = len(p) - 1
    i = 0
    while i < l:
        # Candidate segments
        type, points = intersect(s0, s1, p[i], p[i + 1], True)
        # print(f"segment p{i}-p{i+1}: type {type}")
        if type in {1, 3, 4}:
            index_on_line.append(i)
        elif type == 6:
            p.insert(i + 1, points[0])
            l += 1
        i += 1

    # Now we'll use a linked list to represent the polygon. This is needed to
    # later be able to split off separate parts.
    l = len(p)  # length of the updated polygon
    # Make list of nodes (should not have the closing point)
    node = [Node(i) for i in range(l - 1)]  # Maps the updated polygon
    l -= 1
    # Link them
    for i in range(l):
        node[i - 1].next = node[i]
        node[(i + 1) % l].prev = node[i]
        node[i].prev = node[i - 1]
        node[i].next = node[(i + 1) % l]

    # Sort on distance from start of the dividing line segment.
    def signed_distance(i):
        return (p[i][0] - s0[0]) * (s1[0] - s0[0]) + (p[i][1] - s0[1]) * (s1[1] - s0[1])

    def lor(node):
        """Return L, O, R"""
        i = node.data
        d = (p[i][0] - s0[0]) * (s1[1] - s0[1]) - (p[i][1] - s0[1]) * (s1[0] - s0[0])
        if d < 0:
            return "L"
        if d > 0:
            return "R"
        return "O"

    # Also replace the index stored in index_on_line with a tuple:
    # - node
    # - signed distance
    # - LOR
    on_line_1 = []  # All points on intersection line
    for ndx in index_on_line:
        nod = node[ndx]
        LOR = lor(nod.prev) + "O" + lor(nod.next)
        on_line_1.append((nod, signed_distance(ndx), LOR))
    on_line_2 = []  # after first filter
    for nod, d, LOR in on_line_1:
        if LOR in {"LOR", "ROL", "OOL", "OOR", "LOO", "ROO"}:
            on_line_2.append((nod, d, LOR))

    # Sort on distance.
    on_line_2.sort(key=itemgetter(1))
    # Special treatment for xOO-OOy and OOx-yOO: replace by LOR or ROL.
    # Take care that ponts may be in the wrong order, when they coincide.
    # To account for strange OOO cases, use != "L" in stead of == "R".
    on_line = []  # Store points on interesection line
    done = set()  # Don't store the same point twice
    i = 0
    while i < len(on_line_2):
        nod, d, LOR = on_line_2[i]
        if nod in done:
            i += 1
            continue
        if "OO" not in LOR:
            on_line.append((nod, d, LOR))
        elif LOR == "OOL":
            if lor(nod.prev.prev) != "L":
                on_line.append((nod, d, "ROL"))
                done.update({nod, nod.prev})
        elif LOR == "ROO":
            if lor(nod.next.next) != "R":
                on_line.append((nod, d, "ROL"))
                done.update({nod, nod.next})
        elif LOR == "LOO":
            if lor(nod.next.next) != "L":
                done.update({nod, nod.next})
                # Can't take from list because order might be wrong
                nod = nod.next
                d = signed_distance(nod.data)
                on_line.append((nod, d, "LOR"))
        elif LOR == "OOR":
            if lor(nod.prev.prev) != "R":
                done.update({nod, nod.prev})
                # Can't take from list because order might be wrong
                nod = nod.prev
                d = signed_distance(nod.data)
                on_line.append((nod, d, "LOR"))
        i += 1

    # Only LOR and ROL in on_line and they should be in pairs, LOR-ROL
    lorrol = ["LOR", "ROL"]
    for i in range(len(on_line) - 1):
        node, d, LOR = on_line[i]
        if LOR != lorrol[i % 2]:
            # Not in the right order, d must be the same and swap them
            if d == on_line[i + 1][1]:
                on_line[i], on_line[i + 1] = on_line[i + 1], on_line[i]
            else:
                raise ValueError("Self-intersecting polygon")

    # Now guaranteed to be correct (LOR-ROL) * N.
    for a, b in zip(on_line[::2], on_line[1::2]):
        src_node, dst_node = a[0], b[0]
        # Bridge source and destination
        # print(f"Source {src_node.data}, destination {dst_node.data}")
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
    todo = [p]
    result = []
    if nmax is None:
        nmax = nd.cfg.maxpolygonpoints
    assert 3 < nmax < 8191
    while todo:
        pol = todo.pop()
        npts = len(pol)
        if npts < nmax:
            result.append(pol)
            continue
        # Find a good intersection line:
        # Split horizontal or vertical in half.
        xmax = xmin = pol[0][0]
        ymax = ymin = pol[0][1]
        for x, y in pol:
            if x > xmax:
                xmax = x
            elif x < xmin:
                xmin = x
            if y > ymax:
                ymax = y
            elif y < ymin:
                ymin = y
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
