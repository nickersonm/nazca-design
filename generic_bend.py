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
#
# 2018 (c) X. Leijtens, Ronald Broeke
# -----------------------------------------------------------------------
# -*- coding: utf-8 -*-

"""
Module to generate parametric waveguide curves or 'pcurves' and sine bends.

Polyline approximations of the curve are generated.

The pcurves are based on Francois Ladouceur and Pierre Labeye,
"A New General Approach to Optical Waveguide Path Design",
Journal of Lightwave Technology, vol. 13, no. 3, march 1995.

The algorithm tries to construct a curve that matches the input and output
radius, and position, while maximizing the radius along the curve. Patent
expired.

The radius arguments in this module are sensitive to the sign. A positive
radius indicates a counter-clockwise direction.
"""

from scipy.optimize import fminbound, minimize_scalar
from math import sin, cos, radians, sqrt, log, pi
import nazca as nd
from nazca.util import polyline_length


def distance_point2line(P, A, B):
    """Calculate the shortest distance of point P to a line through A and B.

    Note that the shortest distance is the perpendicular distance
    from P to line AB.

    Args:
        P (float, float): (x, y) target point
        A (float, float): (x, y) first point on line
        B (float, float): (x, y) second point on line

    Returns:
        float: distance
    """
    a = [P[0] - A[0], P[1] - A[1]]
    b = [B[0] - A[0], B[1] - A[1]]
    ab = a[0] * b[0] + a[1] * b[1]
    bb = b[0] * b[0] + b[1] * b[1]
    D = [a[0] - ab / bb * b[0], a[1] - ab / bb * b[1]]
    return sqrt(D[0] * D[0] + D[1] * D[1])


X = list()
Z = list()


def getCurvature(A, B, Lp):
    """Get the local curvature at a given point Lp of the parametrized curve.

    P(Lp) = (A(Lp), B(Lp)) along the P-Curve.

    For local curvature see:
    http://en.wikipedia.org/wiki/Curvature#Local_expressions

    Args:
        A (list of float): six coefficients
        B (list of float): six coefficients
        Lp (float): scaled curve-parameter point

    Returns:
        float: local curvature
    """
    dx = ddx = dy = ddy = 0
    for i in range(1, 6):
        # Because x=A[i]*t^i -> dx/dL = ... dy/dL works similarly
        dx += i * A[i] * Lp ** (i - 1)
        dy += i * B[i] * Lp ** (i - 1)
        # and (d^2 x)/(d^2 L) = ... dy/dt works similarly
        if i > 1:
            ddx += i * (i - 1) * A[i] * Lp ** (i - 2)
            ddy += i * (i - 1) * B[i] * Lp ** (i - 2)

    return (ddx * dy - ddy * dx) * ((dx ** 2 + dy ** 2) ** (-1.5))


def invert_matrix(Lp):
    """Helper function returning the solution matrix as a function of Lp.

    Args:
        Lp(float): scaled curve-parameter point

    Returns:
        tuple of float: 6x6 matrix
    """
    # This is the correct matrix. Please note that there is an error in the
    # original paper (entry [3, 6])
    matrix = \
       ((        1,        0,        0  ,         0,        0,          0),
        (        0,        1,        0  ,         0,        0,          0),
        (        0,        0,        0.5,         0,        0,          0),
        (-10/Lp**3, -6/Lp**2,    -3/2/Lp,  10/Lp**3, -4/Lp**2,     1/2/Lp),
        ( 15/Lp**4,  8/Lp**3,  3/2/Lp**2, -15/Lp**4,  7/Lp**3,   -1/Lp**2),
        ( -6/Lp**5, -3/Lp**4, -1/2/Lp**3,   6/Lp**5, -3/Lp**4,  1/2/Lp**3))
    return matrix


def curve_AB(Lp):
    """Calculate array A and array B for given L.

    Args:
        Lp (float): scaled curve-parameter point

    Returns:
        (list of float, list of float): six elements A, six elements B
    """
    global X, Z
    matrix = invert_matrix(Lp)

    # Compute the coefficient
    A = [0, 0, 0, 0, 0, 0]
    B = [0, 0, 0, 0, 0, 0]
    for i in range(6):
        for t in range(6):
            A[i] = A[i] + matrix[i][t] * X[t]
            B[i] = B[i] + matrix[i][t] * Z[t]
    return (A, B)


def max_curvature(L):
    """Calculate the maximum curvature along L.

    There is an ambiguity when the maximum curvature occurs at the input
    or output position. Many paths in between then satisfy the curvature
    requirement. We therefore multiply by log(length) to favor short
    connections without affecting the result much.

    Args:
        L (float): scaled curve-parameter size

    Returns:
        float: maximum curvature
    """
    # Only this many points for the maximum curvature finding.
    Npoints = 100

    global X, Z
    A, B = curve_AB(L)
    curvature = [getCurvature(A, B, L / (Npoints - 1) * i) for i in range(Npoints)]
    maxcur = abs(max(curvature, key=abs))

    return maxcur * log(L)  # max curvature


def gb_coefficients(xya, radius1=0, radius2=0):
    """Calculate coefficients A, B and L.

    Determine the optimum curve length so that the minimum Radius along
    the curve is maximized. The center of the curve is defined by a polynomial
    with coefficients A[] for x, and B[] for y.
    The independent variable runs from 0 to L.

    Args:
        xya (tuple): (x, y, a) coordinate of final position
        radius1 (float): radius of curvature at start
        radius2 (float): radius of curvature at end

    Returns:
        list of float, list of float, float: array of six A, array of six B, L
    """

    global X, Z
    # In the algorithm, a negative curvature indicates a counterclockwise
    # direction, while Nazca uses the sign of the angle instead (negative
    # angle in Nazca = clockwise direction). We therefore invert the
    # arguments here.
    _radius1 = -radius1
    _radius2 = -radius2
    # Input curvature (when radius larger than a nanometer)
    if abs(_radius1) > 1e-3:
        cin = 1 / _radius1
    else:
        cin = 0
    # Output curvature
    if abs(_radius2) > 1e-3:
        cout = 1 / _radius2
    else:
        cout = 0

    xin = 0
    zin = 0
    thetain = 90

    xout = xya[0]  # final position
    zout = xya[1]
    thetaout = -(xya[2] - 90.0)  # final slope [deg]

    # Search between cartesian distance (straight line between input and
    # output) Note that this length is related to the physical length, but
    # not equal.
    minFindL = sqrt(zout ** 2 + xout ** 2) / 4
    # and 3* the length of a straight line
    maxFindL = sqrt(zout ** 2 + xout ** 2) * 3

    # Calculate initial known terms
    dxin = sin(radians(thetain))
    dzin = cos(radians(thetain))
    ddxin = cin * cos(radians(thetain))
    ddzin = -cin * sin(radians(thetain))

    dxout = sin(radians(thetaout))
    dzout = cos(radians(thetaout))
    ddxout = cout * cos(radians(thetaout))
    ddzout = -cout * sin(radians(thetaout))

    # Vectors of known terms
    X = [xin, dxin, ddxin, xout, dxout, ddxout]
    Z = [zin, dzin, ddzin, zout, dzout, ddzout]

    # Determine the optimum curve length so that the minimum radius along
    # the curve is maximized
    L = fminbound(max_curvature, minFindL, maxFindL, disp=1)
    Rmin = 1 / (max_curvature(L) / log(L))  # minimum radius

    # Get the final matrix and extract the polynomial coefficients from it
    A, B = curve_AB(L)
    # The center of the curve is defined by a polynomial with coefficients
    # A[] for x, and B[] for y and the independent variable runs from 0 to L.
    return A, B, L, Rmin


def gb_point(t, A, B, Lp):
    """Calculate generic-bend point for parameter t.

    Args:
        t (): normalized parameter in range [0, 1].
        A (list of float): array of six.
        B (list of float): array of six.
        Lp (float): scaled curve-parameter point, scales t.

    Returns:
        (float, float): generic bend point
    """
    x = y = 0
    for i in range(6):
        x += A[i] * (t * Lp) ** i
        y += B[i] * (t * Lp) ** i
    return (x, y)


def sinebend_point(t, distance, offset):
    """Calculate sine bend point for parameter t. This bend has zero curvature
    at both ends. It provides an S-shape connection.

    Args:
        t (float): normalized parameter in range [0, 1].
        distance (float): scaling factor of x with t.
        offset (float): scaling factor of y with t.

    Returns:
       (float, float): point (x, y) for given t
    """
    x = t * distance
    if t > 0:
        y = offset * t - offset / (2 * pi) * sin(2 * pi * t)
    else:
        y = 0
    return x, y


def curve2polyline(fie, xya, acc, args=()):
    """Generate a polyline from a parameterized curve.

    Use a sufficient number of points to ensure that the deviation of the
    sampled curve from the real curve is less than the specified accuracy.

    Args:
        fie (function): the curve function that takes one parameter that
            runs from 0 to 1 as first argument. The curve starts at the
            origin at 0 angle.
        xya (tuple): the position (x, y, a) of the end point of the curve.
            Angle a in degrees.
        acc (float): desired accuracy in micrometer.
        args (tuple): additional arguments to be passed to the curve function.

    Returns:
       list of (float, float): list of points (x, y);
           a polyline approximation of the curve
    """
    # As a starting point find the value for t where the y-coordinate is
    # equal to acc. Since the curve always starts horizontal, this gives a
    # good and easy value.
    def fun(x):
        return abs(abs(fie(x, *args)[1]) - acc)

    res = minimize_scalar(fun, bounds=(0, 0.04), method="bounded")
    dt0 = res.x
    t0 = dt0 / 2

    P1 = fie(0, *args)
    p = [P1]
    P1 = fie(t0, *args)
    p.append((P1[0], 0))  # First section has to be horizontal.
    while t0 + dt0 < 1:
        P0 = P1
        P1 = fie(t0 + dt0, *args)
        P2 = fie(t0 + 2 * dt0, *args)
        d = distance_point2line(P1, P0, P2)
        f = sqrt(4 * acc / d)
        dt1 = max(0.8, min(f, 1.2)) * dt0
        if t0 + dt1 < 1:
            P1 = fie(t0 + dt1, *args)
            p.append(P1)
        t0 += dt1
        dt0 = dt1
    if t0 + dt0 - 1 > 0.4 * dt0:
        # Make room for the point with the right direction.
        p.pop()
    t0 = 1 - dt0 / 3
    P1 = fie(t0, *args)
    P2 = fie(1, *args)
    lv = sqrt((P2[0] - P1[0]) ** 2 + (P2[1] - P1[1]) ** 2)
    a = radians(xya[2])
    P1 = (P2[0] - lv * cos(a), P2[1] - lv * sin(a))
    p.append(P1)
    p.append(P2)
    return p


if __name__ == "__main__":

    # From (0,0,0) to (x,y,a)
    xya = (1000, 1000, 20)

    # Curve with different accuracy. Placed at an angle to prevent klayout
    # from removing points which are needed to check the number of points.
    A, B, L, Rmin = gb_coefficients(xya, radius1=0, radius2=0)
    xy = curve2polyline(gb_point, xya, 0.1, (A, B, L))
    polygons = nd.util.polyline2polygons(xy, width=2)
    for polygon in polygons:
        nd.Polygon(layer=202, points=polygon).put(0, 0, 10)
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    polygons = nd.util.polyline2polygons(xy, width=2)
    for polygon in polygons:
        nd.Polygon(layer=203, points=polygon).put(0, 0, 10)

    # Curve with different accuracy and curvature at start and finish.
    A, B, L, Rmin = gb_coefficients(xya, radius1=300, radius2=-100)
    xy = curve2polyline(gb_point, xya, 0.1, (A, B, L))
    nd.Polyline(layer=1, width=2, points=xy).put(0, 0, 10)
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    nd.Polyline(layer=2, width=2, points=xy).put(0, 0, 10)

    # Curves with different "length".
    A, B, L, Rmin = gb_coefficients(xya, radius1=0, radius2=0)  # sets X and Z
    for L in range(300, 2000, 100):
        A, B = curve_AB(L)
        xy = [gb_point(t / 1000, A, B, L) for t in range(1001)]
        nd.Polyline(layer=2, width=2, points=xy).put(500, 0, 10)

    A, B, L, Rmin = gb_coefficients(xya, radius1=300, radius2=-100)
    for L in range(300, 2000, 100):
        A, B = curve_AB(L)
        xy = [gb_point(t / 1000, A, B, L) for t in range(1001)]
        nd.Polyline(layer=3, width=2, points=xy).put(500, 0, 10)

    # Length calculation
    A, B, L, Rmin = gb_coefficients(xya, radius1=300, radius2=-100)
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    nd.Polyline(points=xy, width=2).put(-300, 100, 0)
    length = polyline_length(xy)
    print("Length of curve: {:.3f} µm".format(length))

    # Varying width
    def varwidth(t):
        """Function to define the width of the waveguide along the length, when
        the independent variable runs from 0 to 1.
        """
        return 2 * (1.1 + sin(10 * pi * t))

    A, B, L, Rmin = gb_coefficients(xya, radius1=300, radius2=-100)
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    polygons = nd.util.polyline2polygons(xy, width=varwidth)
    for polygon in polygons:
        nd.Polygon(layer=68, points=polygon).put(-300, 150, 0)

    # Sine bend
    xy = curve2polyline(sinebend_point, (100, 20, 0), 0.001, (100, 20))
    nd.Polyline(layer=2, width=2, points=xy).put(-300, 0, 0)
    polygons = nd.util.polyline2polygons(xy, width=2)
    for polygon in polygons:
        nd.Polygon(layer=203, points=polygon).put(-300, 0, 0)

    nd.export_gds()
