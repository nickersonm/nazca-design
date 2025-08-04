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
# Geometry that defines a cross-section
#
# (c) 2022 Lodovico Rossi, Marco Passoni
#
import copy
from collections import namedtuple
import sys
from functools import lru_cache
from typing import Any
from functools import partial
import math as m
import numpy as np
import yaml
import matplotlib
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
from matplotlib.path import Path
import matplotlib.pyplot as plt
import subprocess
import os

import nazca
from nazca.util import (
    transform_polygon,
    _trapezoid,
    _intersect,
    linvarwidth,
    parabolicvarwidth,
    filter_eval,
)
from nazca.clipper import grow_polygons
from nazca import main_logger


PolygonData = namedtuple("PolygonData", ["vertices", "level", "area"])


def get_polygon_area(x: list, y: list):
    """Returs the area of a polygon

    Args:
        x (list): List of x coordinates of the polygon points
        y (list): List of y coordinates of the polygon points

    Returs:
        float: area of the polygon
    """
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def __run_cmd(command, folder=None):
    """Run command in the terminal and returs the output

    Args:
        command (srt): Command to be run.

    Returns:
        str: output of the command
    """
    command = command.split()
    out = (
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            cwd=folder,
        )
        .stdout.strip()
        .decode("utf-8")
    )
    return out


def _find_args_dict(func: callable) -> dict:
    """Returns dictionary of arguments of a function."""
    try:
        defaults = {i: value for i, value in enumerate(func.__defaults__[::-1])}
    except TypeError:
        defaults = {}
    return {
        arg: defaults.get(i, None)
        for i, arg in enumerate(
            func.__code__.co_varnames[func.__code__.co_argcount - 1 :: -1]
        )
    }


def _get_function_git(func: Any) -> dict:
    """Given a funciton, returs the data of the git repository where it is defined

    Args:
        func (any): python object

    Returns:
        dict: dictionary of git data
    """
    try:
        path, file = os.path.split(sys.modules[func.__module__].__file__)
    except AttributeError:
        return "Cannot access file where function is defined (possible cause: funciton is defined in jupyter notebook)"
    version = __run_cmd("git --version")
    if "version" not in version:
        return "Info not available"
    try:
        git_root = __run_cmd("git rev-parse --show-toplevel", folder=path)
    except subprocess.CalledProcessError:
        return "Function not in a git repository"
    except FileNotFoundError:
        return "Cannot access folder where function is defined (possible cause: funciton is defined in jupyter notebook)"
    git_data = {
        "Root folder": git_root.replace("\n", ""),
        "Root commit sha": __run_cmd(
            "git rev-list --max-parents=0 HEAD", folder=path
        ).replace("\n", ""),
    }
    modifications = __run_cmd("git status --porcelain", folder=path).split("\n")
    modifications = [x.split() for x in modifications]
    git_data["Git Clean"] = not bool([_ for _ in modifications if _ and _[0] == "M"])
    git_data["Untracked Files"] = bool([_ for _ in modifications if _ and _[0] == "??"])
    git_data["Modified function file"] = bool(
        __run_cmd(f"git diff {file}", folder=path)
    )
    try:
        des = __run_cmd("git describe --long", folder=path)
    except subprocess.CalledProcessError:
        des = __run_cmd("git describe --all --long", folder=path)
    git_data["git describe"] = des.replace("\n", "")
    return git_data


def _get_function_data(func: callable) -> dict:
    """Returns summary of function data."""
    data = {
        "Function name": func.__name__,
        "Function file": os.path.split(func.__code__.co_filename)[1],
        "Function path": os.path.split(func.__code__.co_filename)[0],
        "Function line number": func.__code__.co_firstlineno,
        "Arguments": _find_args_dict(func),
        "Git Repository": _get_function_git(func),
    }
    return data


def _check_call(obj: Any, *args, **kwargs) -> Any:
    """Returns call to obj if obj is callable, else object"""
    try:
        return obj(*args, **kwargs)
    except TypeError:
        return obj


def _process_iterable(obj: Any, *args, **kwargs) -> Any:
    """Processes an iterable and call any object inside that is callable"""
    try:
        return type(obj)([_process_iterable(x, *args, **kwargs) for x in obj])
    except TypeError:
        return _check_call(obj, *args, **kwargs)


def _f_x(
    y: float, x1: float, x2: float, local_width: float, xs_width: float, angle: float
) -> float:
    """Calculates the 3D x coordinate given the y coordinate of the cross-section and the plane x points.

    Args:
        y (float): y coordinate of the cross-section
        x1 (float): x coordinate of the first point of the segment
        x2 (float): x coordinate of the second point of the segment
        local_width (float): Local width of the interconnect
        xs_width (float): Width of the cross-section
        angle (float): local angle in radians.

    Returns:
        float: the 3D x coordinate of the cross-section
    """
    _local_width = (local_width - 1.0) * np.sin(angle) + 1.0

    return ((x1 - x2) * y + _local_width * (x1 + x2)) / xs_width


def _f_y(
    y: float, y1: float, y2: float, local_width: float, xs_width: float, angle: float
) -> float:
    """Calculates the 3D y coordinate given the y coordinate of the cross-section.

    Args:
        y (float): y coordinate of the cross-section
        y1 (float): y coordinate of the first point of the segment
        y2 (float): y coordinate of the second point of the segment
        local_width (float): Local width of the interconnect
        xs_width (float): Width of the cross-section
        angle (float): local angle in radians.

    Returns:
        float: the 3D x coordinate of the cross-section
    """
    _local_width = (local_width - 1.0) * np.cos(angle) + 1.0

    return ((y1 - y2) * y + _local_width * (y1 + y2)) / xs_width


def _det(a: tuple, b: tuple) -> float:
    """Calculates the determinant of the 2x2 matrix with row a and b

    Args:
        a (tuple): (float, float) first row of the matrix
        b (tuple): (float, float) second row of the matrix

    Returns:
        float: Determinant of the matrix
    """
    return a[0] * b[1] - a[1] * b[0]


def _check_inside(point: tuple, A: tuple, B: tuple) -> bool:
    """Check if a point lies inside the edge bbox

    Args:
        point (tuple): (float, float), point to be checked
        A (tuple): (float, float), starting point of the edge
        B (tuple): (float, float), ending point if the edge

    Returns:
        bool: True if point is inside the bbox, false otherwise
    """
    tol = 1e-9
    check = (min(A[0], B[0]) - tol) <= point[0] <= (max(A[0], B[0]) + tol) and (
        min(A[1], B[1]) - tol
    ) <= point[1] <= (max(A[1], B[1]) + tol)
    return check


def _get_distance(a: tuple, b: tuple) -> float:
    """Returns the distance between two points

    Args:
        a (tuple): (float, float). First point.
        b (tuple): (float, float). Second point.

    Returns:
        float: distance
    """
    return np.sqrt((a[0] - b[0]) ** 2.0 + (a[1] - b[1]) ** 2.0)


def _order(points: list, a: tuple) -> list:
    """Order a list of point in increasing distance from a reference

    Args:
        points (list): list of (float, float). List of the points to be ordered
        a (tuple): (float, float). Point to be used as a reference

    Returns:
        list: reordered point list
    """

    distance = [_get_distance(point, a) for point in points]
    ind = np.argsort(distance)
    return [points[i] for i in ind]


def _intersect_edges(a: tuple, b: tuple, c: tuple, d: tuple) -> list:
    """Calculates intersection points between two edges

    This function is asymmetric when edges ab and cd partially coincide.
    In that case, the extremities of cd that also lie inside ab are returned.

    Args:
        a (tuple): (float, float). Starting point of the first edge.
        b (tuple): (float, float). Ending point of the first edge.
        c (tuple): (float, float). Starting point of the second edge.
        d (tuple): (float, float). Ending point of the second edge.

    Returns:
        list: list of (float, float). Intersection points between ab and cd edges
    """
    xdiff = (a[0] - b[0], c[0] - d[0])
    ydiff = (a[1] - b[1], c[1] - d[1])
    div = _det(xdiff, ydiff)
    dd = (_det(a, b), _det(c, d))
    d1 = _det(dd, xdiff)
    d2 = _det(dd, ydiff)
    if div == 0:
        return [point for point in [c, d] if _check_inside(point, a, b)]
    else:
        point = (d1 / div, d2 / div)
        if _check_inside(point, a, b) and _check_inside(point, c, d):
            for ref in [a, b, c, d]:
                if _get_distance(point, ref) < 1e-6:
                    return [ref]
            return [point]
        else:
            return []


def _remove_zero_edges(edges: list) -> list:
    """Removes edges of length 0

    Args:
        edges (list): list of ((float, float), (float, float)). List of edges to be filtered.

    Returns:
        list: same as input list, without the length 0 edges.
    """
    return [edge for edge in edges if not np.allclose(edge[0], edge[1])]


def _separate_single_edge(edge: tuple, edge_list: list) -> list:
    """Split and edge into multiple ones.

    The edge is split in the crossing point with a list of other edges.

    Args:
        edge (tuple): ((float, float), (float, float)). Edge to be split.
        edge_list (list): list of ((float, float), (float, float)). List of edges according to which to split.

    Returns:
        list: list of ((float, float), (float, float)). New edges dividing the input one.
    """

    intersection_points = []
    for edge_2 in edge_list:
        intersection_points += _intersect_edges(*edge, *edge_2)
    intersection_points = _order(intersection_points, edge[0])
    intersection_points = [edge[0]] + intersection_points + [edge[1]]
    new_edges = _remove_zero_edges(zip(intersection_points, intersection_points[1:]))
    return new_edges


def _separate_edges(edges: list) -> list:
    """Separate a collection of edges into a collection of non-intersecting ones.

    Args:
        edges (list): list of ((float, float), (float, float)). List of edges to be separated.

    Returns:
        list: Same format as input. List of edges in the new representation.
    """
    new_edges = []
    for i, edge in enumerate(edges):
        _edges = edges.copy()
        _edges.pop(i)
        new_edges += _separate_single_edge(edge, _edges)
    new_edges = _remove_zero_edges(new_edges)
    new_edges = list(set(new_edges))
    return new_edges


def _get_perpendicular_edge(a: tuple, b: tuple, length: float = 0.001) -> tuple:
    """Get and edge cutting the first one in half

    Args:
        a (tuple): (float, float). Starting point of the first edge.
        b (tuple): (float, float). Ending point of the first edge.
        length (float): Length of the cutting edge

    Returns:
        tuple: ((float, float), (float, float)). Perpendicular edges.
    """
    dist = np.sqrt((a[0] - b[0]) ** 2.0 + (a[1] - b[1]) ** 2.0)
    direction = (
        (b[0] - a[0]) / dist,
        (b[1] - a[1]) / dist,
    )
    direction = (direction[1], -direction[0])
    middle = (
        0.5 * (b[0] + a[0]),
        0.5 * (b[1] + a[1]),
    )
    aprime = (
        middle[0] + direction[0] * 0.5 * length,
        middle[1] + direction[1] * 0.5 * length,
    )
    bprime = (
        middle[0] - direction[0] * 0.5 * length,
        middle[1] - direction[1] * 0.5 * length,
    )
    return aprime, bprime


def _clockwiseangle_and_distance(point: tuple):
    """Calculates the angle and the length of a point in the 2D plane

    Args:
        point (tuple): The point (x,y)

    Returns:
        float: angle with respect to vector (1.0, 0.0)
        float: distance  with respect to origin

    """
    lenvector = np.hypot(point[0], point[1])
    if lenvector == 0:
        return -np.pi, 0

    angle = np.arctan2(point[1], point[0])
    if angle < 0:
        return 2 * np.pi + angle, lenvector
    return angle, lenvector


def _argsort(seq, key=None):
    _key = seq.__get_item__ if key is None else (lambda _: key(seq[_]))
    return sorted(range(len(seq)), key=_key)


class Geometry:
    """Class that contains all essential information about the geometry of a cross-section.

    The geometry is a collection of shape classes placed in the exact locations of the yz plane.
    """

    active_geometries = []

    def __init__(
        self,
        name="",
        background_index=1.0,
        xs_width=1.0,
        bbox=None,
        padding=2.0,
        ridge_x_coordinate=0.0,
    ):
        """Initializes the class.

        Args:
            name (str): Name of the structure
            background_index (float | complex): Refractive index of the background
            xs_width (float): Width of the cross-section in [um]. Default is 1.
            bbox (tuple): definition of the bounding box in the form:
                    (x_lower_left, y_lower_left, x_upper_right, y_upper_right).
                    Default is None, meaning the bbox is automatically generated.
            padding (float): distance to add around any defined shape when auto generating the bbox.
                It is ignored if bbox is given. Default to 2.0

        Returns:
            None
        """
        self.xs_width = xs_width
        self.name = name
        if callable(background_index):
            self.background_index = filter_eval(background_index)
        else:
            self.background_index = background_index
        self.shapes = []  # List of Instance
        self.bbox = bbox
        self.padding = padding
        self.ridge_x_coordinate = ridge_x_coordinate

    def __enter__(self):
        Geometry.active_geometries.append(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        Geometry.active_geometries.pop()
        if self.bbox is None:
            self._calculate_bbox()

    def get_material_functions(self):
        functions = set()
        try:
            functions.add(self.background_index.__wrapped__)
        except AttributeError:
            functions.add(self.background_index)
        for shape in self.shapes:
            try:
                functions.add(shape.material_index.__wrapped__)
            except AttributeError:
                functions.add(shape.material_index)
        return functions

    def _calculate_bbox(self):
        """Auto generates the bounding box."""
        points = [p for shape in self.shapes for p in shape.get_points()]
        if len(points) == 0:
            x, y = [0.0], [0.0]
        else:
            x, y = list(zip(*points))
        self.bbox = [
            np.min(x) - self.padding,
            np.min(y) - self.padding,
            np.max(x) + self.padding,
            np.max(y) + self.padding,
        ]

    def process_material_properties(self, **kwargs):
        """Returns a copy of self, processing the material functions with kwargs

        Args:
            **kwargs: keyword arguments to be passed to material functions

        Returns:
            Geometry: Same geometry as self, but with the material indices as floats
        """
        xs = copy.deepcopy(self)
        shapes = set([_.shape for _ in xs.shapes])
        for shape in shapes:
            shape.material_index = shape.get_material_index(**kwargs)
        xs.background_index = xs.get_background_index(**kwargs)
        return xs

    def add_shape(self, instance):
        """Adds a shape to the active_geometries pointer.

        Args:
            instance (Instance): Shape to be added to the geometry.

        Returns:
            None
        """
        self.shapes.append(instance)

    def get_edges_representation(self, remove_unnecessary_edges=True, **kwargs):
        """Return the edges representation of the cross section

        Returns a list of edges. Each edge is a tuple of 3 tuples:
        (n_left, n_right): index on the left and on the right of the edge
        (x_start, y_start): coordinates of the starting point og the edge
        (x_end, y_end): coordinates of the ending point og the edge

        Args:
            remove_unnecessary_edges (bool): if True, remove edges that have the same left and right indexes.
            **kwargs: Extra arguments to the passed to the index functions

        Returns:
            list: list of ((float, float), (float, float), (float, float)). List of the edges.
        """

        edges = self.get_edges()
        edges = _separate_edges(edges)
        perpendicular = [_get_perpendicular_edge(*edge, 0.001) for edge in edges]
        right, left = zip(*perpendicular)
        right = [self._get_index_function(*r) for r in right]
        left = [self._get_index_function(*l) for l in left]
        full_edges = []
        for nout, nin, edge in zip(left, right, edges):
            if remove_unnecessary_edges and nout == nin:
                continue
            full_edges.append(((nout, nin), edge[0], edge[1]))
        return full_edges

    def get_edges(self):
        """Return the list of edges defining the geometry.

        This list is just the collection of the edges of the included shapes

        Returns:
            list: list of ((float, float), (float, float)). List of the edges.
        """
        lines = []
        for shape in self.shapes:
            for x1, x2 in shape.get_edges_representation():
                lines.append((x1, x2))
        return lines

    def show(self, function=np.real, ax=None, cmap="viridis", **kwargs):
        """Plots the geometry.

        Args:
            function (callable): Function that filters the indices in the plot. Default is np.real, showing the real
            part of the cross-section indices.
            ax (matplotlib.Axes): Axes onto which the xsection geometry is plotted. Default is None and the plot is
            shown.
            **kwargs: dictionary of arguments to pass to the background and material index functions.

        Returns:
            matplotlib.Figure: The figure containing the axes.
            matplotlib.Axes: The axis onto which the geometry is plotted.
            matplotlib.Colorbar: The colorbar associated to the axis.
        """
        indices = [self.get_background_index(**kwargs)]
        for shape in self.shapes:
            n = shape.shape.get_material_index(**kwargs)
            if n not in indices:
                indices.append(n)

        # Check if the indices are different
        # TODO: if the index is uniform after the function, plot the uniform, not raise ValueError.
        if np.all(function(indices) == 0):
            raise ValueError(
                f"All the indices evaluate to 0 with the function {function}.\n"
            )

        indices_sorted = sorted(function(indices))
        val = int(
            indices_sorted.index(function(self.get_background_index(**kwargs)))
            * 256
            / (len(indices) - 1)
        )

        cmap = matplotlib.colormaps[cmap]

        if ax is None:
            fig, ax = plt.subplots()
            plot_result = True
        else:
            fig = plt.gcf()
            plot_result = False

        ax.set_facecolor(cmap(val))
        ax.fill(
            [self.bbox[0], self.bbox[2], self.bbox[2], self.bbox[0]],
            [self.bbox[1], self.bbox[1], self.bbox[3], self.bbox[3]],
            color=cmap(val),
        )

        for shape in self.shapes:
            val = int(
                indices_sorted.index(function(shape.shape.get_material_index(**kwargs)))
                * 256
                / (len(indices) - 1)
            )

            _x, _y = zip(*shape.get_points())
            ax.fill(_x, _y, color=cmap(val))

        norm = BoundaryNorm(indices_sorted, cmap.N, extend="both")

        ax.set_xlabel("y [$\mu$m]")
        ax.set_ylabel("z [$\mu$m]")
        ax.set_xlim(self.bbox[0], self.bbox[2])
        ax.set_ylim(self.bbox[1], self.bbox[3])

        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)

        if plot_result:
            plt.tight_layout()

        return fig, ax, cbar

    def show_edges(
        self, ax=None, color="black", linewidth=0.2, remove_unnecessary_edges=True
    ):
        """Show edges of the cross section

        Args:
            ax (matplotlib.Axes): Axes to use for the plot. If None, plot is shown
            color (str): Color for the edges. Default is black.
            linewidth (float): Width of the contour line
            remove_unnecessary_edges (Bool): If True (default) edges between same index are removed.

        Returns:
            matplotlib.Figure: The figure containing the axes.
            matplotlib.Axes: The axis onto which the geometry is plotted.
        """
        if ax is None:
            fig, ax = plt.subplots()
            plot_result = True
        else:
            fig = plt.gcf()
            plot_result = False

        for n, (x1, y1), (x2, y2) in self.get_edges_representation(
            remove_unnecessary_edges
        ):
            ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth)

        ax.set_xlim(self.bbox[0], self.bbox[2])
        ax.set_ylim(self.bbox[1], self.bbox[3])

        if plot_result:
            plt.tight_layout()
            plt.show()

        return fig, ax

    def get_background_index(self, **kwargs):
        """Return background material index.

        If the material index is a function, kwargs can be provided.

        Args:
            **kwargs: dictionary of arguments to pass to the index function.

        Returns:
            float | complex: material index inside the shape
        """
        if callable(self.background_index):
            return self.background_index(**kwargs)
        else:
            return self.background_index

    def get_index(self, y, z, **kwargs):
        """Calculates the index at coordinates y and z of the geometry.

        Args:
            y (numpy.Array): Y coordinate
            z (numpy.Array): Z coordinate
            **kwargs: Optional arguments to pass to the function of the material index.

        Returns:
            np.Array: the index at the desired coordinates.
        """
        if np.shape(z) != np.shape(y):
            return ValueError(
                f"y and z have different shapes {np.shape(y)}!={np.shape(z)}"
            )
        ishape = np.shape(z)
        z = np.reshape(z, -1)
        y = np.reshape(y, -1)
        yz = list(zip(y, z))
        index = np.zeros_like(z, dtype=complex) + self.get_background_index(**kwargs)
        for i, shape in enumerate(self.shapes):
            pol = Path(shape.get_points())
            check = pol.contains_points(yz)
            index[check] = shape.get_material_index(**kwargs)
        return index.reshape(ishape)

    def _get_index_function(self, x, y):
        """Return the index of a specific point

        If the index is a function, the function object is returned.

        Args:
            x (float): x coordinate
            y (float): y coordinate

        Returns:
            Any: index at the coordinates.
        """
        index = self.background_index
        for i, shape in enumerate(self.shapes):
            pol = Path(shape.get_points())
            if pol.contains_point((x, y)):
                index = shape.material_index
        return index

    def polyline2polygon3d(
        self,
        xy,
        width1,
        width2=None,
        parabolic=True,
        miter=0.5,
        anglei=None,
        angleo=None,
    ):
        """Calculates the points of the 3d polygon based on the xs geometry and the polyline.

        Args:
            xy (list): list of (x,y) points that hold the polygon
            width1 (float | list | function): width of the polyline (default 2), or a
                parameterized function w(t) returning the width for the independent
                variable which runs from 0 to 1 from start of the polyline to the
                end, proportional to the length of the polyline segments.
            width2 (float): if width is a number, and if width2 is not None, they
                are interpreted as the start width and end width2 with a parabolic
                change of width vs length.
            parabolic (bool): if begin and end widths are specified as numbers, use
                parabolic tapering, or linear tapering if False (default True).
            miter (float): maximum fraction of the width before an extra point
                is added in outside corners (default 0.5).
            anglei (float): Input angle. (default None).
            angleo (float): Output angle (default None).

        Returns:
            list of (float, float): the 3D object sliced into planes
        """
        n = len(xy)
        if n < 2:
            raise ValueError("polyline2polygon3d: need at least 2 points for polyline.")

        # start angle correction
        N = n
        # 1/3 creates effectively 3 near-equal line segments inside the
        # original 2 of connecting elements: ---|--- -> --.-|-.--
        # with the exact angle in the middle segment and | denoting the element transition
        factor = 1.0 / 3.0
        if anglei is not None:
            l0 = m.hypot(xy[1][0] - xy[0][0], xy[1][1] - xy[0][1])
            xi = xy[0][0] + (l0 / 3.0) * m.cos(m.radians(anglei))
            yi = xy[0][1] + (l0 / 3.0) * m.sin(m.radians(anglei))
        if angleo is not None:
            l0 = m.hypot(xy[n - 2][0] - xy[n - 1][0], xy[n - 2][1] - xy[n - 1][1])
            xo = xy[n - 1][0] + factor * l0 * m.cos(m.radians(angleo + 180))
            yo = xy[n - 1][1] + factor * l0 * m.sin(m.radians(angleo + 180))

        xy2 = []
        for i, point in enumerate(xy):
            if i == N - 1 and angleo is not None:
                n += 1
                xy2.append((xo, yo))
            xy2.append(point)
            if i == 0 and anglei is not None:
                n += 1
                xy2.append((xi, yi))
        xy = xy2

        if isinstance(width1, list):
            width_new = []
            for i, W in enumerate(width1):
                if i == N - 1 and angleo is not None:
                    width_new.append(
                        (width1[N - 1] + factor * (width1[N - 2] - width1[N - 1]))
                    )
                width_new.append(W)
                if i == 0 and anglei is not None:
                    width_new.append((width1[0] + factor * (width1[1] - width1[0])))
            width1 = width_new
        # end angle correction

        # Determine the width at each point of the polyline from the
        # caller-supplied width function (or fixed width)
        # First derive the independent variable from the length of each section.
        length = [0]
        ltot = 0
        for ndx in range(1, n):
            ltot += m.hypot(xy[ndx][0] - xy[ndx - 1][0], xy[ndx][1] - xy[ndx - 1][1])
            length.append(ltot)
        t = [lrun / ltot for lrun in length]  # Normalize to [0, 1]
        # Note: if there is a curvature in the polyline, "length" is a lower boundary.

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
            w = [width1 for k in range(n)]
            dsqrmax = [(miter * width1) ** 2 for k in range(n)]
        elif isinstance(width1, list):
            # List of width: should have the same length as xy
            if len(width1) != n:
                raise ValueError(
                    "polyline2polyhedron: length of xy needs " "to match length of width1."
                )
            # Equal width for left and right edges
            if width2 is None:
                w1 = [_w / 2 for _w in width1]
                w2 = [_w / 2 for _w in width1]
            # Different widths for left and right edges
            elif isinstance(width2, (list, np.ndarray)):
                w1 = width1
                w2 = width2
            else:
                raise ValueError(
                    f"polyline2polyhedron: invalid type {type(width2)}" + \
                    "width2 must be either None, a list or a numpy array"
                )

            dsqrmax = [(miter * w) ** 2 for w in width1]
        elif callable(width1):
            # Equal width for left and right edges
            if width2 is None:
                w1 = [width1(x) / 2 for x in t]
                w2 = [width1(x) / 2 for x in t]
            # Different widths for left and right edges
            elif callable(width2):
                w1 = [width1(x) for x in t]
                w2 = [width2(x) for x in t]
            else:
                raise ValueError(
                    f"polyline2polyhedron: invalid type {type(width2)}" + \
                    "width2 must be either None or a function"
                )

            # the fraction of the width of the line segments that is used to
            # determine if a single point is sufficient to describe the outline, or
            # that two points are needed (miter limit).
            dsqrmax = [(miter * width1(x)) ** 2 for x in t]
        else:
            raise ValueError(
                "polyline2polyhedron: don't know what to do with "
                "this width1 parameter: {}.".format(width1)
            )

        _points_3d = []
        # Loop over the shapes
        for shape in self.shapes:
            # Start with the first two points:
            tr1 = _trapezoid(xy[0], xy[1], w[0], w[1])

            _points_3d.append([])
            _points_3d[-1].append(
                shape.get_3d_points(tr1[3], tr1[0], w[0], self.xs_width)
            )
            # loop over the points in the polyline:
            for i in range(1, n - 1):
                tr0 = tr1  # Current and next trapezoid
                # Get corner points for next segment.
                tr1 = _trapezoid(xy[i], xy[i + 1], w[i], w[i + 1])
                # left or right turn
                lrt = (xy[i + 1][1] - xy[i - 1][1]) * (xy[i][0] - xy[i - 1][0]) - (
                    xy[i + 1][0] - xy[i - 1][0]
                ) * (xy[i][1] - xy[i - 1][1])
                # Distance (squared) between the two points at the kink
                # (Top and bottom have the same distance)
                dsqr = (tr1[0][0] - tr0[1][0]) ** 2 + (tr1[0][1] - tr0[1][1]) ** 2

                if lrt > 0:  # Left turn
                    # Outside corner: use two points, unless these points are close.
                    if dsqr < dsqrmax[i]:
                        _points_3d[-1].append(
                            shape.get_3d_points(
                                _intersect(
                                    tr0[3], tr0[2], tr1[3], tr1[2]
                                ),  # Outside corner
                                _intersect(
                                    tr0[0], tr0[1], tr1[0], tr1[1]
                                ),  # Inside corner: always intersect.
                                w[i],
                                self.xs_width,
                            )
                        )
                    else:
                        _points_3d[-1].append(
                            shape.get_3d_points(
                                tr0[2],  # Outside corner
                                _intersect(
                                    tr0[0], tr0[1], tr1[0], tr1[1]
                                ),  # Inside corner: always intersect.
                                w[i],
                                self.xs_width,
                            )
                        )
                        _points_3d[-1].append(
                            shape.get_3d_points(
                                tr1[3],  # Outside corner
                                _intersect(
                                    tr0[0], tr0[1], tr1[0], tr1[1]
                                ),  # Inside corner: always intersect.
                                w[i],
                                self.xs_width,
                            )
                        )
                elif lrt < 0:  # Right turn
                    # Outside corner: use two points, unless these points are close.
                    if dsqr < dsqrmax[i]:
                        _points_3d[-1].append(
                            shape.get_3d_points(
                                _intersect(
                                    tr0[3], tr0[2], tr1[3], tr1[2]
                                ),  # Inside corner: always intersect.
                                _intersect(
                                    tr0[0], tr0[1], tr1[0], tr1[1]
                                ),  # Outside corner
                                w[i],
                                self.xs_width,
                            )
                        )
                    else:
                        _points_3d[-1].append(
                            shape.get_3d_points(
                                _intersect(
                                    tr0[3], tr0[2], tr1[3], tr1[2]
                                ),  # Inside corner: always intersect.
                                tr0[1],  # Outside corner
                                w[i],
                                self.xs_width,
                            )
                        )
                        _points_3d[-1].append(
                            shape.get_3d_points(
                                _intersect(
                                    tr0[3], tr0[2], tr1[3], tr1[2]
                                ),  # Inside corner: always intersect.
                                tr1[0],  # Outside corner
                                w[i],
                                self.xs_width,
                            )
                        )
                else:
                    _points_3d[-1].append(
                        shape.get_3d_points(tr0[2], tr0[1], w[i], self.xs_width)
                    )

            # Last two points.
            _points_3d[-1].append(
                shape.get_3d_points(tr1[2], tr1[1], w[-1], self.xs_width)
            )

        return _points_3d

    def _get_ridge_layers(self, x):
        """Returns the layer representation along one line at specific x coordinate

        Args:
            x (float): x coordinate where to calculate the layer list
            **kwargs: extra arguments to be passed to the index functions

        Returns:
            list: list of tuples (index, thickness)
        """
        edges = self.get_edges_representation()
        ridge_edge = (
            (x, self.bbox[1]),
            (x, self.bbox[3]),
        )
        intersect_points = []
        for indexes, p1, p2 in edges:
            intersect = _intersect_edges(*ridge_edge, p1, p2)
            if len(intersect) == 2:
                raise ValueError(
                    f"x coordinate {x} is on an edge. Not possible to define a layer stack"
                )
            elif len(intersect) == 1:
                intersect_points.append(intersect[0])
        intersect_points = _order(intersect_points, ridge_edge[0])
        intersect_points = [ridge_edge[0]] + intersect_points + [ridge_edge[1]]
        layers = []
        for p1, p2 in zip(intersect_points, intersect_points[1:]):
            pmiddle = (0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1]))
            thickness = p2[1] - p1[1]
            index = self._get_index_function(*pmiddle)
            if thickness > 1e-9:
                layers.append((index, thickness))
        return layers

    def get_multilayer_representation(self):
        """Return the multiplayer representation of the cross section

        Return a nested list of tuples.
        The order of the slicing is horizontal, the vertical

        Args:
            **kwargs:

        Returns:

        """
        edges = self.get_edges_representation()
        x_coordinates = set()
        for indexes, p1, p2 in edges:
            if np.isclose(p1, p2).any():
                x_coordinates.add(p1[0])
                x_coordinates.add(p2[0])
            else:
                raise ValueError(
                    f"Cross section {self} has slanted edges, no multilayer representation possible"
                )
        x_coordinates = np.sort(list(x_coordinates))
        x_coordinates = [x for x in x_coordinates if self.bbox[0] < x < self.bbox[2]]
        x_coordinates = np.array([self.bbox[0]] + x_coordinates + [self.bbox[2]])
        x_middle = 0.5 * (x_coordinates[:-1] + x_coordinates[1:])
        thickness = x_coordinates[1:] - x_coordinates[:-1]
        layers = [
            (self._get_ridge_layers(x), t)
            for x, t in zip(x_middle, thickness)
            if t > 1e-9
        ]
        return layers

    @classmethod
    def from_function(
        cls,
        function,
        bbox,
        bbox_function=None,
        n_x=51,
        n_y=51,
        levels=11,
        fkwargs=None,
    ):
        """Creates a geometry from a function of 2D variables by creating polygons based on the contour lines.

        Args:
            function (callable): Function of at least x and y
            bbox (array_like): array_like object of 4 elements: the bounding box for the geometry.
                Structure is (x_bottom_left, y_bottom_left, x_top_right, y_top_right)
            bbox_function (array_like): The bounding box used for the evaluaiton of the function. If None, bbox is assumed.
            n_x (int): Number of points along x the function is sampled with
            n_y (int): Number of points along y the function is sampled with
            levels (int): Number of levels the function approximated with.
            fkwargs (dict): Keyword arguments fed to the function

        Returns
            Geometry: The geometry containing the sampled function
        """
        if len(bbox) != 4:
            raise ValueError("bbox needs to be an array like object with lenght of 4")

        bbox_function = bbox if bbox_function is None else bbox_function
        fkwargs = {} if fkwargs is None else fkwargs

        x = np.linspace(bbox_function[0], bbox_function[2], n_x)
        y = np.linspace(bbox_function[1], bbox_function[3], n_y)

        _x, _y = np.meshgrid(x, y, indexing="ij")

        function_map = function(x=_x, y=_y, **fkwargs)
        contours = plt.contour(_x, _y, function_map, levels=levels)
        plt.close()

        levels = contours.levels
        layers = 0.5 * (np.asarray(levels)[:-1] + np.asarray(levels)[1:])
        layers[0] = np.min(function_map)

        bbx, bby = 0.5 * (bbox[0] + bbox[2]), 0.5 * (bbox[1] + bbox[3])
        do_warn = False
        polygons = []

        for i, item in enumerate(contours.collections):
            edge_paths = []
            for path in item.get_paths():
                codes = path.codes
                vertices = path.vertices
                if codes[-1] != Path.CLOSEPOLY and not np.allclose(
                    vertices[0], vertices[-1]
                ):
                    edge_paths.append(vertices)
                    continue
                x, y = zip(*vertices)
                x, y = np.asarray(x), np.asarray(y)
                xc, yc = np.mean(x), np.mean(y)
                x, y = x - xc, y - yc
                area = get_polygon_area(x, y)
                if area < 5.0e-6:
                    continue
                vertices = list(zip(x, y))
                shrink_vertices = grow_polygons([vertices], grow=-0.001, accuracy=0.001)
                x, y = zip(*shrink_vertices[0])
                x, y = np.asarray(x), np.asarray(y)
                x, y = x + xc, y + yc
                n = np.mean(function(x=x, y=y, **fkwargs))
                level = layers[i] if n > levels[i] else layers[max(i - 1, 0)]
                polygons.append(PolygonData(path.vertices, level, area))
            else:
                if len(edge_paths) == 0:
                    continue
                do_warn = True
                start = []
                for vertices in edge_paths:
                    x, y = zip(*vertices)
                    start.append((vertices[0][0] - bbx, vertices[0][1] - bby))
                full_vertices = []
                for j in _argsort(start, _clockwiseangle_and_distance)[::-1]:
                    full_vertices = full_vertices + list(edge_paths[j])
                x, y = zip(*full_vertices)
                x, y = np.asarray(x), np.asarray(y)
                x, y = x - bbx, y - bby
                area = get_polygon_area(x, y)
                shifted_vertices = list(zip(x, y))
                shrink_vertices = grow_polygons(
                    [shifted_vertices], grow=-0.001, accuracy=0.001
                )
                x, y = zip(*shrink_vertices[0])
                x, y = np.asarray(x), np.asarray(y)
                x, y = x + bbx, y + bby
                n = np.mean(function(x=x, y=y, **fkwargs))
                level = layers[i] if n > levels[i] else layers[max(i - 1, 0)]
                polygons.append(PolygonData(full_vertices, level, area))

        if do_warn:
            nazca.logger.warning(
                "Not all the contour lines are inside the bbox for Geometry.from_function. Correctness cannot be guaranteed, plese verify index profile manually"
            )

        polygons.sort(key=lambda tup: tup.area)

        with cls(background_index=levels[0], bbox=bbox) as geo:
            for pol in polygons[::-1]:
                Polygon(points=pol.vertices, material_index=pol.level).put()

        return geo

    def put(
        self,
        y=0.0,
        z=0.0,
        a=0.0,
        scale=1.0,
        flip=False,
        flop=False,
    ):
        """Puts the geometry into the Geometry class in the location provided.

        Args:
            y (float): y translation in um (default = 0.0).
            z (float): z translation in um (default = 0.0).
            a (float): a translation in deg (default = 0.0).
            scale (float): scaling factor (default = 1.0).
            flip (bool): flip x coordinate x -> -x (default = False).
            flop (bool): flip y coordinate y -> -y (default = False).

        Returns:
            self
        """

        for shape in self.shapes:
            _instance = Instance(
                shape.shape,
                y=y + shape.y,
                z=z + shape.z,
                a=a + shape.a,
                scale=scale,
                flip=flip,
                flop=flop,
            )

            Geometry.active_geometries[-1].add_shape(instance=_instance)

        return self


class Shape:
    """Geometrical shape class"""

    def __init__(self, material_index):
        if callable(material_index):
            self.material_index = filter_eval(material_index)
        else:
            self.material_index = material_index

    def put(
        self,
        y=0.0,
        z=0.0,
        a=0.0,
        scale=1.0,
        flip=False,
        flop=False,
    ):
        """Puts the shape into the Geometry class in the location provided.

        Args:
            y (float): y translation in um (default = 0.0).
            z (float): z translation in um (default = 0.0).
            a (float): a translation in deg (default = 0.0).
            scale (float): scaling factor (default = 1.0).
            flip (bool): flip x coordinate x -> -x (default = False).
            flop (bool): flip y coordinate y -> -y (default = False).

        Returns:
            self
        """

        _instance = Instance(
            self,
            y=y,
            z=z,
            a=a,
            scale=scale,
            flip=flip,
            flop=flop,
        )

        Geometry.active_geometries[-1].add_shape(instance=_instance)

        return _instance

    def get_points(self):
        raise NotImplementedError(
            "You are trying to get the points from the parent class. Call it from the child one"
            "instead\n"
        )

    def get_material_index(self, **kwargs):
        """Return material index.

        If material index is a function, kwargs can be provided

        Args:
            **kwargs: dictionary of arguments to pass to material index function

        Returns:
            float | complex: the material index of the shape
        """
        if callable(self.material_index):
            return self.material_index(**kwargs)
        else:
            return self.material_index


class Instance:
    def __init__(
        self,
        shape,
        y=0.0,
        z=0.0,
        a=0.0,
        scale=1.0,
        flip=False,
        flop=False,
    ):
        """Initializes the class.

        Args:
            shape (Shape): Shape that generates the instance.
            y (float): y translation in um (default = 0.0).
            z (float): z translation in um (default = 0.0).
            a (float): a translation in deg (default = 0.0).
            scale (float): scaling factor (default = 1.0).
            flip (bool): flip x coordinate x -> -x (default = False).
            flop (bool): flip y coordinate y -> -y (default = False).

        Returns:
            None
        """
        self.shape = shape
        self.y = y
        self.z = z
        self.a = a
        self.scale = scale
        self.flip = flip
        self.flop = flop

    def get_points(self):
        """Calculates the points of the instance.

        Returns:
            list[tuple]: The modified points of the instance.
        """
        return transform_polygon(
            self.shape.get_points(),
            dx=self.y,
            dy=self.z,
            da=self.a,
            scale=self.scale,
            flipx=self.flip,
            flipy=self.flop,
        )

    @property
    def material_index(self):
        """Returns the material index of the corresponding shape"""
        return self.shape.material_index

    def get_material_index(self, **kwargs):
        """Return material index of the parent shape"""
        return self.shape.get_material_index(**kwargs)

    def get_edges_representation(self, **kwargs):
        """Return list of edges of the shape

        Returns:
            list[tuple]: list of tuples of edges ((x1,y1),(x2,y2))
        """
        points = self.get_points()
        return list(zip(points, points[1:] + points[:1]))

    def get_3d_points(
        self,
        xyb,
        xyt,
        local_width,
        xs_width,
    ):
        """Calculates the 3D points of the cross-section given the spine points.

        The routine already takes into account the requested polygon transformations specified in the put, so that
        the shapes are placed in the correct way.

        Args:
            xyb (tuple): x and y coordinates of the point at the bottom of the spine.
            xyt (tuple): x and y coordinates of the point at the top of the spine.
            local_width (float): Local width of the waveguide in [um].
            xs_width (float): Width of the cross-section in [um].

        Returns:
            list of tuple: The 3D coordinates of each point of the cross-section.
        """
        _pts = self.get_points()

        _pts_3d = []

        if xyb[1] == xyt[1]:  # Same y => theta=pi/2
            angle = np.pi / 2
        else:
            angle = np.arctan2(xyt[0] - xyb[0], xyb[1] - xyt[1])

        _y_mid = (xyb[1] + xyt[1]) / 2
        _y_1 = xyb[1] - _y_mid
        _y_2 = xyt[1] - _y_mid
        _x_mid = (xyb[0] + xyt[0]) / 2
        _x_1 = xyb[0] - _x_mid
        _x_2 = xyt[0] - _x_mid

        for pt in _pts:
            _pts_3d.append(
                (
                    _f_x(pt[0], _x_1, _x_2, local_width, xs_width, angle) + _x_mid,
                    _f_y(pt[0], _y_1, _y_2, local_width, xs_width, angle) + _y_mid,
                    pt[1],
                )
            )

        return _pts_3d


class Circle(Shape):
    def __init__(self, radius, N=50, material_index=None):
        """Initializes the class.

        Args:
            radius (float): Radius of the circle.
            N (int): Number of points used to represent the circle.
            material_index (float or complex): Refractive index of the circle.

        Returns:
            None
        """
        super().__init__(material_index=material_index)
        self.radius = radius
        self.N = N

    def get_points(self):
        """Calculates the points of the circle.

        Returns:
            list[tuple]: The points of the circle.
        """
        angles = np.linspace(0, 360, self.N + 1)
        return [
            (
                self.radius * m.cos(m.radians(angle)),
                self.radius * m.sin(m.radians(angle)),
            )
            for angle in angles
        ]


class Polygon(Shape):
    def __init__(self, points=None, material_index=None):
        """Initializes the class.

        Args:
            points (list of tuple): List containing the coordinates of the points of the polygon.
            material_index (float or complex): Refractive index

        Returns:
            None
        """
        super().__init__(material_index=material_index)
        self.points = points

    def get_points(self):
        """Returns the points of the polygon.

        Returns:
            list[tuple]: The points of the polygon.
        """
        return self.points


class Trapezoid(Polygon):
    def __init__(self, top_width, height, base_angle=90.0, material_index=None):
        """Initializes the class.

        Args:
            top_width (float): Width of the upper side.
            height (float): Height of the trapezoid.
            base_angle (float): Angle in degrees between the bottom and top row. 90 degrees for a rectangle.
            material_index (float or complex): Refractive index

        Returns:
            None
        """
        assert 0 <= base_angle <= 90

        super().__init__(material_index=material_index)
        self.top_width = top_width
        self.height = height
        self.base_angle = base_angle

        add = height * m.tan(m.radians(90 - self.base_angle))
        self.base_width = self.top_width + 2 * add

        self.points = [
            (-self.base_width / 2, -self.height / 2),
            (self.base_width / 2, -self.height / 2),
            (self.top_width / 2, self.height / 2),
            (-self.top_width / 2, self.height / 2),
        ]


class IndexModel:
    def __init__(
        self, xs_function, solver, ridge_solver=None, bbox=None, model_name=None
    ):
        """Initializes the class.

        Args:
            xs_function (callable): The function that returns the cross-section geometry.
            solver (any): The desired solver to solve the cross-section.
            ridge_solver (any): The solver that extracts the ridge index. Default is None.
            bbox (list): list of float of length 4. If provided, overwrites the xs native bbox
                Meaning is (x bottom left, y bottom left, x top right, y top right)

        Returns:
            None
        """
        if not callable(xs_function):
            raise TypeError("xs_function is not callable. Please provide a function.\n")

        self.xs = None
        self.xs_function = filter_eval(xs_function)
        self.solver = solver
        self.ridge_solver = solver if ridge_solver is None else ridge_solver
        self.delta_wl_grp = 0.001
        self.bbox = bbox
        self.model_name = model_name
        self._create_metadata(model_name=model_name)

    def _create_metadata(self, model_name=None):
        xs_original = self.xs_function.__wrapped__
        self.meta_data = {
            "Model Name": model_name,
            "Metadata_creation": "Auto",
            "Author": subprocess.check_output(
                "git config user.name",
                shell=True,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            ).replace("\n", ""),
            "Author Email": subprocess.check_output(
                "git config user.email",
                shell=True,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            ).replace("\n", ""),
            "Nazca version": nazca.__version__,
            "XS function": _get_function_data(xs_original),
        }
        xs = self.xs_function()
        self.meta_data["Material functions"] = [
            _get_function_data(mat)
            if callable(mat)
            else {"Constant material index": mat}
            for mat in xs.get_material_functions()
        ]
        self.meta_data["Solver"] = {
            "Type": f"{self.solver.__module__}.{self.solver.__class__.__name__}",
            "file": os.path.split(sys.modules[self.solver.__module__].__file__)[1],
            "file path": os.path.split(sys.modules[self.solver.__module__].__file__)[0],
            "Git Repository": _get_function_git(self.solver),
        }
        try:
            self.meta_data["Solver"]["Data"] = self.solver.get_metadata()
        except AttributeError:
            pass
        return self.meta_data

    def log_metadata(self, filename=None):
        """Logs the metadata of the index model

        Args:
            filename (str): Name of the log file

        Returns:
            None
        """
        filename = (
            f"{self.model_name}.yaml" if self.model_name is not None else filename
        )
        filename = filename if filename is not None else "index_model.yaml"
        with open(filename, "w") as f:
            yaml.dump(self.meta_data, f)

    def process_xsection(self, **kwargs):
        """Calculates the Geometry object and send it to the solvers

        The Geometry object is created calling self.xs_funciton wiht the given keyword arguments

        Args:
            **kwargs: Keyword arguments for xs_function

        Returns:
            nazca.xsection_geometry.Geometry: The obtained cross section.
        """
        self.xs = self.xs_function(**kwargs)
        if self.bbox is not None:
            self.xs.bbox = self.bbox
        self.solver.xs = self.xs
        try:
            self.ridge_solver.xs = self.xs
        except AttributeError:
            pass
        return self.xs

    @lru_cache(128)
    def Neff(self, **kwargs):
        """Effective index function

        Args:
            kwargs: Keyword arguments

        Returns:
            float | complex: The effective index of the waveguide.
        """
        self.process_xsection(**kwargs)
        return self.solver.Neff(**kwargs)

    @lru_cache(128)
    def Nridge(self, **kwargs):
        """Returns the Nridge for the index

        Args:
            **kwargs: Keyword arguments

        Returns:
            float | complex: Ridge effective index
        """
        self.process_xsection(**kwargs)
        return self.ridge_solver.Nridge(**kwargs)

    @lru_cache(128)
    def Nbg(self, **kwargs):
        """Returns the background index of the cross section

        Args:
            **kwargs: Keyword arguments

        Returns:
            float | complex: Background index
        """
        xs = self.process_xsection(**kwargs)
        return xs.get_background_index(**kwargs)

    def Ngrp(self, wl=None, **kwargs):
        """Calculates the group index by central difference

        The step used can be changed setting the delta_wl_grp attribute

        Args:
            wl (float): Wavelenght
            **kwargs: Other keyword arguments

        Returns:
            float | complex: The group index
        """
        d = self.delta_wl_grp
        dn_dwl = (
            0.5 * (self.Neff(wl=wl + d, **kwargs) - self.Neff(wl=wl - d, **kwargs)) / d
        )
        return self.Neff(wl=wl, **kwargs) - wl * dn_dwl

    def Loss(self, **kwargs):
        """Bending loss function

        Args:
            kwargs: Keyword arguments

        Returns:
             float: The bending loss
        """
        self.process_xsection(**kwargs)
        return self.solver.Loss(**kwargs)

    def te_fraction(self, **kwargs):
        """te fraction function

        Args:
            kwargs: Keyword arguments

        Returns:
             float: The te fraction of the mode
        """
        self.process_xsection(**kwargs)
        return self.solver.te_fraction(**kwargs)

    def field1D(self, **kwargs):
        """1D field function.

        Args:
            kwargs: Keyword arguments

        Returns:
            any: The 1 dimensional field.
        """
        self.process_xsection(**kwargs)
        return self.solver.field1D(**kwargs)

    def field2D(self, **kwargs):
        """2D field function.

        Args:
            kwargs: Keyword arguments

        Returns:
            any: The 2 dimensional field.
        """
        self.process_xsection(**kwargs)
        return self.solver.field2D(**kwargs)

    def index2D(self, **kwargs):
        """2D index function.

        Args:
            kwargs: Keyword arguments

        Returns:
            any: The 2 dimensional refractive index.
        """
        self.process_xsection(**kwargs)
        return self.solver.index2D(**kwargs)

    def modes(self, **kwargs):
        """Returns the number of supported modes

        Args:
            **kwargs: Keyword arguments

        Returns:
            float: number of supported modes
        """
        kwargs.pop("mode", None)
        self.process_xsection(**kwargs)
        return self.solver.modes(**kwargs)

    def show_edges(
        self,
        ax=None,
        color="black",
        linewidth=0.2,
        remove_unnecessary_edges=True,
        **kwargs,
    ):
        """Show the edges between shapes in the cross-section

        Args:
            ax (matplotlib.Axes): Axes to use for the plot. If None, plot is shown
            color (str): Color for the edges. Default is black.
            linewidth (float): Width of the contour line
            remove_unnecessary_edges (Bool): If True (default) edges between same index are removed.
            **kwargs: Other arguments to be passed to xs_function. The arguments for the index functions are not needed in this function.

        Returns:
            matplotlib.Figure: The figure containing the axes.
            matplotlib.Axes: The axis onto which the geometry is plotted.
        """
        xs = self.process_xsection(**kwargs)
        return xs.show_edges(
            ax=ax,
            color=color,
            linewidth=linewidth,
            remove_unnecessary_edges=remove_unnecessary_edges,
        )

    def L100(self, wl, distance, mode0=0, mode1=1, **kwargs):
        """Returns the 100% coupling length between

        Args:
            wl (float): Wavelength in um.
            distance (float): Ditance between the coupled waveguides
            mode0 (int): first mode to use for the calculation
            mode1 (int): second mode to use for the calculation
            **kwargs: Any other keyword argument.

        Returns:
            float: 100% coupling length of the system
        """
        mode = kwargs.pop("mode", None)
        if mode is not None:
            nazca.logger.error(
                "Called {type(self).__name__}.L100 with keyword argument 'mode'. This is not allowed. See 'mode0' and 'mode1' instead"
            )
        n1 = self.Neff(wl=wl, distance=distance, mode=mode0, **kwargs)
        n2 = self.Neff(wl=wl, distance=distance, mode=mode1, **kwargs)
        return wl / (2.0 * abs(n1 - n2))


if __name__ == "__main__":
    with Geometry() as xs_g:
        Trapezoid(1.0, 1.0, 90.0, material_index=4.0).put()

    xs_g.show()

    trapezoid = Trapezoid(top_width=1, height=1, base_angle=82, material_index=1.9)

    # x = np.linspace(-10.0, 10.0, 501)
    # X, Y = np.meshgrid(x, x, indexing='ij')
    # a = xs_g.get_index(X, Y)
    # plt.contourf(X, Y, a)
