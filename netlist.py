#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------
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
#
# @author: Ronald Broeke (c) 2016-2020, 2011 (c) Xaveer Leijtens,
# @email: ronald.broeke@brightphotonics.eu
#
"""
Nazca classes to construct Nodes, Pointers and Cells,
and from those the graphs/netlists of the design.
"""

import os
from itertools import count
from collections import defaultdict, OrderedDict, namedtuple
import copy as COPY
from math import sin, cos, acos, pi, atan2, hypot

from scipy.spatial import ConvexHull
from numpy.linalg import inv
from numpy import dot
import numpy as np
from pprint import pprint
from .clipper import grow_polygons
import yaml

from . import cfg
from nazca.logging import logger
from .mask_layers import get_layer, get_xsection
from .geometries import ring, transform
from . import gds_import as gdsimp
from . import gds_base as gb
import nazca.trace as trace
import nazca.bb_util as bbu
import nazca.font as font
import nazca.util as util


Et = count(0) # counter for default cell names:
elmlist = []
Pt = count(0)

# message counters:
cfg.ND_msgcnt = 0  # counter for cummulative errors
cfg.NL_msgcnt = 0  # counter for netlist errors
cfg.IC_msgcnt = 0  # counter for interconnect messages
cfg.drc_msgcnt = 0  # count DRC messages.drc (None raises all)

cfg.ND_raise = False
cfg.NL_raise = False
cfg.IC_raise = False

# gds loading flags:
cfg.gdsload = False  # flag that contains the gds file name if a file is loading. False otherwise
cfg.gdsloadstore = ""  # store name of gds if a layer mapping conflict occurs file while loading a gds to only message a warning once.


def to_mat(chain=0, x=0, y=0, a=0, inverse=False):
    """Transform vector into matrix."""
    a = (a+chain*180.0) / 180. * pi
    c, s = cos(a), sin(a)
    if not inverse:
        return np.array([[ c, s, 0.0],
                         [-s, c, 0.0],
                         [x, y, 1.0]])
    else:
        return np.array([[c, -s, 0.0],
                         [s,  c, 0.0],
                         [-x*c-y*s, x*s-y*c, 1.0]])


def inverse(mat):
    """Return inverse matrix of <mat>"""
    c, s = mat[0, 0], mat[0, 1]
    x, y = mat[2, 0], mat[2, 1]
    return np.array([[c, -s, 0.0],
                     [s,  c, 0.0],
                     [-x*c-y*s, x*s-y*c, 1.0]])
    # return inv(mat) # ~4x slower


def set_db_unit():
    pass


def set_db_user():
    pass


def drc_instance_angle(rule, cell=None):
    """Set design rule for the angle of an instantiated cell.
    """
    if cell is None:
        cell = cfg.cells[-1]
    cfg.drc_instance_angle['angle'][cell.cell_basename] = rule


def pin2pin_drc(on=True):
    """Switch pin2pin DRC on.

    NOTE: the drc error maybe in a connection that is not exported to (gds) output
        and, therefor, not show up in that specific layout view.

    Args:
        on (bool): raise exception for trace back if True (default=True)
        num (int): if provided, raise exception on drc number <num> (default=None)
            if num is set then raise=True unless explicitly set to False

    Returns:
        None
    """
    if on:
        if not cfg.drc:
            cfg.drc = True
            logger.info("pin2pin drc: {}".format(cfg.drc))
    else:
        cfg.drc = False
        logger.info("pin2pin drc: {}".format(cfg.drc))

def pin2pin_drc_on(exception=None, num=None):
    pin2pin_drc(on=True)
def pin2pin_drc_off():
    pin2pin_drc(on=False)


def pin2pin_drc_raise(num=0):
    if num > 0:
        logger.info("Set to raise pin2pin drc on message: {}".format(num))
    if not cfg.drc:
        pin2pin_drc(on=True)
    cfg.drc_raisenum = num
    cfg.drc_raise = True
raise_pin2pin_drc = pin2pin_drc_raise


def netlist_raise(num=0, on=True):
    """Switch on or off exceptions for netlist errors.

    The traceback can be used to find the cause of the netlist error.

    Args:
        num (int): if provided, raise exception on drc number <num> (default=None)
            if num is set then raise=True unless explicitly set to False
        on (bool): raise exception for trace back if True (default=True)

    Returns:
        None
    """
    if on:
        cfg.NL_raisenum = num
        cfg.NL_raise = True
    else:
        cfg.NL_raise = False

raise_netlist = netlist_raise


def interconnect_raise(num=0, on=True):
    """Switch on or off exceptions for interconnect warnings and errors.

    The traceback can be used to find the cause of the netlist error.

    Args:
        num (int): if provided, raise exception on drc number <num> (default=None)
            if num is set then raise=True unless explicitly set to False
        on (bool): raise exception for trace back if True (default=True)

    Returns:
        None
    """
    if on:
        cfg.IC_raisenum = num
        cfg.IC_raise = True
    else:
        cfg.IC_raise = False

raise_interconnect = interconnect_raise


def nazca_raise(num=0, on=True):
    """Raise an exception for any design related Nazca error.

    Args:
        num (int): if provided, raise exception on drc number <num> (default=None)
            if num is set then raise=True unless explicitly set to False
        on (bool): raise exception for trace back if True (default=True)

    Returns:
        None
    """
    if on:
        cfg.ND_raisenum = num
        cfg.ND_raise = True
    else:
       cfg.ND_raise = False
raise_nazca = nazca_raise


# =============================================================================
# message loggers to be used in combination with raising Exceptions.
# =============================================================================
def netlist_logger(msg='', level='info'):
    """Raise netlist exception if conditions apply.

    Args:
        msg (str): log message
        level (str): log level: 'debug', 'info' (default), 'error', or 'exception'.

    Returns:
        None
    """
    fullmsg = f"NL-{cfg.NL_msgcnt}/ND-{cfg.ND_msgcnt}: {msg}"
    if cfg.NL_raise and cfg.NL_msgcnt == cfg.NL_raisenum:
        raise Exception(fullmsg)
        # TIP: use the traceback (2 steps up) to find the cause of the exception.

    elif cfg.ND_raise:
        main_logger(msg=msg, level=level)

    cfg.NL_msgcnt += 1
    cfg.ND_msgcnt += 1

    if level == 'warning':
        logger.warning(fullmsg)
    elif level == 'error':
       logger.error(fullmsg)
    elif level == 'info':
       logger.info(fullmsg)
    else:
        raise Exception(f"Unknown interconnect_logger level given: '{level}'")
        # TIP: use the traceback (2 steps up) to find the cause of the exception.


def interconnect_logger(msg='', level='info'):
    """Log or raise an interconnect exception if conditions apply.

    Args:
        msg (str): log message
        level (str): log level: 'debug', 'info' (default), 'error', or 'exception'.

    Raises:
        Exception if the "exception number" is set and matching this log entry.
        Exception if unknown level is provided

    Returns:
        None
    """
    cfg.IC_msgcnt += 1
    cfg.ND_msgcnt += 1
    fullmsg = f"IC-{cfg.IC_msgcnt}/ND-{cfg.ND_msgcnt}: {msg}"
    if level == 'warning':
        logger.warning(fullmsg)
    elif level == 'error':
       logger.error(fullmsg)
    elif level == 'info':
       logger.info(fullmsg)
    else:
        raise Exception(f"Unknown interconnect_logger level given: '{level}'")
            # TIP: use the traceback (2 steps up) to find the cause of the exception.
    if cfg.IC_raise:
        if cfg.IC_msgcnt == cfg.IC_raisenum:
            raise Exception(fullmsg)
            # TIP: use the traceback (2 steps up) to find the cause of the exception.
    elif cfg.ND_raise:
        main_logger(msg=msg, level=level, _redirect=True)


# TODO: check out logging.LoggerAdapter and logging.Filter
def main_logger(msg, level='info', _redirect=False):
    """Log or raise an exception if conditions apply.

    This message handler will take care of loggin any error or warning.
    It will raise an exception if the error counter = errnum
    (if errnum has been set) for debugging/trace back purposes.

    Arg:
        msg (str): log message
        level (str): log level: 'debug', 'info' (default), 'error', or 'exception'.
        _redirect (bool): internal variable for counting messages.

    Returns:
        None
    """
    if not _redirect:
        cfg.ND_msgcnt += 1

    fullmsg = f"ND-{cfg.ND_msgcnt}: {msg}"
    if cfg.ND_raise and cfg.ND_msgcnt == cfg.ND_raisenum:
        raise Exception(fullmsg)
        # TIP: use the traceback (2 steps up) to find the cause of the exception.

    if level == 'warning':
        logger.warning(fullmsg)
    elif level == 'error':
       logger.error(fullmsg)
    elif level == 'info':
       logger.info(fullmsg)
    else:
        raise Exception(f"Unknown interconnect_logger level given: '{level}'")
        # TIP: use the traceback (2 steps up) to find the cause of the exception.


def varinfo(min, default, max, type, unit, doc):
    """Set parameter value info.

    Experimental function.

    Returns:
        dict: var info
    """
    return {
        'min': min,
        'default': default,
        'max': max,
        'type': type,
        'unit': unit,
        'doc': doc,
    }


class Pointer():
    """A pointer with state information.

    Note: For normal Nazca usage avoid using this class directly and use the methods
    provided in the Node class instead.

    The pointer has the positional information of a node.
    The positional information is kept in a matrix that holds the homogeneous
    coordinates. Structures are drawn in local coordinates and are easily
    converted to the global coordinate system by a matrix multiplication. By
    using homogeneous coordinates, translations are also described by a
    matrix multiplication.
    """

    def __init__(self, x=0.0, y=0.0, a=0.0):
        """Constructor.

        Args:
            x (float): x coordinate in um
            y (float): y coordinate in um
            a (float): angle coordinate un degrees
        """
        self.flipstate = False
        self.mat = None
        self.chain = 0
        self.goto(x, y, a)


    def __str__(self):
        return "{0}:\t(x = {2:.2f}, y = {3:.2f}, a = {1:.2f}°)".\
            format(self.__class__.__name__, self.a, *self.get_xy())


    def goto(self, x, y, a=0.0):  # Absolute
        """Move pointer to a coordinate relative to the org.

        Args:
            x (float): x-coordinate of the pointer position [µm]
            y (float): y-coordinate of the pointer position [µm]
            a (float): angle of the pointer [degrees]

        Returns:
            None
        """
        a = pi*a/180
        self.mat = np.array([
            [cos(a), sin(a), 0],
            [-sin(a), cos(a), 0],
            [x, y, 1]
            ])
        self.flipstate = False


    def move(self, x=0, y=0, a=0.0):  # Relative
        """Move pointer relative to current location.

        Args:
            x (float): x-coordinate of the pointer position [µm]
            y (float): y-coordinate of the pointer position [µm]
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        a = pi*a/180
        mat = np.array([[cos(a), sin(a), 0.0],
                        [-sin(a), cos(a), 0.0],
                        [x, y, 1.0]])
        self.mat = dot(mat, self.mat)
        return self


    def set_a(self, a=0.0):
        """Set angle absolute. Keep x, y position

        Args:
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        a = pi*a/180
        x = self.x
        y = self.y
        self.mat = np.array([
            [cos(a), sin(a), 0.0],
            [-sin(a), cos(a), 0.0],
            [x, y, 1.0]
        ])
        # self.mat = dot(mat, self.mat)
        return self


    def move_ptr(self, ptr):  # Relative
        """Move pointer relative by a pointer 'ptr' w.r.t. current location.

        Args:
            ptr (Pointer): move by value in pointer

        Returns:
            Pointer: self
        """
        x, y, a = ptr.get_xya()
        a = pi*a/180
        mat = np.array([
            [cos(a), sin(a), 0.0],
            [-sin(a), cos(a), 0.0],
            [x, y, 1.0]
        ])
        self.mat = dot(mat, self.mat)
        return self


    def multiply_ptr(self, ptr):  # Relative
        """Multiply the pointer by the matrix in <ptr>.

        Args:
            ptr (Pointer): multiply by <ptr>

        Returns:
            Pointer: self
        """
        self.mat = dot(self.mat, ptr.mat)
        return self


    def inv(self):  # Relative
        """Inverse the matrix in the pointer. Returns the pointer.

        Returns:
            Pointer: self
        """
        self.mat = inv(self.mat)
        return self


    def trans(self, mat):  # Relative
        """Translate pointer by matrix <mat>.

        Return:
            Pointer: translated by matrix <mat>
        """
        return dot(mat, self.mat)


    def chain180(self, chain):
        """Set the chain property of the pointer.

        Args:
            chain (bool): True for chain connecting pins.

        Returns:
            None
        """
        self.chain = chain


    def set_mat(self, t):  # Relative
        """Fill pointer matrix based on vector <t> = (x, y, a).

        Args:
            t (tuple): (x, y, a)
        Returns:
           None
        """
        self.mat = t


    def skip(self, s):  # move forward (negative: backward)
        """Translate a pointer in the direction it is pointing in.

        Args:
            s (float): move pointer along its direction by <s>

        Returns:
            Pointer: self, moved
        """
        self.move(s, 0, 0)
        return self


    def offset(self, s):  # offset positive/negative to left/right.
        """Translate a pointer perpendicular to the direction it is pointing in.

        Note that in a cartesian coordinate system while pointing/looking
        along the positive x-axis, a positive offset is in the left direction,
        i.e. along the positive y-axis.

        Args:
            s (float): move pointer along perpendicular direction by <s>

        Returns:
            Pointer: self, moved pointer
        """
        self.move(0, s, 0)
        return self


    def flip(self):
        """Flip (mirror) the state of the pointer.

        Returns:
            Pointer: self, flipped
        """
        self.mat = dot(self.mat, [[1, 0, 0], [0, -1, 0], [0, 0, 1]])
        self.flipstate = not self.flipstate
        return self


    def rotate(self, a):
        """Rotate pointer by angle a.

        Args:
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self, rotated
        """
        self.move(0, 0, a)
        return self


    @property
    def x(self):
        """Return position x of the pointer."""
        return self.mat[2, 0]


    @property
    def y(self):
        """Return the position y of the pointer."""
        return self.mat[2, 1]


    @property
    def a(self):
        """Return angle a of the pointer."""
        a = min(1, max(self.mat[0, 0], -1))  # clip to [-1, 1]
        if self.mat[0, 1] >= 0:
            return acos(a) * 180/pi
        else:
            return 360 - acos(a) * 180/pi


    def get_xy(self, x=0, y=0):
        """Return the (x, y) position of the pointer."""
        return tuple(dot([x, y, 1], self.mat))[0:2]


    def get_xya(self):
        """Return the pointer position as (x, y, a)."""
        return (self.x, self.y, self.a)


    def get_xy_a(self):
        """Return the pointer position as (x, y, a)."""
        return ((self.x, self.y), self.a)


    def xya(self):
        """Return the pointer position as (x, y, a)."""
        return (self.x, self.y, self.a)


    def fxya(self, digits=3):
        """Get formatted pointer position (x, y, a) of the Node.

        Returns:
            str: formatted pointer position (x, y, a)
        """
        return "({0:0.{prec}f}, {1:0.{prec}f}, {2:0.{prec}f})".format(self.x, self.y, self.a, prec=digits)


    def xy(self, x=0, y=0):
        """Return the (x, y) position of the pointer."""
        return tuple(dot([x, y, 1], self.mat))[0:2]


    def copy(self, flip=False):
        """Copy pointer.

        Returns:
            Pointer: copy of the Pointer
        """
        # One level deep copy. The state information is still copied by
        # reference only. This is probably the desired way of working.
        copy = COPY.copy(self)
        if flip:
            copy.y = -copy.y
            copy.a = -copy.a
        return copy


class smatrix:
    """S-matrix tbd."""
    #TE-TE
    #TM-TM
    #TE-TM
    #TM-TE
    #phase
    #amplitude
    def __init__(self, tete=(0, 0), tmtm=(0, 0), tetm=(0, 0), tmte=(0, 0)):
        self.tete = tete
        self.tmtm = tmtm
        self.tetm = tetm
        self.tmte = tmte


class Node():
    """Node class for creating node objects.

    Note that pins in the Nazca design terminology are actually Node objects.

    The netlist of Nodes and connections between the nodes (edges) construct
    the photonic circuit and/or layout.
    """

    def __init__(self, name=None):
        """Contruct a Node.

        Args:
            name (str): optional Node name (default is an integer counter)

        Returns:
            None
        """
        #print('______________________new node', name)
        self.id = next(Pt)
        self.nb_geo = []  #nearest neighbors geometry
        self.nb_opt = []  #nearest neighbors, optical connection (in Smatrix)
        self.nb_ele = []  #nearest neighbors, electronic connection
        self.nb_cnode = []  #cell tree, only for used for cnodes
        self.cnode = None  # cnode of the node. A cnode is its own cnode
        self.parent_cnode = None
        self.cell = None  # TODO: in use? Cell object of the cnode
        self.pointer = None  # Pointer of the Node. Contains the node position after solving.
        self.xs = None  # String name reference to the xs of the node.
        self.width = None  # width of the node
        self.radius = None  # radius of the curvature at the pin.
        self.instance = False  # store if node is in a cell (False) or instance (True)
        self.show = False  # show pin in layout.
        self.remark = None  # add info to the pin.
        self.io = None   # indicate of a port is an explicitly in or out in case of directional ports
        self.type = None  #

        if name is None:
            self.name = str(self.id)
        else:
            self.name = str(name)


    def __repr__(self):
        coor = "({t[0]:0.3f}, {t[1]:0.3f}, {t[2]:0.3f})".format(t=self.xya())
        if self.cnode.instance is not None:
            cell = "instance of cell '{}'".format(self.cnode.cell.cell_name)
        else:
            cell = "cell '{}'".format(self.cnode.cell.cell_name)
        return "<Node(id={}, name='{}') object in {}, xs='{}', width={}, xya={}, type={}, io={}, remark='{}'>".\
            format(self.id, self.name, cell, self.xs, self.width, coor, self.type, self.io, self.remark)


    def copy(self, inplace=False):
        """Copy a Node object into the active cell or inplace.

        Args:
            inplace (bool): Copy the node into the same cell if False (True) or the active cell (default).

        Returns:
            Node: copy of the node with the same attributes as the original

        Example::

            import nazca as nd

            n1 = nd.Node()
            n2 = n1.copy()
        """
        # do NOT use deepcopy: extremely slow.
        node = Node()
        node.pointer = self.pointer.copy()
        node.xs = self.xs
        node.width = self.width
        node.instance = self.instance
        node.type = self.type
        node.io = self.io
        node.remark = self.remark

        if inplace:
            cnode = self.cnode
        else:
            cnode = cfg.cells[-1].cnode
        if self.cnode is not cnode and self.cnode.parent_cnode is not cnode:
            raise Exception('It is not allowed to create a node outside the scope of max one level deep:'\
                ' You cannot connect to a cell that is not placed inside the cell you are creating.')

        node.cnode = cfg.cells[-1].cnode #place node in active cell.
        node.cnode.parent_cnode = cfg.cells[-1].parent_cnode
        return node


    @property
    def chain(self):
        """Attribute to indicate if a Node is a chain connector."""
        return self.pointer.chain


    def goto(self, x=0, y=0, a=0):
        """Create a new Node at position (<x>, <y>, <a>) relative to the cell origin.

        Args:
            x (float): x position in um
            y (float): y position in um
            a (float): angle a in degrees

        Returns:
            Node: a new Node at (<x>, <y>, <a>) w.r.t. the cell's origin
        """
        node = Node()
        node.pointer = COPY.copy(self.pointer)
        node.xs = self.xs
        node.width = self.width
        node.type = self.type
        
        #node.pointer.chain = 0
        node.cnode = cfg.cells[-1].cnode
        node.parent_cnode = cfg.cells[-1].parent_cnode
        connect_geo(cfg.cells[-1].cnode, node, (x, y, a), 1, 0)
        #print('cfg.cells[-1].cnode.cell.cell_name:', cfg.cells[-1].cnode.cell.cell_name)
        cfg.cp = node
        return node


    def xya(self):
        """Get pointer position of the Node.

        Returns:
            tuple: pointer position (x, y, a) with respect to the cell origin
        """
        return self.pointer.get_xya()


    def fxya(self, digits=3):
        """Get formatted pointer position (x, y, a) of the Node.

        Returns:
            str: formatted pointer position (x, y, a) with respect to the cell origin
        """
        return self.pointer.fxya(digits=digits)


    def fxy(self, digits=3):
        """Get formatted pointer position (x, y) of the Node.

        Returns:
            str: formatted pointer position (x, y) with respect to the cell origin
        """
        p = self.pointer.get_xya()
        return "({:.3f}, {:.3f})".format(p[0], p[1])


    @property
    def x(self):
        """Get pointer x-position of the Node.

        Returns:
            tuple: pointer x-position with respect to the cell origin
        """
        return self.pointer.x


    @property
    def y(self):
        """Get pointer y-position of the Node.

        Returns:
            tuple: pointer y-position with respect to the cell origin
        """
        return self.pointer.y

    @property
    def a(self):
        """Get pointer a-position of the Node.

        Returns:
            tuple: pointer a-position with respect to the cell origin
        """
        return self.pointer.a

    def angle(self, angle):
        """Get pointer a-position of the Node.

        Returns:
            tuple: pointer a-position with respect to the cell origin
        """
        self.pointer.set_a(angle)
        return self



#==============================================================================
#     Create pointer operations on Node level for syntax convenience
#==============================================================================
    def move(self, x=0, y=0, a=0, xs=0, width=0, radius=0, drc=True):
        """Create a new Node translated by (<x>, <y>, <a>) with respect to self.

        Also 'width' and 'xs' attributes of the can be 'moved'. This is for
        controling DRC on pin2pin connections; It allows for explicitly
        connecting pins having a different width and or xs without raising a
        DRC error.

        Args:
            x (float): move in x-direction
            y (float): move in y-direction
            a (float): angle of rotation in degree
            xs (str): move to xs
            width (float): move to width
            radius (float): move to radius
            drc (bool): switch off drc if False (default=True)

        Returns:
            Node: a new Node translated by vector (<x>, <y>, <a>)
        """
        drc = drc
        node = self.copy()
        if cfg.clear_pin_on_move:
            if max(abs(x), abs(y)) > cfg. drc_max_distance:
                node.xs = None
                node.width = None
        if xs != 0:
            drc = False
            node.xs = xs
#           if xs is not None:
#               node.width = get_xsection(xs).width
        if width != 0:
            drc = False
            node.width = width
        if radius != 0:
            drc = False
            node.radius = radius
        connect_geo(self, node, (x, y, a), 0, 0, drc=drc)
        return node


    def rotate(self, a=0, xs=0, width=0, radius=0, drc=True):
        """Create a new Node rotate pointer by <a>.

        Also 'width' and 'xs' attributes of the can be 'moved'. This is for
        controling DRC on pin2pin connections; It allows for explicitly
        connecting pins having a different width and or xs without raising a
        DRC error.

        Args:
            a (float): angle of rotation in degree
            xs (str): move to xs
            width (float): move to width
            radius (float): move to radius
            drc (bool): switch off drc if False (default=True)

        Returns:
            Node: a new Node translated by vector (0, 0, <a>) w.r.t. self
        """
        drc = drc
        node = self.copy()
        if xs != 0:
            drc = False
            node.xs = xs
#            if xs is not None:
#                node.width = get_xsection(xs).width
        if width != 0:
            drc = False
            node.width = width
        if radius != 0:
            drc = False
            node.radius = radius
        connect_geo(self, node, (0, 0, a), 0, 0, drc=drc)
        return node
    rot = rotate


    def rot2ref(self, angle=0.0, ref=None, xs=0, width=0, radius=0, drc=True):
        """Create a new Node rotate pointer by <a> with respect to node ref.

        Also 'width' and 'xs' attributes of the can be 'moved'. This is for
        controling DRC on pin2pin connections; It allows for explicitly
        connecting pins having a different width and or xs without raising a
        DRC error.

        Args:
            angle (float): angle of rotation in degree w.r.t. ref (default is 0)
            ref (Node): if None, org of current cell is assumed
            xs (str): move to xs
            width (float): move to width
            radius (float): move to radius
            drc (bool): switch off drc if False (default=True)

        Returns:
            Node: a new Node translated by vector (0, 0, <a>) w.r.t. ref
        """
        drc = drc
        node = self.copy()
        if ref is None:
            ref = cfg.cells[-1].pin['org']
        x, y, a = diff(ref.rot(angle), self)
        if xs != 0:
            drc = False
            node.xs = xs
#            if xs is not None:
#                node.width = get_xsection(xs).width
        if width != 0:
            drc = False
            node.width = width
        if radius != 0:
            drc = False
            node.radius = radius
        connect_geo(self, node, (0, 0, -a), 0, 0, drc=drc)
        return node


    def skip(self, x=0, xs=0, width=0, radius=0, drc=True):
        """Create a new Node skipping <x> in direction of the pointer.

        Also 'width' and 'xs' attributes of the can be 'moved'. This is for
        controling DRC on pin2pin connections; It allows for explicitly
        connecting pins having a different width and or xs without raising a
        DRC error.

        Args:
            x (float): move in x-direction
            xs (str): move to xs
            width (float): move to width
            radius (float): move to radius
            drc (bool): switch off drc if False (default=True)

        Returns:
            Node: a new Node translated by vector (<x>, 0, 0) w.r.t. self
        """
        drc = drc
        node = self.copy()
        if cfg.clear_pin_on_move:
            if abs(x) > cfg. drc_max_distance:
                node.xs = None
                node.width = None
        if xs != 0:
            drc = False
            node.xs = xs
#            if xs is not None:
#                node.width = get_xsection(xs).width
        if width != 0:
            drc = False
            node.width = width
        if radius != 0:
            drc = False
            node.radius = radius
        connect_geo(self, node, (x, 0, 0), 0, 0, drc=drc)
        return node


    def shift(self, x=0, y=0, xs=0, width=0, radius=0, drc=True):
        """Create a new Node shifting pointer (<x>, <y>), keeping orientation.

        Also 'width' and 'xs' attributes of the can be 'moved'. This is for
        controling DRC on pin2pin connections; It allows for explicitly
        connecting pins having a different width and or xs without raising a
        DRC error.

        Args:
            x (float): move in x-direction
            y (float): move in y-direction
            xs (str): move to xs
            width (float): move to width
            radius (float): move to radius
            drc (bool): switch off drc if False (default=True)

        Returns:
            Node: a new Node translated by vector (<x>, <y,> 0) w.r.t. self
        """
        drc = drc
        node = self.copy()
        if cfg.clear_pin_on_move:
            if max(abs(x), abs(y)) > cfg. drc_max_distance:
                node.xs = None
                node.width = None
        if xs != 0:
            drc = False
            node.xs = xs
#            if xs is not None:
#                node.width = get_xsection(xs).width
        if width != 0:
            drc = False
            node.width = width
        if radius != 0:
            drc = False
            node.radius = radius
        connect_geo(self, node, (x, 0, 0), 0, 0, drc=drc)
        return node


    def offset(self, y=0, xs=0, width=0, radius=0, drc=True):
        """Create a new Node that is offset by <y>.

        Also 'width' and 'xs' attributes of the can be 'moved'. This is for
        controling DRC on pin2pin connections; It allows for explicitly
        connecting pins having a different width and or xs without raising a
        DRC error.

        Args:
            y (float): move in y-direction
            xs (str): move to xs
            width (float): move to width
            radius (float): move to radius
            drc (bool): switch off drc if False (default=True)

        Returns:
            Node: a new Node translated by vector (0, <y>, 0) w.r.t. self
        """
        drc = drc
        node = self.copy()
        if cfg.clear_pin_on_move:
            if abs(y) > cfg. drc_max_distance:
                node.xs = None
                node.width = None
        if xs != 0:
            drc = False
            node.xs = xs
#            if xs is not None:
#                node.width = get_xsection(xs).width
        if width != 0:
            drc = False
            node.width = width
        if radius != 0:
            drc = False
            node.radius = radius
        connect_geo(self, node, (0, y, 0), 0, 0, drc=drc)
        return node
    os = offset


    def flip(self):
        """Create a new Node that is flipped.

        A 'flip' mirrors the Node in the line of the pointer.
        Note that a flip state of a Node doesn't matter in Nazca.
        Flipping is handled by 'put', hence this funtion is equivalent to
        a Node copy in place (pointer not effected).

        Returns:
            Node: a new Node flipped w.r.t. self
        """
        node = self.copy()
        connect_geo(self, node, (0, 0, 0), 0, 0 )
        return node


    def flop(self):
        """Create a new Node that is flopped.

        A 'flop' mirrors the Node pointer 'backward'.
        A flop is effecively a flip + 180 degree rotation.
        Note that the flip state of a Node or not doesn't matter in Nazca,
        hence this flop is equivalent to rot(180).

        Returns:
            Node: a new Node flopped w.r.t. self
        """
        node = self.copy()
        connect_geo(self, node, (0, 0, 0), 0, 0 )
        return node


#==============================================================================
# netlist generators
#==============================================================================
    def path_nb_iter(self, sigtype):
        """Get path neighbours.

        Args:
            sigtype (str): type of neighbour connections to traverse

        Yields:
            Node: next neighbor
        """
        for nn in self.nb_opt:
            if nn[3].split('_')[0] == sigtype:
                yield nn[0], nn[1], nn[2], nn[3].split('_')[-1], nn[4]


    def geo_nb_iter(self):
        """Get geometrical neighbours.

        Yields:
            Node: next neighbour
        """
        for nn in self.nb_geo:
            #if nn[2]==1:
            yield nn


    def cnode_nb_iter(self, end=None, auxiliary=None):
        """Get cnodes at the cell level below the cell level of the node.

        Args:
            end (bool or None): include last element
            auxiliary (bool or None): filter on auxiliary flag. Ignore flag if None

        Yields:
            Node: next cell node
        """
        if end is None:
            cnodes = self.nb_cnode
        else:
            cnodes = self.nb_cnode[:end]
        for nn in cnodes:
            if nn[2] == 1:
                cnode = nn[0]
                if auxiliary is None:
                    yield cnode
                elif auxiliary is cnode.cell.auxiliary:
                    yield cnode
    instance_iter = cnode_nb_iter


    def print_neighbours(self):
        """Print neighbours of the node in the netlists.

        Returns:
            str: information on neighbors of this node
        """
        out = '\n'.join([
            '---neighbours geo:', '\n'.join(map(str, self.nb_geo)),
            '---neighbours opt:', '\n'.join(map(str, self.nb_opt)),
            '---neighbours ele:', '\n'.join(map(str, self.nb_ele))]
                       )
        return out


    def get_info(self):
        """Get information string.

        Returns:
            str: information on the node.
        """
        s = 'Node---------------------:\n  name: {}\n  cell: {}\n  node: {}'.format(
            self.name,
            self.cnode.cell.cell_name,
            self)
        s += '\n--cnode:\n    name: {}\n    cnode: {}'.format(
            self.cnode.name,
            self.cnode)
        #s += '\n--parent_cnode:\n    cell: {1}\n    cnode: {0}'.format(
        #     self.cnode.parent_cnode,
        #     self.cnode.parent_cnode.cell.cell_name)
        return s


#==============================================================================
# Functions to defining an edge (geometrical, optical, electrical) between nodes.
#==============================================================================

def pin_drc(pin1, pin2, t):
    """Perform DRC on pin-pin connections.

    Pin2pin DRC can be turned on and off.
    Included DRC are xsection, width and angle DRC.

    Args:
        pin1 (Node): first pin
        pin2 (Node): second pin
        t (tuple): (x, y, a) distance between pins, needed for angle.

    Returns:
        None
    """
    global drc_msgcnt
    if not (pin1.xs is not None and pin2.xs is not None):
        return None
    # else continue if both pins have a xs ->

    if not (pin1 is not pin1.cnode and pin2 is not pin2.cnode):
        return None
    # else continue if pins are not a cnode ->

    DRCxsection = DRCangle = DRCwidth = DRCsymmetry = False

    # check for xsection DRC:
    xs2 = cfg.drc_rule_xs.get(pin2.xs, [pin2.xs])
    if pin1.xs not in xs2:
        DRCxsection = True

    # check for symmetry DRC:
    elif not get_xsection(pin1.xs).symmetry:
        sameside = pin1.io == pin2.io
        sameflip = pin1.cnode.flip == pin2.cnode.flip
        sameparent = pin1.cnode.parent_cnode == pin2.cnode.parent_cnode
        if not (sameside != sameflip) and sameparent:
            DRCsymmetry = True

    # check for width DRC:
    drc_width1 = True
    drc_width2 = True
    if hasattr(get_xsection(pin1.xs), 'drc_width'):
        drc_width1 = get_xsection(pin1.xs).drc_width
    if hasattr(get_xsection(pin2.xs), 'drc_width'):
        drc_width2 = get_xsection(pin2.xs).drc_width
    if not (drc_width1 is False or drc_width2 is False):
        if pin1.width != pin2.width:
            if pin1.width is not None and pin2.width is not None:
                DRCwidth = True

    # check for angle DRC:
    drc_angle1 = True
    drc_angle2 = True
    if hasattr(get_xsection(pin1.xs), 'drc_angle'):
        drc_angle1 = get_xsection(pin1.xs).drc_angle
    if hasattr(get_xsection(pin2.xs), 'drc_angle'):
        drc_angle2 = get_xsection(pin2.xs).drc_angle
    if not (drc_angle1 is False or drc_angle2 is False):
        delta = cfg.drc_max_angle
        if abs(t[2]) > delta\
                and abs(abs(t[2])-180) > delta\
                and abs(t[2]-360) > delta\
                and abs(t[2]+360) > delta:
            DRCangle = True

    # Act on DRC violations:
    if DRCxsection or DRCangle or DRCwidth or DRCsymmetry:
        cp = cfg.cp
        # create pin reference strings for logging:
        if pin1.cnode.instance is False:
            pin1name = "{}.pin['{}']".format(pin1.cnode.cell.cell_name, pin1.up.name)
        elif pin1.cnode.instance is None: # None in case of the topcell.
            pin1name = "{}.pin['{}']".format(pin1.cnode.cell.cell_name, pin1.name)
        else:
            pin1name = "{}.pin['{}']".format(pin1.cnode.cell.cell_name, pin1.up.name)
        if pin2.cnode.instance is False:
            pin2name = "{}.pin['{}']".format(pin2.cnode.cell.cell_name, pin2.up.name)
        elif pin2.cnode.instance is None:
            pin2name = "{}.pin['{}']".format(pin2.cnode.cell.cell_name, pin2.name)
        else:
            pin2name = "{}.pin['{}']".format(pin2.cnode.cell.cell_name, pin2.up.name)

        if DRCxsection:
            cfg.drc_msgcnt += 1
            cfg.ND_msgcnt += 1
            msg = "DRC-{}/ND-{}: xsection mismatch in cell '{}': {} != {} : {} - {} @ {}".\
                format(cfg.drc_msgcnt, cfg.ND_msgcnt, cfg.cells[-1].cell_name, pin1.xs, pin2.xs,
                    pin1name, pin2name, pin1.fxy())
            Polygon(points=ring(*cfg.drc_ring_xs), layer=cfg.drc_layer_xs).put(pin1)
            font.text('x'+str(cfg.drc_msgcnt), layer=cfg.drc_layer_xs, height=5, align='cc').put(pin1)
            if (
                (cfg.drc_raise and (cfg.drc_msgcnt == cfg.drc_raisenum)) or
                (cfg.ND_raise and (cfg.ND_msgcnt == cfg.ND_raisenum))
            ):
                logger.error(msg, exc_info=True)
                raise Exception(msg)
            else:
                logger.error(msg)

        if DRCsymmetry:
            cfg.drc_msgcnt += 1
            cfg.ND_msgcnt += 1
            msg = "DRC-{}/ND-{}: symmetry mismatch in cell '{}' for xs = {} on pins {} to {} @ {}".\
                format(cfg.drc_msgcnt, cfg.ND_msgcnt, cfg.cells[-1].cell_name, pin1.xs,
                    pin1name, pin2name, pin1.fxy())
            Polygon(points=ring(*cfg.drc_ring_sym), layer=cfg.drc_layer_sym).put(pin1)
            font.text('s'+str(cfg.drc_msgcnt), layer=cfg.drc_layer_sym, height=5, align='cc').put(pin1)
            if (
                (cfg.drc_raise and (cfg.drc_msgcnt == cfg.drc_raisenum)) or
                (cfg.ND_raise and (cfg.ND_msgcnt == cfg.ND_raisenum))
            ):
                logger.error(msg, exc_info=True)
                raise Exception(msg)
            else:
                logger.error(msg)

        if DRCangle:
            cfg.drc_msgcnt += 1
            cfg.ND_msgcnt += 1
            msg = "DRC-{}/ND-{}: angle mismatch in cell '{}': {:.3f} != 0 deg between pins {} and {} @ {}".\
                format(cfg.drc_msgcnt, cfg.ND_msgcnt, cfg.cells[-1].cell_name, t[2], pin1name, pin1name, pin1.fxy())
            Polygon(points=ring(*cfg.drc_ring_angle), layer=cfg.drc_layer_angle).put(pin1)
            font.text('a'+str(cfg.drc_msgcnt), layer=cfg.drc_layer_angle, height=5, align='cc').put(pin1)
            if (
                (cfg.drc_raise and (cfg.drc_msgcnt == cfg.drc_raisenum)) or
                (cfg.ND_raise and (cfg.ND_msgcnt == cfg.ND_raisenum))
            ):
                logger.error(msg, exc_info=True)
                raise Exception(msg)
            else:
                logger.error(msg)

        if DRCwidth and not DRCxsection:
            cfg.drc_msgcnt += 1
            cfg.ND_msgcnt += 1
            msg = "DRC-{}/ND-{}: width mismatch in cell '{}': {} != {} : {} - {} @ {}".\
                format(cfg.drc_msgcnt, cfg.ND_msgcnt, cfg.cells[-1].cell_name, pin1.width, pin2.width,
                    pin1name, pin2name, pin1.fxy())
            Polygon(points=ring(*cfg.drc_ring_width), layer=cfg.drc_layer_width).put(pin1)
            font.text('w'+str(cfg.drc_msgcnt), layer=cfg.drc_layer_width, height=5, align='cc').put(pin1)
            if (
                (cfg.drc_raise and (cfg.drc_msgcnt == cfg.drc_raisenum)) or
                (cfg.ND_raise and (cfg.ND_msgcnt == cfg.ND_raisenum))
            ):
                logger.error(msg, exc_info=True)
                raise Exception(msg)
            else:
                logger.error(msg)

        cfg.cp = cp


def connect_geo(pin1, pin2, t=(0, 0, 0), chain_node=0, chain_connect=None,
        drc=True, solve=True):
    """Generate a geometrical edge between two Nodes.

    Connect pin1 and pin2 with geometrical transformation t.
    Pins are only allowed to connect if they the have the same parent cell,
    same cell or are part of an instance placed in the same parent cell.

    Args:
        pin1 (Node): 1st pin
        pin2 (Node): 2nd pin
        t (tuple): edge value: transformation (x, y, a) to get from pin1 to pin2
        chain_node (int): 1 if a chain node, 0 if not (default=0)
        chain_connect: overrules chain connect rules if None.
            chain_connect = 0 -> chain pin,
            chain_connect = 1 -> no-chain pin,
        (default=None)
        DRC (bool): apply DRC on the edge (default=True)

    Returns:
        None

    Exceptions:
        Pin connection out of scope.
    """
    #create a pointer in the nodes if none exists.
    assert pin1.pointer is not None
    #print("1: {}, 2: {}".format(pin1.pointer, pin2.pointer))
    if pin1.pointer is None:
        pin1.pointer = Pointer()
        pin1.pointer.chain = chain_node
    if pin2.pointer is None:
        pin2.pointer = Pointer()
        pin2.pointer.chain = chain_node

    if chain_connect is None:
        if pin1.cnode is pin2.cnode: # pins in same cell
            chain_connect = 0
        else:
            chain_connect = pin1.pointer.chain or pin2.pointer.chain

    mat = to_mat(chain_connect, *t)
    pin1.nb_geo.append((pin2, mat))
    pin2.nb_geo.append((pin1, inverse(mat)))

    if cfg.solve_direct and solve:
        if pin2.cnode is pin1.cnode:
            # pins placed in the same parent
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        elif pin2.cnode.parent_cnode is pin1.cnode:
            # pin2 in instance placed in the same parent as pin1
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        elif pin1.cnode.parent_cnode is pin2.cnode:
            # pin1 in instance placed in the same parent as pin2
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        elif pin1.cnode.parent_cnode is pin2.cnode.parent_cnode:
            # both are pins of instances placed in the same parent
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        else:
            raise Exception("Error: pin connection not in scope."\
                " pins reside in cells '{}' and '{}'.".\
                format(pin1.cnode.cell.cell_name, pin2.cnode.cell.cell_name))
    if cfg.drc and drc:
        if abs(t[0]) < cfg.drc_max_distance and abs(t[1]) < cfg.drc_max_distance: # points at same location
            pin_drc(pin1, pin2, t)


def connect_path(n1, n2, s, sigtype=None, path=None):
    """Generate a connection/edge between n1 and n1 with edge s.

    n1 (Node): first node
    n2 (Node): second node
    s (float): connection value
    sigtype (str): signal type, e.g. optical ('opt') or electrical ('dc')
    path(): list of pins for visualization

    Returns:
        None
    """
    # TODO: flag that active cell has an optical/sigtype connection
    if sigtype is None:
        sigtype = 'dis'
    if n1 is not None:
        n1.nb_opt.append((n2, s, 1, sigtype, path))
    if n2 is not None:
        n2.nb_opt.append((n1, s, -1, sigtype, path))
    # None neighbours represent terminations.



def connect_optical_path(n1, n2, s, sigtype=None, path=None):
    """Generate a connection/edge between n1 and n1 with edge s.

    n1 (Node): first node
    n2 (Node): second node
    s (float): connection value
    sigtype (str): signal type, e.g. optical ('opt') or electrical ('dc')
    path(): list of pins for visualization

    Returns:
        None
    """
    # TODO: flag that active cell has an optical/sigtype connection
    if sigtype is None:
        sigtype = 'dis'
    sig1=sigtype.split('_')[0]
    sig2=sigtype.split('_')[-1]
    if n1 is not None:
        n1.nb_opt.append((n2, s, 1, sig1+'_'+sig2, path))
    if n2 is not None:
        n2.nb_opt.append((n1, s, -1, sig2+'_'+sig1, path))
    # None neighbours represent terminations.




def connect_cnode(n1, n2, s):
    """Generate cell connection.

    Returns:
        None
    """
    n1.nb_cnode.append((n2, s, 1, 0))
    n2.nb_cnode.append((n1, s, -1, 0))


def _connect_ribbon(pin1, pin2, N0):
    """Helper function to connect all pins from ribbons rib0 and rib1.

    This function is called by the ribbon_connect function specifically to
    connect the ribbon pins.

    Args:
        rib0 (Ribbon pin): first Ribbon that is put
        rib1 (Ribbon pin): second Ribbon that is put
        N0 (int): amount of waveguides of both ribbons

    Returns:
        None
    """
    order1 = list(range(N0))
    order0 = range(N0)

    flip = pin1.cnode.flip != pin2.cnode.flip
    sameside = pin1.up.type == pin2.up.type

    if (sameside and not flip) or (not sameside and flip):
        order1.reverse()

    # TODO: store prefix in instance for flexibility instead of assuming 'a', 'b'
    if pin1.up.type == 3:
        prefix1 = "a"
    else: # TODO: check explicitly for 4
        prefix1 = "b"

    if pin2.up.type == 3:
        prefix2 = "a"
    else: # TODO: check explicitly for 4
        prefix2 = "b"

    #print(f"ribbons: {pin1.cnode.cell.cell_name}, {pin2.cnode.cell.cell_name}")
    for ord0, ord1 in zip(order0, order1):
        #print(f"{prefix1}{ord0}", f"{prefix2}{ord1}") # print connection made
        connect_geo(pin1.cnode.instance.pin[f"{prefix1}{ord0}"], pin2.cnode.instance.pin[f"{prefix2}{ord1}"])


def parse_pin(C1=None, C2=None, C3=None, C4=None, instance=None, rot=False,
        default='out'):
    """Parse the part of the put statement after a potential first string arg.

    Args:
        C1: tuple of floats, len<=3; C2 = pinname
        C1 to C[i]: 1<=i<=3 are floats; C[i+1] = pinname
        instance (Instance object): Instance object to connect to
        rot (bool): Add an extra rotation for end pins in p2p connections if True

    Returns:
        pin, (x, y, a): describes a point at 'translation' from 'pin'.
            default='org', (0, 0, 0).

    Examples:
        Allowed input formats:

        * ()
        * (Node)
        * ((0))
        * ((0, 0))
        * ((0, 0, 0))
        * (0,)
        * (0, 0)
        * (0, 0, 0)

        * 'pin' can be referred to by a Node, Instance or string name.
        * (pin)
        * (None, pin)
        * ((0), pin)
        * ((0, 0), pin)
        * ((0, 0, 0), pin)
        * (0, pin)
        * (0, 0, pin)
        * (0, 0, 0, pin)
    """
    if C1 is None:
        return cfg.cp, (0, 0, 0)

    xya = [0, 0, 0]
    if isinstance(C1, tuple):
        P = C2
        for i, elm in enumerate(C1):
            if isinstance(elm, (int, float)):
                xya[i] = float(elm)
    else:
        P = C4
        for i, elm in enumerate([C1, C2, C3]):
            if isinstance(elm, (int, float)):
                xya[i] = float(elm)
            else:
                P = elm
                break

    if P is None: #position w.r.t. 'org'
        if instance is None:
            if rot:
                xya[2] += 180
            return cfg.cells[-1].cnode, xya
        else:
            return instance.cnode.parent_cnode, xya
    elif isinstance(P, str):
        return cfg.cells[-1].pin[P], xya
    elif isinstance(P, Node):
        return P, xya
    elif isinstance(P, Instance):
        if default == 'out':
            return P.pin[P.cnode.cell.default_out], xya
        else:
            return P.pin[P.cnode.cell.default_in], xya
    else:
        raise Exception('pin could not be parsed {}, {}, {}, {}'.\
            format(C1, C2, C3, C4))


def validate_basename(basename):
    """Check if basename is unique on substring level.

    This is only relevant for pcell replacement where no mistake between
    different basenames is allowed.

    Args:
        basename (str): basename of a cell's name

    Returns:
        bool: True if basename is valid
    """
    size1 = len(basename)
    for bname in cfg.basenames:
        size2 = len(bname)
        size = min(size1, size2)
        #print(bname[0:size], basename[0:size])
        if bname[0:size] == basename[0:size]:
            main_logger(f"basename overlap {bname} - {basename}", "error")
            return False
    return True


class Instance():
    """Class to store pin information and properties of cell instantiations.
    """
    id = 0

    def __init__(self, name=None):
        """Construct an Instance object.

        Returns:
            None
        """
        self.name = name
        self.cnode = None
        self.cell = None
        self.parent_cnode = None
        self.pin = dict()
        self.default_in = None
        self.default_out = None
        self.org = None
        self.flip = None
        self.array = None
        self.id = Instance.id
        self.bbox = None # bbox position of the instance in the parent cell
        Instance.id += 1 #keep instance counter unique

    @property
    def length_geo(self):
        try:
            lgeo = self.cnode.cell.length_geo
        except AttributeError as E:
            #cfg.print_except("bare-except".)
            logger.exception(E)
            lgeo = 0

        return lgeo

    def __repr__(self):
        return "<Instance() object of cell '{}' in cell '{}'>".format(
            self.cnode.cell.cell_name,  self.cnode.parent_cnode.cell.cell_name)

    def ic_pins(self):
        """Generator over interconnect pins, filter out org and bounding box.

        Yields:
            str, Pin: iterator over pins in cell
        """
        for name, p in self.pin.items():
            if (name != 'org') and (name not in bbu.bbox_pinnames):
                yield name, p


    def raise_pins(self, namesin=None, namesout=None, show=True):
        """Copy pins of a cell instance to its parent cell.

        If a cell A is "put" in cell P, then an instance of cell A is created
        inside "parent" cell P. Each instance of cell A, say A1, A2, ...,
        contains a copy of all pins of cell A to represent cell A in P at a
        particular location, rotation, flip state and/or scaling. raise_pins
        automatically copies all pin attributes.

        The instance pins are themselves not pins of cell P.
        Pins have to be explicitly set in P with nd.Pin.put().
        Alternatively, all or part of the pins in an instance
        can be set at once in P by "raise_pins", avoiding many nd.Pin().put()
        statements. This can very useful when an instance has a large number
        of pins that have to be "raised".

        Note:
            raise_pins is a method of Class Instance, not of Cell. Hence,
            it operates on a cell that has been put.

        Args:
            namesin (list of str): list of pin names to raise
                (default=all pins of the instance)
            namesout (list of str): list of new names of the raised pins
                (default=<namesin>)
            show (bool): show the pins in the layout with an arrow (default=True)

        Returns:
            None

        Example:
            Raise (copy) all pins of an instance to its parent cell.
            In the example a 90 degree bend is connected to port 'a0'
            of the new_cell.

            Create pins without the `raise_pins` method looks like this::

                # create and put Pins explicitly one by one:
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example_cell().put(100, 200)
                    nd.Pin('a0').put(instance1.pin['a0'])
                    nd.Pin('b0').put(instance1.pin['b0'])

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['a0'])
                nd.export_plt()

            Now **raise** pins with `raise_pins` on <instance1> in the new_cell.
            This also copies all pin attributes::

                # raise pins
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example_cell().put(100, 200)
                    instance1.raise_pins()

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['a0'])
                nd.export_plt()

            Now **raise and rename** pins 'a0' and 'b0' of <instance1>
            with method `raise_pins`. Note that the
            cell default pins, 'a0' and 'b0', are automatically added to
            <new_cell>::

                # raise and rename pins
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example_cell().put(0)
                    instance1.raise_pins(['a0', 'b0'], ['new_a0', 'new_b0'])

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['new_a0'])
                nd.bend(angle=-90).put(instance2.pin['a0'])
                nd.export_plt()
        """
        do_not_raise = ['org']
        if namesin is None:
            for name, node in self.pin.items():
                if name not in do_not_raise:
                    Pin(name, width=node.width, xs=node.xs, type=node.type,
                        chain=node.chain, show=show, remark=node.remark).put(node)
        else:
            if namesout is None:
                namesout = namesin
            for namein, nameout in zip(namesin, namesout):
                if namein in self.pin.keys():
                    node = self.pin[namein]
                    Pin(nameout, width=node.width, xs=node.xs, type=node.type,
                        chain=node.chain, show=show, remark=node.remark).put(node)
                else:
                     main_logger(f"'raise_pins': pin '{namein}' not defined.", "warning")

num = 0
class Cell():
    """Cell object, corresponding to a GDS cell when exported to GDS.

    A cell has default pins:
        'org': cell origin (non-chain)
        'a0' : input pin   (chain)
        'b0' : output pin  (chain)

    It is up to the user to place 'a0' and 'b0' at sensible positions in the cell.
    If not provided, they will be placed at the origin automatically.

    When starting a new cell: cp = C.pin['org']

    A cell typically contains one or more of the following objects:
        1. pins: A position w.r.t. the cell origin that can be used to connect to.
            method: Pin().put()

        2. cell instances: A reference to another cell
            method: Cell().put()

        3. polygons: A polygon assigned to a specific GDS layer.
            method: Polygon().put()

        4. polylines: A polyline/path assigned to a specific GDS layer.
            method: Polyline().put()

        5. annotations:
            method: Annotation().put()
    """
    id = 0 # cell counter/ID

    def __init__(
        self,
        name='cell',
        celltype='element',
        instantiate=True,
        autobbox=False,
        hull=None,
        store_pins=None,
        hashme=False,
        cnt=False,
        params=False,
        show_hull=False,
        hull_based_bbox=None,
    ):
        #TODO: params unused?
        """Construct a Cell object.

        Cell definitions can be nested.

        If the cell is defined as the first cell inside a function that has
        a @hashme decorator, the cell obtains its (unique) name from hashme.

        Args:
            name (str): cell name (default='cell')
            celltype (str): type of cell, option 'element', 'mask' (default='element')
            instantiate (bool): flag if the cell is instantiated (default=True)
            autobbox (bool): create a cell bounding box (default=False)
            hull (bool): calculate bbox as a convex hull (may be slower)
            store_pins (bool): store cell pins in annotation inside GDS
            hashme (bool): optionally set True if cell obtains information from
                hashme decorator, but this is autodetected as well.
            cnt (bool): append ordinal counter to the cellname (default=False)
            params (bool): add parameter annotation in the cell (default=False)
                Note that autobbox = True will already add parameter annotation.

        Returns:
            None

        Example:
            Create a cell with cell name 'mycell' and variable name C.
            Put the cell in the mask layout and export the layout::

                import nazca as nd

                with nd.Cell('mycell') as C:
                    nd.strt(length=40).put(0)
                    nd.bend(angle=15).put()

                C.put(0)

                nd.export_plt()
        """
        global num
        num += 1
        self.properties = {}
        self.hashme = hashme
        try: # detect if hashme is on, if so, use it for the cell name via True
            hme = cfg.hashme
            if hme:
                self.hashme = True
                # reset because there may be a cell in a cell in hashme
                cfg.hashme = False
        except:
            pass
        if self.hashme:
            name = cfg.hash_name
            if name == '':
                raise Exception("Detected a call nd.Cell(hashme=True) without first calling the @hashme('cellname') decorator with a non-empty 'cellname'.")
            basename = cfg.hash_basename
            paramsname = cfg.hash_paramsname
            self.func_id = cfg.hash_func_id
            self.parameters = cfg.hash_params
            self.properties['parameters'] = cfg.hash_params
            self.properties['cellnameparameters'] = [key for key in cfg.hash_cellnameparams]
        else:
            basename = cfg.gds_cellname_cleanup(name)
            if cfg.validate_basename:
                validate_basename(basename)
            paramsname = cfg.gds_cellname_cleanup(name)
            if cnt:
                name = '{}_{}'.format(name, num)
                name = cfg.gds_cellname_cleanup(name)
            self.parameters = None
            self.properties['parameters'] = None
            self.func_id = None

        name = cfg.gds_cellname_cleanup(name)
        if name in cfg.cellnames.keys():
            #TODO: what if new name also exists?
            name2 = '{}_{}'.format(name, num)
            cfg.cellnames[name2] = self
            if instantiate is True:
                gdsmsg = ""
                if cfg.gdsload != "" and cfg.gdsload is not False:
                    gdsmsg = f" Occured during gds load of '{cfg.gdsload}'. Use load_gds(cellsreused=['{name}']) to reuse existing cells instead."
                main_logger(
                    f"Duplicate cell name in nd.Cell(name='{name}') renamed to '{name2}'." + gdsmsg,
                    "warning"
                )
            name = name2
        else:
            cfg.cellnames[name] = self
        self.properties['cell_name'] = name

        self.id = Cell.id
        Cell.id += 1
        cfg.cells.append(self)
        cfg.self = cfg.cells[-1]
        self.bbinfo = dict() # store metadata on the cell
        self.closed = False
        #self.id = next(Et)
        self.pin2 = None # store pinout for p2p interconnect cells
        self.cnode = Node('org')
        self.cnode.cnode = self.cnode
        self.cnode.cell = self
        self.cnode.pointer = Pointer(0, 0, 0)
        self.cnode.flip = False
        self.cnode.instance = None
        self.org = self.cnode
        self.parent_cnode = None

        self.polygons = []
        self.polylines = []
        self.gdsfiles = []
        self.annotations = []
        self.name = self.id
        self.cell_basename = basename     # name: base
        self.cell_paramsname = paramsname # name: base + params
        self.cell_name = name             # name: base + params + hash
        self.instantiate = instantiate
        self.group = ''  # (str) optional value for grouping/sorting cells
        self.name_layer = None  # (str) custom layer for displaying BB name polygon. None defaults to 'bbox_name'.
        self.module = '' # optional str value to store the module name the cells is defined in.
        self.ranges = {} #  for storing parameter information like the allowed range
        self.bbox = None
        self.autobbox = autobbox
        self.bboxbuf= 0
        self.hull = None  # contain the actual hull data if any.
        self.userbbox = None
        self.add_bbox_flag = False  # store if bbox has been calculated already
        self.userbbox_include = True  # make userbbox part of the bbox if True
        self.protectbox = False  # copy the bbox to the userbbox if no userbbox is defined.
        self.hull_based_bbox = None  # hull based bbox if True
        self.ibbox = {}  # {inode, ibbox} store instance bounding boxes for group_connect
        self.params = params
        self.version = None
        self.auxiliary = False  # indicate if the cell is not real layout, e.g. pins, and stubs
        self.bbox_polygon = []

        if store_pins is None:
            self.store_pins = cfg.store_pins
        else:
            self.store_pins = store_pins

        if show_hull is None:
            self.show_hull = cfg.show_hull
        else:
            self.show_hull = show_hull  # export the hull as layout


        # Note that if use_hull is True:
        # - calculate a hull for this cell
        # - AND use hulls of instances in this cell if available (otherwise the bbox)
        if hull is None:
            self.use_hull = cfg.use_hull
        elif hull is True:
            self.use_hull = True
        else:
            self.use_hull = False

        if hull_based_bbox is None:
            self.hull_based_bbox = cfg.hull_based_bbox
        else:
            self.hull_based_bbox = hull_based_bbox

        if hull is True and hull_based_bbox:
            self.hull_based_bbox = True

        #pin handling:
        self.oldcp = cfg.cp
        self.pin = {self.cnode.name: self.cnode}
        if celltype == 'mask':
            self.default_in = self.cnode.name # 'org'
            self.default_out = self.cnode.name # 'org'
        else: #'element':
            self.default_in = 'a0'
            self.default_out = 'b0'

        cfg.cp = self.org


    def filter(layers):
        """Create a new cell after applying a layer filter.

        Args:
            layers (list of layer): layers to keep

        Returns:
            Cell: new cell with subset of original layers

        To be implemented
        """


    def __repr__(self):
        return "<Cell(name='{}', instantiate={})>".\
            format(self.cell_name, self.instantiate)


    def default_pins(self, pinin=None, pinout=None):
        """Set default input and output pin.

        Args:
            pinin (Node | str): pin to set as default cell input pin
            pinout (Node | str): pin to set as default cell output pin

        Returns:
            None
        """
        if pinin is None:
            pinin = self.default_in
        else:
            self.default_in = pinin

        if pinout is None:
            pinout = self.default_out
        else:
            self.default_out = pinout


    def add_property(self, item):
        """Add item to the cell properties dictionary.

        item (dict): new property

        Returns:
            None
        """

        def merge_properties(D1, D2):
            """Merge dictionaries recursively.

            D1 (dict): existing properties
            D2 (dict): properties to add

            Returns:
                dict: merged level of the dictionary
            """
            update = {}
            if isinstance(D1, dict):
                D1keys = list(D1.keys())
            else:
                D1keys = []
            if isinstance(D2, dict):
                D2keys = list(D2.keys())
            else:
                D2keys = []

            if not D1keys or not D2keys:
                main_logger(
                    f"In cell '{self.cell_name}': Overwriting properties: {D1} to {D2}",
                    "warning"
                )
                return D2
            else:
                for e in set(D1keys + D2keys):
                    if e not in D1keys:
                        update[e] = D2[e]
                    elif e not in D2keys:
                        update[e] = D1[e]
                    else:
                        update[e] = merge_properties(D1[e], D2[e])
                #print("update", update)
                return update

        self.properties = merge_properties(self.properties, item)
        #print("\n", self.properties)
        return None


    def ic_pins(self):
        """Generator over interconnect pins, filter out org and bounding box.

        Yields:
            str, Pin: iterator over pins in cell
        """
        for name, p in self.pin.items():
            if (name != 'org') and (name not in bbu.bbox_pinnames):
                yield name, p


    def __enter__(self):
        return self


    def __exit__(self, type, value, traceback):
        if value is None:
            self.close()
        return False


    def parse_instance_connection(self, C1, C2, C3, C4, C5, instance=None):
        """Connect a (closed) cell as an Instance object into the active cell.

        Parse up to five connection arguments to extract: pin1, pin2 and (x,y,a):

        Args:
            C1 (float | Node | Instance): connection information
            C2 (float | Node | Instance): connection information
            C3 (float | Node | Instance): connection information
            C4 (float | Node | Instance): connection information
            C5 (float | Node | Instance): connection information
            instance (Instance): instance object to connect to

        Returns:
            Node, Node, (x, y, a): instance pin, cell pin, connection

        Notes:
            The connection has the following syntax:

            [pin_inst [, translation]]

            or

            [pin_inst, ] translation  [, pin_cell]

            * translation = x | x,y | x, y, a | (x) | (x, y) | (x, y, z), where
              translation is w.r.t. <pin_cell>
            * pin_inst: pinname of instance (default='a0')
            * pin_cell: pinname of cell (default='org')

            Options:

            * ()
            * (node)
            * (instance)
            * (0,)
            * (0, 0)
            * (0, 0, 0)
            * ((0))
            * ((0, 0))
            * ((0, 0, 0))
            * (0, pin)
            * (0, 0, pin_cell)
            * (0, 0, 0, pin_cell)
            * ((0), pin_cell)
            * ((0, 0), pin_cell)
            * ((0, 0, 0), pin_cell)

            * (pin_inst)
            * (pin_inst, 0)
            * (pin_inst, 0, 0)
            * (pin_inst, 0, 0, 0)
            * (pin_inst, (0))
            * (pin_inst, (0, 0))
            * (pin_inst, (0, 0, 0))

            * (pin_inst, pin_cell)
            * (pin_inst, 0, pin)
            * (pin_inst, 0, 0, pin_cell)
            * (pin_inst, 0, 0, 0, pin_cell)
            * (pin_inst, (0), pin_cell)
            * (pin_inst, (0, 0), pin_cell)
            * (pin_inst, (0, 0, 0), pin_cell)
        """
        if C1 is None:
            return instance.pin[self.default_in], cfg.cp, (0, 0, 0)
        elif isinstance(C1, Node):
            if C1.cnode.instance is None and C1.cnode.cell.closed:
                raise Exception("You are trying to connect to a closed Cell object:\n"\
                    "  The construction found is similar to\n"\
                    "  $ foo.put(cell.pin['a0'])\n"\
                    "  Connect not to the Cell but to an instance of the cell instead\n"\
                    "  $ instance = cell.put()\n"\
                    "  $ foo.put(instance.pin['a0'])\n")
            return instance.pin[self.default_in], C1, (0, 0, 0)
        elif isinstance(C1, Instance):
            if C1.cnode.instance is None and C1.cnode.cell.closed:
                raise Exception("Trying to connect to a closed cell. "\
                    "Try using an instance instead.")
            return instance.pin[self.default_in], C1.pin[C1.cnode.cell.default_out], (0, 0, 0)

        elif isinstance(C1, str):
            pin2, T = parse_pin(C2, C3, C4, C5, instance)
            return instance.pin.get(C1), pin2, T
        else:
            pin2, T = parse_pin(C1, C2, C3, C4, instance)
            return instance.pin[self.default_in], pin2, T


    def parse_pin_connection(self, C1, C2, C3, C4):
        """Parse pin connection for polygon, polyline, annotation and gds methods.

        Args:
             C1, C2, C3, C4 (float | Node | instance): Parsing information

        Returns:
            Node, (x, y, a): pin, translation

        Example:
            Options:

            * ()
            * (node)
            * (0,)
            * (0, 0)
            * (0, 0 ,0)
            * ((0))
            * ((0, 0))
            * ((0, 0 ,0))
            * (0, 'b0')
            * (0, 0, 'b0')
            * (0, 0 ,0, 'b0')
            * ((0), 'b0')
            * ((0, 0), 'b0')
            * ((0, 0 ,0), 'b0')
        """
        #TODO: check if pin connection is in scope (parent or sibling)
        if C1 is None or (isinstance(C1, tuple) and len(C1) == 0):
            return cfg.cp, (0, 0, 0)
        elif isinstance(C1, Node):
            if C1.cnode.instance is None and C1.cnode.cell.closed and not cfg.patchcell:
                raise Exception("It is not allowed to connect to a closed cell:"
                    " from cell '{2}'"
                    " to pin '{1}' of the already *closed* cell '{0}'."
                    " Solution: instantiate cell '{0}' first inside cell '{2}'"
                    " using put(),"
                    " and connect to a pin of the resulting instance instead.".
                    format(C1.cnode.cell.cell_name, C1.name, self.cell_name))
            return C1, (0, 0, 0)
        elif isinstance(C1, Instance):
            if C1.cnode.instance is None and C1.cnode.cell.closed and not cfg.patchcell:
                raise Exception("It is not allowed to connect from cell '{2}'"
                    " to pin '{1}' of the already *closed* cell '{0}'."
                    " Solution: instantiate cell '{0}' first inside cell '{2}'"
                    " using put(),"
                    " and connect to a pin of the resulting instance instead.".
                    format(C1.cnode.cell.cell_name, C1.name, self.cell_name))
            return C1.pin[C1.cnode.cell.default_out], (0, 0, 0)

        xya = [0, 0, 0]
        if isinstance(C1, tuple):
            for i, elm in enumerate(C1):
                # isinstance(elm, int) may fail on numpy.int64 type
                try:
                    elm = 1.0*elm
                except TypeError as E:
                    raise Exception("Invalid value provided pin({}). {}".\
                        format(repr(elm), E))
                xya[i] = elm
            if isinstance(C2, str):
                return self.pin[C2], xya
            else:
                return self.org, xya
        else:
            P = C4
            for i, elm in enumerate([C1, C2, C3]):
                # isinstance(elm, int) may fail on numpy.int64 type
                if elm is None:
                    P = elm
                    break
                else:
                    try:
                        elm = 1.0*elm
                    except TypeError as E:
                        raise Exception("Invalid value provided: put({})."\
                            "This may be caused by using a pin name string in a Polygon put().\n {}".\
                            format(repr(elm), E))
                    xya[i] = elm
            if isinstance(P, str):
                return self.pin[P], xya
            else:
                return self.org, xya


    def _put_pin(self, name=None, connect=None, xs=None, width=0, type=None,
            io=None, chain=1, show=False, remark=None):
        """Add new pin and connect it to an existing pin.

        Returns:
            Node: pin being put
        """
        node_new = Node(name) #next(Pt)
        node_new.pointer = Pointer()
        node_new.pointer.chain = chain
        node_new.cnode = self.cnode
        node_new.xs = xs
        node_new.width = width
        node_new.io = io
        node_new.type = type
        node_new.show = show
        node_new.remark = remark

        if name is not None: # add a node
            self.pin[name] = node_new

        C = [None] * 4
        if isinstance(connect, tuple):
            for i, conn in enumerate(connect):
                C[i] = conn
        else:
            C[0] = connect

        node_exist, T = self.parse_pin_connection(*C)
        connect_geo(node_exist, node_new, t=T, chain_node=1, chain_connect=0)
        return node_new


    def get_pin(self, pinname=None, connect=None):
        """"Parse pinname and connect for put_polygon and put_gds methods.

        Returns:
            Node: pin position, new if needed
        """
        if pinname is None:
            pinname = 'org'
        if pinname not in cfg.cells[-1].pin.keys(): #add new pin w.r.t. org
            if connect is None: # org
                node = self._put_pin(pinname, cfg.cells[-1].org)
            elif isinstance(connect, tuple): #(x,y,a)
                node = self._put_pin(pinname, connect)
                #node = cfg.cells[-1].org.move(*connect)
                #node.name = pinname
            elif isinstance(connect, Node):  #Node
                self._put_pin(pinname, connect)
                node = connect
            return node
        else: #existing pin
            if connect is None: # org
                node = cfg.cells[-1].org
            elif isinstance(connect, tuple): #(x,y,a)
                node = cfg.cells[-1].pin[pinname].move(*connect)
            elif isinstance(connect, Node):  #Node
                node = connect #maybe set parent node !!!!
            return node


    def _put_polygon(self, connect=None, polygon=None, flip=False, flop=False,
            scale=1.0):
        """Add a polygon to the cell.

        Args:
            connect (float, float, float): position (x, y, a) of polygon in cell.
            polygon (Polygon)

        Returns:
            None
        """
        node = self._put_pin(None, connect)
        node.flip = flip
        node.flop = flop
        node.scale = scale
        if polygon is not None:
            self.polygons.append((node, polygon))
        return None


    def _put_polyline(self, connect=None, flip=False, flop=False, polyline=None,
        scale=1.0):
        """Add a polyline to the cell.

        Args:
            connect (float, float, float): position (x, y, a) of polygon in cell.
            polyline (Polyline)

        Returns:
            None"""
        node = self._put_pin(None, connect)
        node.flip = flip
        node.flop = flop
        node.scale = scale
        if polyline is not None:
            self.polylines.append((node, polyline))
        return None


    def _put_annotation(self, connect=None, annotation=None):
        """Add an annotation to the cell.

        Args:
            pinname (str): pinname to attach the annotation to.

        Returns:
            None
        """
        node = self._put_pin(None, connect)
        if annotation is not None:
            self.annotations.append((node, annotation))
        return None


    def _put_gds(self, connect=None, filename=None, cellname=None,
            newcellname=None, layermap=None, cellmap=None, scale=1.0, strm=None):
        """Add a GDS to the cell. Anchor it to a pin by the pinname (string).

        Example::

            cell_1.put_gds(filename='path_to_file.gds',
                cellname='name', newcellname='new')
        """
        node = self._put_pin(connect=connect)
        self.gdsfiles.append((node,
            (filename, cellname, newcellname, layermap, cellmap, scale, strm)))
        #print('put_gds node {}, pinname {} = '.format(node, pinname))
        return None


    def _copy_cell2instance(self, parent_cnode, flip, flop, array=None, scale=1.0):
        """Copy all nodes from the Cell into an Instance.

        Args:
            parent_node (Node): main cell node of the cell to copy.
            flip (bool): flip state of the instance.
            array (list): creates an array instance upon GDS output if not None.
                Format the array as [col#, [dx1, dy1], row#, [dx2, dy2]]
                dx1, dy1 is the column vector measured between adjacent elements.
                dx2, dy2 is the row vector measured between adjacent elements.
                (default=None)

        Returns:
            Instance: instance created from <parent_node>
        """
        instance = Instance(name=f"I_{self.cell_name}")

        flipflop = flip != flop
        if array is not None:
            try:
                A = array
                array = [ A[0], [ A[1][0], A[1][1] ], \
                          A[2], [ A[3][0], A[3][1] ] ]
                if not self.instantiate:
                    main_logger(
f"Trying to set array on cell '{self.cell_name}' with instantiate=False"
f" in parent cell '{self.cell.parent_cnode.cell.cell_name}'."
f" Only a single instance will be visible."
f" Set instantiate=True to obtain the full array in GDS.",
"warning"
                    )
                    array = None
            except Exception as excp:
                logger.exception(
f"Error: incompatible array format {array} when placing "
f"cell '{self.cell_name}' in '{parent_cnode.cell.cell_name}'. "
f"Array will be ignored.\n{excp}",
"exception"
                )
                array = None

        # cnode:
        node = Node()
        #node.name = '{}_cnode_{}'.format(self.cell_name, Cell.I)
        node.name = '{}_i{}'.format(instance.id, 'org')
        node.up = self.cnode
        node.instance = instance
        node.pointer = self.cnode.pointer.copy()
        node.cnode = node
        node.xs = self.cnode.xs
        node.width = self.cnode.width
        node.type = self.cnode.type
        node.cell = self
        node.flip = flipflop
        node.array = array
        node.scale = scale
        instance.cnode = node
        instance.cell = node.cell
        instance.cnode.parent_cnode = parent_cnode
        instance.pin['org'] = node
        instance.org = node
        instance.properties = self.properties

        connect_cnode(parent_cnode, instance.cnode, None)
        # all other nodes:
        for key, pin in self.pin.items():
            if True: #the name 'org' maybe outside the cnode. #key is not 'org':
                # TODO: do all properties have to be copied?
                #       maybe better to get them from the cell?
                node = Node()
                node.up = pin  # link to original node in cell
                node.xs = pin.xs
                node.io = pin.io
                node.width = pin.width
                node.type = pin.type
                node.show = pin.show
                node.remark = pin.remark
                # node.name = '{}_{}_{}'.format(self.cell_name, key, Cell.I)
                node.name = '{}_{}'.format(instance.id, key)
                node.cnode = instance.cnode
                instance.pin[key] = node
                # reconstruct connectivity
                # 1. org-to-node connectivity
                #    no solving needed, this will be done in put statement.
                x, y, a = pin.pointer.xya()
                if flipflop:
                    y, a = -y, -a
                connect_geo(instance.cnode, node, (x, y, a), solve=False)
                # 2. copy pointer properties (stored position xya is irrelevant)
                node.pointer = pin.pointer.copy()

        instance.pinin = instance.pin[self.default_in]
        instance.pinout = instance.pin[self.default_out]
        return instance


    # TODO: add instantiate to put to flatten at put-time
    def put(self, *args, flip=False, flop=False, array=None, scale=1.0,
            drc=None, **kwargs):
        """Instantiate a Cell object and put it in the active cell.

        Args:
            *args (float | Node | Instance): a set of maximum 5 unnamed parameters
                interpreted as a connection between two pins, where
                put() connects to the current pin cp
                For more details see method 'parse_instance_connection'.
            flip (bool): mirror the instance in the vector line (default=False)
                flip is preformed after any translation and/or rotation.
                Hint: for connecting waveguides you only need flip (not flop).
            flop (bool): mirror the instance in the line perpendicular to the vector
                (default=False)
                flop is performed after any translation and/or rotation.
                A flop is the same as rot(180) + a flip
                Hint: for mask assembly or connecting bbox pins you may find both
                flip and flop useful.
            cp (str): Set current pin cp to a specific pin after putting
                (default='b0'). Example: cp='b1'
            array (list): creates an array instance upon GDS output if not None.
                Format the array as [col#, [dx1, dy1], row#, [dx2, dy2]]
                dx1, dy1 is the column vector measured between adjacent elements.
                dx2, dy2 is the row vector measured between adjacent elements.
                (default=None)
            scale (float): set scaling factor of the instance (default=1.0).
                In gds the scaling will be set as instance attribute. In case
                of flattening the instance, the scaling factor will be applied
                to the instance content.

        Returns:
            Instance: reference to the instance of the cell being put

        Examples:
            Below are different syntax examples for put.

            In the examples we use the following objects:
                * nd.example: a pre-defined cell in Nazca
                * 'a0': default input pin name of "example"
                * 'b0': default output pin name of "example"


            1. connect the default pin of cell "example" to current pin (cp)::

                import nazca as nd
                nd.example_cell().put() # -> put('a0', cp)
                nd.export_plt()

            2. connect pin 'b0' of cell "example_cell()" to current pin (cp)::

                import nazca as nd
                nd.example_cell().put('b0') # -> put('b0', cp)
                nd.export_plt()

            3. connect default pin of cell "example_cell()" to its parent cell at org + (10, 20, 30)::

                import nazca as nd
                nd.example_cell().put(10, 20, 30) # -> put('a0').move(10, 20, 30)
                nd.export_plt()

            4. connect default pin of cell "example_cell()" to instance C at pin 'b0'::

                import nazca as nd
                C = nd.example_cell().put(0)
                nd.example.put(C.pin['b0']) # -> put('a0', C.pin['b0'])
                nd.export_plt()

            5. connect pin 'b0' of cell "example_cell()" to instance C at pin 'b0'::

                import nazca as nd
                C = nd.example_cell().put(0)
                nd.example_cell().put('b0', C.pin['b0']) # -> put('b0', C.pin['b0'])
                nd.export_plt()

            6. connect pin 'b0' of cell "example_cell()" to instance C at its default out pin, i.e.
               refering to instance C only without pin attribute is interpreted as
               connecting to the default output pin of C.::

                import nazca as nd
                C = nd.example_cell().put(0)
                nd.example.put('b0', C) # -> put('b0', C.pin['b0'])
                nd.export_plt()
        """
        if drc is None:
            drc = cfg.drc
        C = [None] * 5
        for i, arg in enumerate(args):
            C[i] = arg

        active_cell = cfg.cells[-1]
        parent_cnode = active_cell.cnode
        nbs = len(parent_cnode.nb_cnode)
        if self.cnode is parent_cnode:
            main_logger(
                f"can not connect cell '{active_cell.cell_name}' to itself.",
                "error"
            )

        instance = self._copy_cell2instance(parent_cnode, flip, flop, array, scale)
        nodeI, nodeC, T = self.parse_instance_connection(*C, instance=instance)
        if flop:
            T = (T[0], T[1], T[2]+180)
        connect_geo(nodeC, nodeI, T, drc=drc)

        if cfg.solve_direct:
            # solve position of cnodeI via 'mat' of solved node nodeI
            mat = nodeI.nb_geo[0][1]
            instance.cnode.pointer.set_mat(nodeI.pointer.trans(mat))
            # solve all instance nodes using cnodeI (pp[0])
            for pp in instance.cnode.nb_geo:
                pp[0].pointer.set_mat(instance.cnode.pointer.trans(pp[1]))

#            if False: # print connection info for debugging.
#                for name, P in instance.pin.items():
#                    ni = ''
#                    if P == nodeI:
#                        ni = "*"
#                    print(" {}{} nb:{}".format(ni, name, len(P.nb_geo)))
#                    for pp in P.nb_geo:
#                        if pp[0] == instance.cnode:
#                            print('    cnodeI', pp[0].name)
#                        elif pp[0] == nodeI:
#                            print('    nodeI', pp[0].name)
#                        else:
#                            print('    other', pp[0].name)

        # Calculate bbox position of the instance in the parent

        # p2p connect:
        if self.pin2 is not None and not cfg.group_connect:
            if flip:
                msg = f"Can not close the p2p connection with flip is True in cell '{self.cell_name}'."
                netlist_logger(msg, "error")
            elif C[0] is not None:
                msg = f"Can not close the p2p connection with non-empty put() call on cell '{self.cell_name}'."
                netlist_logger(msg, "error")
            else:
                connect_geo(instance.pinout, self.pin2, diff(instance.pinout, self.pin2), solve=False, drc=drc)

        # add the new instance bbox information to the parent
        parent = parent_cnode.cell
        newi = parent_cnode.nb_cnode[nbs][0]  # new instance node
         # TODO: check here for auxiliary?
        self._calculate_instance_bbox(newi)
        ibboxloc = newi.instance.bbox

        # group connect
        if cfg.group_connect and not instance.cell.auxiliary:

            # create list of pins in the parent (self) that are available for connections:
            names = []
            for name, pin in parent.pin.items():
                # parent pins:
                if pin.type in [0, 1]:
                     names.append(f"{name} - {pin.type}")
            if len(names) > 0:
                pass
                #print(names)

            # select pins in new instance available for connections:
            newipins = {}
            for name, pin in newi.instance.pin.items():
                # TODO: pre-store eligible pins (types)
                if pin.type != 'bbox':
                    newipins[name] = pin

            # Note that adding DRC items to self may change parent.ibbox in the iteration if not auxiliary cells
            ibboxes = parent.ibbox.copy()
            for oldinode, oldibbox in ibboxes.items(): # does not yet include the new instance
                if oldinode.instance.cell.auxiliary:
                    continue
                if ibboxloc == [] or ibboxloc == [] or ibboxloc is None or oldibbox is None:
                    #print(f"SKIP: {parent_cnode.cell.cell_name:30}")
                    continue
                if ibboxloc[2] < oldibbox[0]-1:
                    continue
                if ibboxloc[0] > oldibbox[2]+1:
                    continue
                if ibboxloc[3] < oldibbox[1]-1:
                    continue
                if ibboxloc[1] > oldibbox[3]+1:
                    continue

                # instance pins:
                # create dict of existing instance pins available for connections
                # TODO: store this to not regenerate every time?
                opindict = {}
                for name, pin in oldinode.instance.pin.items():
                     # if pin.type in [0, 1]:
                    if pin.type != 'bbox':
                        opindict[name] = pin
                        #names.append(f"{name} - {pin.type}")
                if len(names) >= 0:
#                    print(
#                        f"  P:{parent_cnode.cell.cell_name:20}"\
#                        f" Inew:{instance.id:4}"\
#                        f" {instance.name:20}"\
#                        f" Iold:{inode.instance.id:4}"\
#                        f" {inode.instance.name:20}"\
                       # f"  {inode.instance.bbox}"\
                       # )
                    #print("  oldI:", old.keys())

                    for oname, opin in opindict.items():
                        if opin.xs is not None:  # do not scan further for pins with xs = None
                            for nname, npin in newipins.items():
                                if npin.xs is not None:  # do not scan further for pins with xs = None
                                    if nodeI.up.name == nname and nodeC.instance is not None:
                                        continue
                                    if nname != 'org' and oname != 'org':
                                        x, y, a = diff(opin, npin)
                                        d = hypot(x, y)
                                        if d <= cfg.autoconnectdistance:
                                            #print(f"{nodeI.up.name} {nodeC.instance} {oname} - {nname}: {d:0.1f}")
                                            # may include a new put statement
                                            connect_geo(opin, npin, (x, y, a), solve=False, drc=drc)
                pass
        parent.ibbox[newi] = ibboxloc

        # ribbon connect:
        if nodeC.type is not None and nodeI.type is not None and not cfg.group_connect:
            if nodeC.type in [3, 4] and nodeI.type in [3, 4]:
                #print("Ribbon connection detected!")
                N0 = nodeC.cnode.cell.properties.get("N", 0)
                N1 = nodeI.cnode.cell.properties.get("N", 0)

                assert (N0 == N1),\
                    "Can not connect two ribbons with different amounts of waveguides ({}, {}).".format(N0, N1)
                assert N0 != 0, "Error: Ribbon has N=0"

                pitch0 = nodeC.cnode.cell.properties.get("pitch", None)
                pitch1 = nodeI.cnode.cell.properties.get("pitch", None)
                assert (pitch0 == pitch1),\
                    "Can't connect two ribbons with different pitches ({}, {}).".format(pitch0, pitch1)
                assert pitch0 != 0, "Error, ribbon pitch=0"

                dx, dy, da = diff(nodeC, nodeI)
                dis = hypot(dx, dy)
                dismax = 10e-11
                damax = 0,.001
                if (dis < dismax) or (180-damax < da < 180+damax):
                    _connect_ribbon(nodeC, nodeI, N0)
                else:
                    print("WARNING: pins do not smoothly connect: distance {:0.3f} > {:0.3f} and/or angle |{:0.4f} -180|;> {:0.4f}.".\
                        format(dis, dismax, da, damax))

        # update cp:
        p = kwargs.get('newcp', None)
        if p is None:
            p = kwargs.get('cp', None) # backward compatible
        if p is None:
            cfg.cp = instance.pin[self.default_out]
        else:
            cfg.cp = instance.pin[p]

        # for ID in trace.trace_id_list:
        if trace.trace_id_list:
            trace.trace_append(instance)

        return instance


    def rebuild(self, instantiate=True, flat=False, layermap=None,
        layermapmode=None, cellmap=None, infolevel=0, clear=False):
        """Flatten, rename, relayer, reshape and/or filter a Nazca Cell object.

        The original cell(tree) remains unchanged.
        This method operates at mask-element level. It contructs the output
        cell(tree) as geometrical info from the ground up, i.e.
        any circuit-type netlist information will not be copied.
        Rebuild can be used, for example, to place
        each layer of a cell at a different position in the mask layout.

        Args:
            instantiate (bool): instantiate setting of returned (top)cell.
                (default=True)
            flat (bool): flatten cell(tree) (default=False)
            layermap (dict): {oldlayer, newlayer} mapping
            cellmap (dict): to be implemented
            layermapmode ('all' | 'none'): start mapping with all layers
                included in the map: 'all', or an empty map 'none'.
                (default='all')

        Returns:
            Cell: rebuilt input cell
        """
        return layout.rebuild(
            cell=self,
            instantiate=instantiate,
            flat=flat,
            layermap=layermap,
            layermapmode=layermapmode,
            clear=clear,
            infolevel=infolevel,
        )


    def __str__(self):
        """Print which structures are contained in the cell."""
        s = 'cell info of \'{}\''.format(self.cell_name)
        s += '\n--pin     {:2}x              = {}'.\
            format(len(self.pin.keys()), self.pin.keys())
        s += '\n--polygon {:2}x (layer, N)   = '.format(len(self.polygons))

        L = []
        for g in self.polygons:
            L.append('({g[1]}, {N})'.format(g=g, N=len(g[1].points)))
        s += ', '.join(L)
        s += '\n--gds     {:2}x (file, cell) = '.format(len(self.gdsfiles))

        L = []
        for pointer, g in self.gdsfiles:
            L.append('({g[1]}, {g[2]})'.format(g=g))
        s += ', '.join(L)

        L = []
        #s += '\n--inst    {:2}x (x,y, cell)  = '.format(len(self.instances))
        #for i in self.instances:
        #    L.append('({}, {})'.format(i[0], i[1].cell_name))
        #s += ', '.join(L)
        return s + '\n'


    def _polyshape_iter(self, trans=Pointer(0), flip=False, scale=1.0):
        """Generator to iterate over all polylines and polygons.

        Yields:
            Polygon | Polyline, list(float, float): poly object, transformed points
        """
        pgon_iter = Netlist().polygon_transflip_iter(
            self.cnode,
            trans,
            flip=flip,
            scale=scale,
            apply=True,
            hierarchy='apply',
            hull=True,
        )
        for pgon, xy, bbox in pgon_iter:
            if pgon.layer not in cfg.bbox_layers_ignore:
                yield pgon, xy

        pline_iter = Netlist().polyline_transflip_iter(
            self.cnode,
            trans,
            flip=flip,
            scale=scale,
            apply=True,
            hierarchy='apply',
            hull=True,
            polygon=True,
        )
        for pline, xy, bbox in pline_iter:
            if pline.layer not in cfg.bbox_layers_ignore:
                yield pline, xy


    def _calculate_instance_bbox(self, inode):
        """Calculate bbox of the instance inode with respect to its parent.

        inode (Node): instance node of self

        Returns:
            (float, float, float, float): bbox (x1, y1, x2, y2)
        """
        org = inode.pointer.copy()
        if inode.flip:
            s = -1
        else:
            s = 1
        x0, y0, a0 = org.xya()
        a = np.radians(a0)
        S = inode.scale

        # add hull_points
        bbhull = []
        basepoints = []
        hullarrs = []
        if self.hull_based_bbox:
            if self.hull is not None:
                basepoints = self.hull
        if len(basepoints) == 0:
            basepoints = self.bbox_polygon
            # this may validly be an empty bbox: []
            #if len(basepoints) == 0:
            #    print(f"Empty basepoints in {self.cell_name}'")

        for x, y in basepoints:
            bbhull.append((
                x0 + S * (cos(a)*x - s*sin(a)*y),
                y0 + S * (sin(a)*x + s*cos(a)*y)
            ))

        # add points in case of array's
        if inode.array is not None:
            hullarrs = []
            Nx, (x1, y1), Ny, (x2, y2) = inode.array
            Dx1 = (Nx-1)*x1
            Dy1 = (Nx-1)*y1
            Dx2 = (Ny-1)*x2
            Dy2 = (Ny-1)*y2
            for x, y in bbhull:
                hullarrs.extend([
                    (x, y),
                    (Dx1 + x, Dy1 + y),
                    (Dx2 + x, Dy2 + y),
                    (Dx1 + Dx2 + x, Dy1 + Dy2 + y)
                ])
        else:
             hullarrs = bbhull

        if len(hullarrs) > 0:
            xy = list(zip(*hullarrs))
            inode.instance.bbox = (min(xy[0]), min(xy[1]), max(xy[0]), max(xy[1]))
        inode.instance.hullarrs = hullarrs
        return hullarrs


    def _calculate_bbox(self, instance=True):
        """Calculate the bounding box, userbox and hull of a cell.

        Args:
            instance (bool): instantiation flag

        The bounding box (bbox) is based on the cell content and rectangular.
        The userbbox is a polygon and can have any number of points > 2.
        The hull is the smallest convex polygon that still encloses all the
        points in the cell.

        The userbbox can by choice either be part of the bbox and hull calculation or
        stay separate and possible extend beyond the bbox.

        The bbox can be calculated based on the hull of its instances,
        if it exists, or based on the bbox of the instance.

        Returns:
            None
        """
        self.bbox_complete = True # flag to indicate if all bbox info was availble
        limit = 1e8
        self.bbox = (limit, limit, -limit, -limit)

        elements = 0
        hullarrs = []
        if self.userbbox is not None and self.userbbox_include:
            elements += 1
            hullarrs.append(self.userbbox)

        # add polygon dimensions
        for polyobj, xy in self._polyshape_iter():
            elements += 1
            hullarrs.append(xy)

        # add instance its bbox dimensions to hullarrs:
        for inode in self.cnode.instance_iter(auxiliary=False):
            elements += 1
            ha = inode.instance.hullarrs
            if ha:
                hullarrs.append(ha)

        if not hullarrs:
            self.bbox = None
            self.bbox_polygon = []
            #logger.warning("Empty bounding box for cell '{}'".format(self.cell_name))
            return None

        else:
            cat_hullarrs = np.concatenate(hullarrs)

        hull_created = False
        if self.use_hull: # make a hull for this cell (self)
            try: # ConvexHull may fail
                H = ConvexHull(cat_hullarrs)
                self.bbox = (
                    min(self.bbox[0], H.min_bound[0]), min(self.bbox[1], H.min_bound[1]),
                    max(self.bbox[2], H.max_bound[0]), max(self.bbox[3], H.max_bound[1]))
                self.hull = np.take(H.points, H.vertices, axis=0)
                if self.show_hull:
                    Polygon(self.hull, layer=cfg.default_layers['hull']).put(0)
                hull_created = True
            except Exception as excp:
                self.hull = None
                hull_created = False
                main_logger(
                    f"No hull resolved in cell '{self.cell_name}'. Using square bbox instead.\n{excp}",
                    "exception"
                ) # probably 1D polygons, like strt(length=0)
        if not hull_created:  # if no convex hull has been created, make a normal bbox
            xy = list(zip(*cat_hullarrs))
            self.bbox = (min(xy[0]), min(xy[1]), max(xy[0]), max(xy[1]))

        if elements == 0:
            self.bbox_complete = False
            self.bbox = None
            self.bbox_polygon = []
        else:
            self.bbox = (
            self.bbox[0] - self.bboxbuf,
            self.bbox[1] - self.bboxbuf,
            self.bbox[2] + self.bboxbuf,
            self.bbox[3] + self.bboxbuf)
            #print("\n-- bbox: '{}': ({:.3f}, {:.3f}, {:.3f}, {:.3f})".\
            #    format(self.cell_name, *self.bbox ))
            self.bbox_polygon = np.array([
                [self.bbox[0], self.bbox[1]],
                [self.bbox[0], self.bbox[3]],
                [self.bbox[2], self.bbox[3]],
                [self.bbox[2], self.bbox[1]]])

        area = None
        if self.bbox is not None:
           area = abs((self.bbox[2] - self.bbox[0]) * (self.bbox[3] - self.bbox[1]))

        self.add_property(
            {'dimension':
                {'bbox': self.bbox,
                 'bbox_area': area},
            }
        )
        return None


    def _add_bbox(self, draw_bbox=True):
        """Add the bbox to the cell.

        Returns:
            None
        """
        self.add_bbox_flag = True
        self._calculate_bbox()

        if draw_bbox:
            if len(self.gdsfiles) > 0:
                if self.bbox is None:
                    main_logger(
                        f"Can not determine a bounding box for cell '{self.cell_name}',"
                        "  e.g. due to only non-native GDS subcell(s) and/or empty subcell(s)."
                        "  In case of gds loading causing this message:"
                        "  Try keyword native=True in  to get the complete bbox"
                        "  or set keyword bbox=False to get rid of this message.",
                        "warning"
                    )
                else:
                    nd.logger.main(
                        f"Not all instances in cell '{self.cell_name}' have a known bbox,"
                        "  e.g. due to a non-native GDS subcell or an empty subcell."
                        "  Therefore, the bbox maybe smaller than the actual structures."
                        "  Try option native=True in gds loading to get the complete bbox"
                        "  or set bbox=False to get rid of this message.",
                        "warning"
                    )
            if self.bbox is not None:
                length = self.bbox[2] - self.bbox[0]
                width = self.bbox[3] - self.bbox[1]
                bbu.put_boundingbox(
                    'org',
                    length=length,
                    width=width,
                    move=(self.bbox[0], self.bbox[1]+0.5*width, 0),
                    params=True,
                    name_layer=self.name_layer,
                )
                self.length = length
                self.width = width
        if not cfg.solve_direct:
            Netlist().solvecell(self)


    def _store_pins(self):
        """Store pin information in annotation.

        Returns:
            None
        """
        pintxt = 'pins:\n'
        for name, node in sorted(self.pin.items()):
            if node.xs is None:
                xs = 'None'
            else:
                xs = None
            if node.width is None:
                w = 'None'
            else:
                w = "{:.3f}".format(float(node.width))
            pintxt += '{0}: {c[0]:.5f}, {c[1]:.5f}, {c[2]:.5f}, {xs}, {w}\n'.\
                format(name, xs=xs, w=w, c=node.xya())
        Annotation(text=pintxt, layer=cfg.default_layers['pin_text']).put(0)

    
    def close(self):
        """Close the cell.

        Solve the geometry of the nodes in the Cell.
        Set the cp back to the position before opening the cell.
        If default input and/or output pins have not been set yet they will be
        set on 'org', pointing in opposite directions.

        Returns:
            Cell: self
        """
        if self.closed:
            main_logger(f"Cell '{self.cell_name}' already closed.",
                "warning"
            )

        #add default ports if missing:
        if self.default_in not in self.pin:
            self._put_pin(self.default_in, (0, 0, 180))
            #print("Warning: no default pin '{}' set in '{}', setting it to {}.".\
            #    format(self.default_in, self.cell_name, (0, 0, 180)))
        if self.default_out not in self.pin:
            self._put_pin(self.default_out, (0))
            #print("Warning: no default pin '{}' set in '{}', setting it to {}.".\
            #    format(self.default_out, self.cell_name, (0, 0, 0)))
            #self.put_pin(self.default_out, cfg.cp)

        self.pinin = self.pin[self.default_in]
        self.pinout = self.pin[self.default_out]

        if not cfg.solve_direct:
            Netlist().solvecell(self)

        if not self.add_bbox_flag:
            if self.autobbox:
                self._add_bbox(draw_bbox=True)
            else:
                self._add_bbox(draw_bbox=False)
            if self.params:
                bbu.put_parameters(parameters=self.parameters)

        if self.userbbox:
            Polygon(points=self.userbbox, layer='userbbox').put(0)
        elif self.protectbox:
            Polygon(points=self.bbox_polygon, layer='userbbox').put(0)


        if self.store_pins:
            self._store_pins()

        # Add version annotation to cell.
        if cfg.store_bbname and self.version is not None:
            version = self.version.copy()
            if self.instantiate:
                if isinstance(self.version, dict):
                    if version.get('cellname', None) == 'auto':
                        version['cellname'] = self.cell_basename
                    if cfg.force_pdk_version is not None:
                        version = cfg.force_pdk_version
                    bbtxt = yaml.dump(version)
                    Annotation(text=bbtxt, layer=get_layer('bb_name')).put(0)
                else:
                    main_logger(f"In cell '{self.cell_name}' version is not of type dict: {self.version}",
                        "error"
                    )
                self.version = version
            else:
                raise Exception("Cell must be instantiated to add a name and version to it"
                    " Cell: '{self.cell_name}'"
                )
        self.closed = True

        cfg.cp = self.oldcp
        cfg.cells.pop()
        if cfg.cells:
            cfg.self = cfg.cells[-1]
        return self


def diff(node1, node2):
    """Calculate the geometrical difference between two nodes.

    Args:
        node1 (Node): start node
        node2 (Node): end node

    Example:
        Obtain dx, dy, da from p1 to p2::

            import nazca as nd

            p1 = nd.Pin().put(10, 10, 0)
            p2 = nd.Pin().put(20, 30, 45)
            dx, dy, da = nd.diff(p1, p2)
            print(dx, dy, da)
            # (10.0, 20.0, 45.0)

    Returns:
        (dx, dy, da): difference vector <node2> - <node1> between their pointers
    """
    if not cfg.solve_direct:
        cfg.cells[-1]._solve()
    ptr1 = node1.pointer.copy()
    ptr2 = node2.pointer.copy()
    ptr2.multiply_ptr(ptr1.inv())
    if not cfg.solve_direct:
        cfg.cells[-1].closeflag = False
    #print('nd.diff:', ptr2.get_xya())
    dx, dy, da = ptr2.get_xya()
    if da == 360:
        da = 0
    return dx, dy, da


def midpoint(node1, node2):
    """Calculate the geometrical center between two nodes.

    See also function midpointer.

    Args:
        node1 (Node): start node
        node2 (Node): end node

    Returns:
        (x, y, a): mid-point between pointers of <node1> and <node2>

    Example::

        import nazca as nd

        p1 = nd.Pin().put(10, 20, 0)
        p2 = nd.Pin().put(20, 40, 90)

        x, y, a = nd.midpoint(p1, p2)
        print(x, y, a)
        # 15.0, 30.0, 45.0
    """
    x1, y1, a1 = node1.xya()
    x2, y2, a2 = node2.xya()
    da = a2 - a1
    if da > 180:
        da -= 360
    a = round(a1 + da /2, 12) # avoid small negative values
    if a < 0:
        a += 360
    return x1+(x2-x1)/2, y1+(y2-y1)/2, a


def midpointer(node1, node2):
    """Calculate the geometrical center between two nodes.

    See also midpoint.

    Args:
        node1 (Node): start node
        node2 (Node): end node

    Returns:
        Pointer: pointer in mid-point between pointers of <node1> and <node2>

    Example::

        import nazca as nd

        p1 = nd.Pin().put(10, 20, 0)
        p2 = nd.Pin().put(20, 40, 90)

        p3 = nd.midpointer(p1, p2)
        print(p3.xya())
        # (15.0, 30.0, 45.0)
    """
    x1, y1, a1 = node1.xya()
    x2, y2, a2 = node2.xya()
    da = a2 - a1
    if da > 180:
        da -= 360
    a = round(a1 + da /2, 12) # avoid small negative values
    if a < 0:
        a += 360
    return Pointer(x1+(x2-x1)/2, y1+(y2-y1)/2, a)


def bbinfo(item=None, data=None):
    """Attach metadata to the Active cell in dictionary bbinfo.'

    Args:
        item (hashable object): dictionary key
        data: The value of the metadata

    Returns:
        None
    """
    cfg.cells[-1].bbinfo[item] = data


def get_transformation(start=None, end=None):
    """Return the transformation matrix that connects two pointers.

    Returns:
        float[:]: translation matrix
    """
    translate = dot(inverse(start.pointer.mat), end.pointer.rotate(180).mat)
    return translate


#==============================================================================
# Class structures to represent mask structures that can be put.
#==============================================================================
class Pin():
    """Pin class"""
    def __init__(
        self,
        name=None,
        xs=None,
        layer=None,
        width=None,
        radius=None,
        type=None,
        pin=None,
        io=None,
        flip=False,
        show=None,
        remark=None,
        chain=None
    ):
        """Construct a Pin object.

        A Pin object (capital P to denote the class object) holds information
        to describe a Node object. Nodes are the objects that are connected in
        a graph to connect layout elements. The Node object is created inside
        the active cell via the put() method of the Pin.
        When a Node is created in the active cell it is stored in
        the pin attribute of the active cell (not to be confused with the Pin
        object). The pin attribute is a dictionary
        with the string name of the Pin object as key and the Node object as
        the value.

        Note that Pin represenation in a layout is defined by the pinstyle
        of the xsection of the pin. The pinstyle maps the xs to the
        visualisation layer.

        Args:
            name (str): name of the Pin to refer to it.
            xs (str): xsection name assigned to the pin.
            layer (str or number): optional layer to represent the pin.
                This will add an xs pinstyle. If no xs is provided it will use
                the default xs. Note that layer is not a pin property.
            width (float): width in um of the connection the pin represents
            radius (float): radius of curvature at the pin in um.
            type (str): extra property for the pin
            pin (Node):
            chain (int):indicate if pin is a chain or connector (1) or not (0)
                (default=1)
            remark (str): short optional documentation of the Pin (default=None)
            io (int): connectivity info in bit form
            flip (bool): symmetry state of the pin. This flip is *only* relevant
                for stub representation and DRC checks of assymmetric interconnects,
               default=False.

        Returns:
            None

        Example:
            Create a Pin object and put it in the active cell to create a Node.
            Subsequently, retreive the Node by name 'a0' via cell pin attribute
            and print properties as set in the Pin object::

            import nazca as nd

            nd.Cell(name='test') as C:
                nd.Pin(name='a0', width=2.0).put()

            p = C.pin['a0']
            print(p.name)
            print(p.width)
            # 'a0'
            # 2.0
        """
        # set defaults:
        self.xs = None
        self.width = None
        self.radius = 0
        self.name = 'noname'
        self.type = None
        self.chain = 1
        self.show = False
        self.remark = remark
        self.io = 0
        self.pin = pin

        # copy pin properties if a pin is provided:
        if pin is not None:
            if pin.cnode.instance is False:  # check if not by mistake referencing to a cell
                if pin.cnode is not cfg.cells[-1].cnode:  # allow a copy from pin in same call
                    interconnect_logger(
                        msg=f"Connecting to cell pin {cfg.cells[-1].cell_name}.pin['{pin.name}'] rather than an instance of this cell.",
                        level="warning"
                    )
            self.name = pin.name
            self.xs = pin.xs
            self.width = pin.width
            self.radius = pin.radius
            self.io = pin.io
            self.type = pin.type
            self.chain = pin.chain
            self.remark = pin.remark

        if name is not None:
            self.name = name

        # override pin properties if set explicitly:
        if xs is not None:
            # TODO: Make it possible to overrule the pin with a None value (for DRC)?
            self.xs = xs
            if pin is None:
                # take default width setting from xsection
                XS = cfg.XSdict.get(xs, None)
                if XS is not None:
                    self.width = XS.width
                    if layer is not None: # add a new layer (only for new pinstyles)
                        layer = get_layer(layer)
                        presentstyle = cfg.pinstyles.get(xs, None)
                        if presentstyle is None:
                            newstyle = cfg.pinstyles['default'].copy()
                            newstyle['layer'] = layer
                            bbu.add_pinstyle(xs, newstyle)
                            XS.pinstyle = xs
                        else:
                            main_logger(f"In setting pin '{cfg.cells[-1].cell_name}'.'{self.name}',"
                                f" can not add layer {layer} to existing pinstyle '{xs}'"
                                f" having already layer {presentstyle['layer']}."
                                f" Alternatively, set an explicit with bb_util.add_pinstyle() for xs '{xs}'.",
                                "error"
                            )
        if width is not None:
            self.width = width
        if type is not None:
            self.type = type
        if radius is not None:
            self.radius = radius
        if chain is not None:
            self.chain = chain
        if show is not None:
            self.show = show
        if remark is not None:
            self.remark = remark
        if flip:
            self.io = 1  # store flip in io.
        if io is not None:
            self.io = io


    def __repr__(self):
        return "<Pin() object, name='{}', xs='{}', width={}, type={}, "\
            "chain={}, show={}, io={}, remark={}>".\
            format(self.name, self.xs, self.width, self.type, self.io, self.chain,
                self.show, self.io, self.remark)


    def put(self, *args):
        """Put a Pin object in the layout."""
        if self.pin is not None:
            cfg.cp = self.pin
        return cfg.cells[-1]._put_pin(self.name, connect=args, xs=self.xs,
            width=self.width, type=self.type, io=self.io, chain=self.chain, show=self.show,
            remark=self.remark)


class Polygon():
    """Polygon class."""
    def __init__(self, points=None, layer=None):
        """Construct a Polygon object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int | (int, int) | str): layer number, (layer, datatype) tuple or layer name to put the Polygon in.
            points (list): list of points [(x1, y1), (x2, y2), ...]

        Returns:
            None
        """
        self.layer = get_layer(layer)
        if points is None:
            self.points = []
        else:
            self.points = points
        hull = False  # flag
        if cfg.use_hull:
            try:
                hull = True
                H = ConvexHull(points)
                self.hull = np.take(H.points, H.vertices, axis=0)
                self.bbox = [H.min_bound[0], H.min_bound[1], H.max_bound[0], H.max_bound[1]]
            except Exception as excp:
                #less than 2D polygon, e.g. a zero length straight
                main_logger(
                    f"CovexHull Error Polygon: {excp}\n{points}.",
                    "error"
                )
                hull = False
        if not hull:
            xy = list(zip(*points))
            self.bbox = (min(xy[0]), min(xy[1]), max(xy[0]), max(xy[1]))


    def __repr__(self):
        size = len(self.points)
        return "<Polygon() object, layer={}, points={}, bbox={}>".\
            format(self.layer, size, self.bbox)


    def transform(self, center=(0, 0, 0), scale=1.0,
        flipx=False, flipy=False, move=(0, 0, 0), layer=None):
        """Transform the points in the Polygon object.

        See nazca.geometries.transform().

        Returns:
            Polygon: new Polygon with transformed points.
        """
        if layer is None:
            layer = self.layer
        points = transform(points=self.points, center=center, scale=scale,
            flipx=flipx, flipy=flipy, move=move)
        newPolygon = Polygon(points=points, layer=layer)
        return newPolygon


    def grow(self, grow=5, accuracy=0.1, jointype='square', layer=None):
        """Grow the polygon.

        Returns:
            Polygon: new Polygon with growth applied.
        """
        if layer is None:
            layer = self.layer
        points = grow_polygons(paths=[self.points], grow=grow, accuracy=accuracy,
            jointype=jointype)
        newPolygon = Polygon(points=points[0], layer=layer)
        return newPolygon


    def put(self, *args, flip=False, flop=False, scale=1.0):
        """Put a Polygon object in the layout.

        Args:
            flip (bool): flip state of the polygon put
            flop (bool): flop state of the polygon put
            scale (float): scaling factor of the Polygon. This scaling looks
                as if the whole cell the polygon would have scaled.
                Note that it may be more appropiate to use
                the 'transform' method to scale a Polygon.

        Returns:

        """
        return cfg.cells[-1]._put_polygon(connect=args, flip=flip, flop=flop,
            scale=scale, polygon=self)


class Polyline():
    """Polyline (path) class."""

    def __init__(self, points=None, layer=None, width=None, pathtype=0, miter=0.5):
        """Construct a Polyline object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int | (int, int) | str): layer number, (layer, datatype) tuple or layer name to put the Polygon in.
            width (float): width of the polyline
            pathtype (int): gds type of path: 0 (flush), 1 (round), 2, 3
            points (list): list of points [(x1, y1), (x2, y2), ...]
            miter (float): length at joint of two straight sections (default=0.5)

        Returns:
            None
        """
        self.layer = get_layer(layer)
        if points is None:
            self.points = []
        else:
            self.points = points
        if width is not None:
            self.width = float(width)
        else:
            self.width = 0.20
        self.pathtype = pathtype
        self.points = points
        self.miter = miter
        self.polygon = util.polyline2polygon(self.points, width=self.width,
                miter=self.miter)
        hull = False
        if cfg.use_hull:
            try:
                hull = True
                H = ConvexHull(self.polygon)
                self.hull = np.take(H.points, H.vertices, axis=0)
                self.bbox = [H.min_bound[0], H.min_bound[1], H.max_bound[0], H.max_bound[1]]
            except Exception as excp:
                #less than 2D polyline, e.g. a zero length straight
                main_logger(f"CovexHull Error Polygon: {excp}\n{points}.",
                    "error"
                )
                hull = False
        if not hull:
            xy = list(zip(*self.polygon))
            self.bbox = (min(xy[0]), min(xy[1]), max(xy[0]), max(xy[1]))


    def __repr__(self):
        size = len(self.points)
        return "<Polyline() object, layer={}, width={}, pathtype={}, "\
            "points={}, bbox={}>".\
            format(self.layer, self.width, self.pathtype, size, self.bbox)


    def transform(self, center=(0, 0, 0), scale=1.0,
        flipx=False, flipy=False, move=(0, 0, 0), width=None, pathtype=None,
            layer=None):
        """Transform the points in the Polyline object.

        See nazca.geometries.transform().

        Returns:
            Polyline: new Polyline with transformed points.
        """
        if layer is None:
            layer = self.layer
        if width is None:
            width = self.width
        if pathtype is None:
            pathtype = self.pathtype
        points = transform(points=self.points, center=center, scale=scale,
            flipx=flipx, flipy=flipy, move=move)
        newPolyline = Polyline(points=points, width=width, pathtype=pathtype,
            layer=layer)
        return newPolyline


    def grow(self, grow=5, accuracy=0.1, jointype='square', width=None,
            layer=None):
        """Grow the polyline.

        Returns:
            Polyline: new Polyline with growth applied.
        """
        if layer is None:
            layer = self.layer
        points = grow_polygons(paths=[self.points], grow=grow, accuracy=accuracy,
            jointype=jointype)
        newPolyline = Polyline(points=points[0], layer=layer)
        return newPolyline


    def put(self, *args, flip=False, flop=False, scale=1.0):
        """Put a Polyline object in the layout.

        Args:
            flip (bool): flip state of the polygon put
            flop (bool): flop state of the polygon put
        scale (float): scaling factor of the Polygon. This scaling looks
                as if the whole cell the polygon would have scaled.
                Note that it may be more appropiate to use
                the 'transform' method to scale a Polygon.

        Returns:

        """
        return cfg.cells[-1]._put_polyline(connect=args, flip=flip, flop=flop,
            scale=scale, polyline=self)


class Annotation():
    """Annotation class."""
    def __init__(self, layer=None, text=''):
        """Construct an annotation object.

        Args:
            layer (int | (int, int) | str): layer number, (layer, datatype) tuple or layer name to put the Polygon in.
            text (str): annotation text (default = '')

        Returns:
            None
        """
        self.layer = get_layer(layer)
        self.text = text

    def __repr__(self):
        if len(self.text) >= 10:
            text = self.text[:7]+'...'
        else:
            text = self.text
        return "<Annotation() object, layer={}, text='{}'>".format(
            self.layer, text)

    def put(self, *args):
        """Put an Annotation object in the layout."""
        return cfg.cells[-1]._put_annotation(connect=args, annotation=self)


def _scan_branch(strm, cellname, celllist=None, level=0):
    """Generator for cell names in a branch. Bottom-up, deep first."""

    cell_rec = strm.cells[cellname]
    snames = cell_rec.snames
    #print("cell:{}, snames:{}".format(cellname, snames))
    if snames is not None:
        level += 1
        for sname in snames:
            yield from _scan_branch(strm, sname, level=level)
        level -= 1
    yield(cellname)


def _string_to_pins(string):
    """Convert a pin annotation string to a pin dictionary.

    Expected format of <string>:

    "pins:
    <name>: <x, y, a xs, w>
    <name>: <x, y, a xs, w>
    ..."

    Header 'pins:' is case incensitive and spaces will be stripped.

    Returns:
        ordered dict: {<parameter_name>: <parameter_value>}
    """
    pins = None
    lines = string.split('\n')
    pins = OrderedDict()
    if (lines[0].lower() == 'pins:'):
        for line in lines[1:]:
            if line == '':
                continue
            name, prop = line.split(':', 1)
            x, y, a, xs, w = prop.replace(' ', '').split(',')
            x = float(x)
            y = float(y)
            a = float(a)
            if w.lower() != 'none':
                w = float(w)
            pins[name.strip()] = (x, y, a, xs, w)
    return pins


def _gds2native(
    strm,
    topcellname=None,
    cellmap=None,
    bbox=False,
    bboxbuf=0,
    hull=None,
    show_hull=False,
    hull_based_bbox=None,
    scale=None,
    flat=False,
    instantiate=True,
    asdict=False,
    topcellsonly=True,
    select='',
    reuse=False,
    cellsnotreused=None,
    cellsreused=None,
):
    """Translate a GDS stream into native native Nazca cell structure.

    If no cellname is provided it will select the topcell if there is only
    one topcell.

    Args:
        strm (bytestr): gds layout
        topcellname (str): name of topcell to process
        cellmap (dict): cellname (and under) to process
        bbox (bool): Add bounding box if True (default=False)
        bboxbuf (float): extra buffer space added to the bbox (default=0)
        scale (float): scaling of the cells (default=1.0)
        flat (bool): set instantiation off all subcells (default=False)
        instantiate (bool): instantiation setting of the topcell (default=True)
        reuse (bool): default is False. If True, load gds assumes cells with the
            same name are the same cell and it resolves conficts by using
            only the first cell occurance. Setting reuse=True is dangerous.
        cellsnotreused (list of str): list of cellnames to ignore for
            reuse if reuse=True
        cellsreused (list of str): list of cellnames to reuse. If set all other cells will
            not be reused and reuse and cellsnotreused will be ignored.

    Returns:
        list of str: list of cellnames in <filename> under cell name <cellname>
    """
    annotated_cells = []
    um0 = 1e-6 / gb.gds_db_unit
    um = um0 / scale # scale factor between gds integers and micro-meters.
    gds_elements_unknown = set()

    NC = {}
    rename_map = {}

    if cellsnotreused is None:
        cellsnotreused = []
    if cellsreused is None:
        cellsreused = []
    if cellsreused:
        reuse = True
    if reuse:
        if cellsreused:
           cells_forreuse = set(cellsreused) - set(cellsnotreused)  # only reuse explicit cells from the list for reuse and cells in this load call
        else:
           cells_forreuse = set(cfg.cellnames.keys()) - set(cellsnotreused)   # reuse all cells in memory
    else:
        cells_forreuse = []

    cells_visited = set()

    if asdict:
        topcells = strm.topcell()
    else:
        topcells = [topcellname]

    for topcellname in topcells:
        branch_iter = _scan_branch(strm, topcellname)
        C = None  # stays None if only resued cells are encountered.
        for cellname in branch_iter: # bottom up
            if cellname in cells_visited:
                continue
            if reuse:
                if cellname in set(cfg.cellnames):
                    if cellname in cells_forreuse:  # cells visited during program execution
                        main_logger(f"reusing cell '{cellname}' when loading file '{cfg.gdsload}'", "info")
                        continue

            cells_visited.add(cellname)
            cellrecord = strm.cells[cellname]

            if cellname == topcellname:
                _instantiate = instantiate
            else:
                if flat:
                    _instantiate = False
                else:
                    _instantiate = True
            #print(cellname)
            with Cell(
                    cellname,
                    instantiate=_instantiate,
                    hull=hull,
                    show_hull=show_hull,
                    hull_based_bbox=hull_based_bbox
                ) as C:
                C.default_pins('org', 'org')
                for elem in cellrecord.elements:

                    if elem.etype == gb.GDS_record.BOUNDARY:
                        LD, XY = elem.polygon
                        LD = strm.layermap.get(LD, LD)
                        XY = [(x/um, y/um) for x, y, in zip(*[iter(XY)]*2)]
                        Polygon(points=XY, layer=LD).put(0)

                    elif elem.etype == gb.GDS_record.PATH:
                        LD, XY, pathtype, width = elem.polyline
                        LD = strm.layermap.get(LD, LD)
                        if width is not None:
                            width = width/um
                        if pathtype is None:
                            # TODO: add debug?
                            pathtype = 0
                        XY = [(x/um, y/um) for x, y, in zip(*[iter(XY)]*2)]
                        Polyline(points=XY, layer=LD, pathtype=pathtype,
                            width=width).put(0)

                    elif elem.etype == gb.GDS_record.TEXT:
                        LD, XY, TEXT = elem.annotation
                        LD = strm.layermap.get(LD, LD)
                        XY = (XY[0]/um, XY[1]/um)
                        Annotation(layer=LD, text=TEXT).put(XY)
                        #TODO: also check for the layer
                        if TEXT[:4].lower() == 'pins':
                            pins = _string_to_pins(TEXT)
                            for name, (x, y, a, xs, w) in pins.items():
                                if name == 'org':
                                    continue
                                if isinstance(xs, str):
                                    if xs.lower() == 'none':
                                        xs = None
                                if isinstance(w, str):
                                    if w.lower() == 'none':
                                        w = None
                                if cellname == topcellname:
                                    x, y = x*scale, y*scale
                                Pin(name=name, xs=xs, width=w).put(x, y, a)
                                #print(name, x, y, a, xs, w)
                        if LD == cfg.default_layers[cfg.bb_version_layer]:
                            Annotation(layer=LD, text=TEXT).put(XY)
                            annotated_cells.append(C)
                            try:
                                bbname, bbversion = TEXT.split('\n')
                                C.annotated_name = bbname    # store cell name
                                C.annotated_version = bbversion # store version
                                #print('BB ID found:\n{}'.format(TEXT))
                            except ValueError as e:
                                if not isinstance(e, ValueError):
                                    logger.exception("e")
                                #print("ERROR: (netist): No valid black box"\
                                #    " annotation in layer {}: {}".format(LD, e))

                    elif elem.etype == gb.GDS_record.BOX:
                        gds_elements_unknown.add('BOX')

                    elif elem.etype == gb.GDS_record.SREF:
                        name, trans, mag, angle, XY = elem.instance
                        if trans == 1:
                            flip = True
                        else:
                            flip = False
                        name = rename_map.get(name, name)
                        cell = NC.get(name, cfg.cellnames[name])
                        cell.put('org', XY[0]/um, XY[1]/um, angle, flip=flip, scale=mag)

                    elif elem.etype == gb.GDS_record.AREF:
                        name, trans, mag, angle, col, row, XY = elem.array
                        if trans == 1:
                            flip = True
                        else:
                            flip = False
                        x, y = XY[0]/um, XY[1]/um
                        x1, y1 = XY[2]/um, XY[3]/um
                        x2, y2 = XY[4]/um, XY[5]/um
                        dx1, dy1 = (x1-x)/col, (y1-y)/col
                        dx2, dy2 = (x2-x)/row, (y2-y)/row
                        name = rename_map.get(name, name)
                        cell = NC.get(name, cfg.cellnames[name])
                        cell.put('org', x, y, angle,
                            array=(col, (dx1, dy1), row, (dx2, dy2)),
                            flip=flip, scale=mag)
                    else:
                        gds_elements_unknown.add(elem.etype)
                if cellname == topcellname:
                    # only apply bbox and/or hull to topcell
                    if bbox:
                        C.autobbox = True
                        C.bboxbuf = bboxbuf
                    if hull:
                        C.hull_based_bbox = True

            rename_map[cellname] = C.cell_name
            NC[C.cell_name] = C  # 'cellname' after cell generation

    if gds_elements_unknown and cfg.show_unknown_GDS_records:
        elist= ["{}({})".format(gb.GDS_record.name[e], e) for e in gds_elements_unknown]
        strlist = ", ".join(elist)
        main_logger(
            f"gds2nazca: Ignoring elements in cell (branch) '{cellname}': {strlist}",
            "warning"
        )
    cfg.annotated_cells = annotated_cells
    if C is None:
        #print(f"Only reused cells found when loading {topcells}")
        C = cfg.cellnames[topcellname]
    if not asdict:
        return C
    else:
        NCfilter = dict()
        if topcellsonly:
            renamed_topcells = {name:rename_map[name] for name in topcells}
            for cellname, cell in NC.items():
                if cellname in renamed_topcells.values():
                    NCfilter[cellname] = cell
        else:
            NCfilter = NC
        return NCfilter


def gds_filter(filename, cellmap=None, layermap=None):
    """Filter layers and cells from gds and export a copy to file."""
    #TODO filter out cellnames by pattern.
    print('Filtering \'{}\'...'.format(filename))
    strm = gdsimp.GDSII_stream(filename=filename, cellmap=cellmap, layermap=layermap)
    newfilename = '{}_filtered.gds'.format(filename[:-4])
    strm.GDSII_write(newfilename)
    print('...Wrote \'{}\''.format(newfilename))


def print_structure(filename, cellname, level=0):
    """Print gds structure in ascii art.

    Returns:
        None
    """
    strm = gdsimp.GDSII_stream(filename)
    strm.print_structure(cellname, level)


def load_gds(
    filename,
    cellname=None,
    newcellname=None,
    asdict=False,
    topcellsonly=True,
    select='top',
    layermap=None,
    cellmap=None,
    scale=None,
    prefix='',
    instantiate=True,
    native=True,
    bbox=False,
    bboxbuf=0,
    hull=None,
    show_hull=False,
    hull_based_bbox=None,
    connect=None,
    flat=False,
    reuse=False,
    cellsnotreused=None,
    cellsreused=None,
    layermapmode=None,
):
    """Load one or more GDS cells (and its instances) from <filename> into a Nazca cell.

    By default, load_gds returns a single cell in the specified gds file in
    <filename>, as specified in cellname, or the topcell, if only one exists.

    Alternatively, when using a gds file as a library it is more
    efficient and convenient to call, load the GDS file only once. In that case
    use asdict=True in load_gds, which returns a dictionary with the keys referencing
    each cell by cellname, while keeping a distinction between topcells and subcells:

    {'topcells': {<cellname>: Cell}, 'subcells': {<cellname>: Cell}}

    load_gds checks for cellname clashes, because a cell(tree)
    loaded into the layout can not contain any cell name that already exists
    in the layout. There are three ways to avoid name clashes in case they occur:

    1. Set a <newcellname>: changes only the name of the top cell;
    Note <newcellname> is ignored if a <cellmap> is provided as cellmap is then
    supposed to contain the renaming.

    2. Provide a <cellmap>: maps original cell names one-by-one to new cell
    names. Cells names omitted will not be renamed.

    3. Provide a cell name <prefix>: applies to all cells in a branch,
    except those in <newcellname> or <cellmap>.

    Note, if you need to create a building block from GDS with pins and stubs,
    and you have the pin and xsection info available already (for example in a
    file) then you may prefer to use method 'load_gdsBB' instead.

    Args:
        filename (str|bytes): filename of gds to load or bytes object that
            contains the gds byte stream.
        cellname (str): cellname in filename to load
        newcellname (str): new name for <cellname>
        asdict (bool): return a dict of cells if True (default=False).
            If True, the cellname keyword will be ignored.
        select (str): set the structure of the dictionary returned by this function.
            Only used when asdict is True.
        layermap (dict): layer mapping {old_layernumber: new_layernumber}
        cellmap (dict): cellname mapping {old_cellname: new_cellname}
        scale (float): scaling factor of the cell (default=1.0).
            Scaling in this function will be applied directly to the polygon
            structures, not as a modifier attribute like scaling/magnification
            on a gds instance. For the latter use the scale keyword in 'put()'.
            Warning, use scaling with care and basically only on structures
            like logos and text. Avoid (do not) use on library building
            blocks and/or optical components as their functionality may be
            compromised by scaling.
        prefix (str): optional string to avoid name clashes (default='')
        instantiate (bool): instantiate the GDS (default=False)
        native (bool): create native Nazca cells from GDS (default=True).
            Native cells provide full access to all elements in the gds.
        bbox (bool): add bounding box if True (default=False)
        bboxbuf (float): add a buffer to the bounding box (default=0)
        connect ((float, float, float)): Only for native=False,
            move GDS origin by (x, y, a) when putting it into the layout.
        moveorg ((float, float, float)): alternative name for connect
        reuse (bool): default is False. If True, load gds assumes cells with the
            same name are the same cell and it resolves conficts by using
            only the first cell occurance. Setting reuse=True is dangerous.
        cellsnotreused (list of str): list of cellnames to ignore for
            reuse if reuse=True
        cellsreused (list of str): list of cellnames to reuse. If set all other cells will
            not be reused and reuse and cellsnotreused will be ignored.

        topcellsonly (bool): default=True

    Returns:
        Cell | {cellname: Cell}: Nazca Cell containing the loaded gds cell(tree).
            If asdict is True return a dictionary.

    Example::

        cell1 = load_gds(filename='path_to_file.gds', cellname='name')
        cell1.put()

        # using asdict = True:
        cellAll = load_gds(filename='path_to_file.gds', asdict=True)
        cellAll[<cellname>].put()

    """

    # TODO: The load_GDS method loads the gds in a stream object to check the cellnames
    #    After that, the stream is discarded. Only the filename, cell reference
    #    are stored. At mask export time the GDS is loaded again. For large GDS
    #    it may be wiser to store the stream in memory for reuse at export time.
    #print('load_gds: {} -- {}'.format(filename, cellname))
    global num
    num += 1
    cfg.gdsload = filename

# =============================================================================
# Load gds and select topcell(s) to return
# =============================================================================
    strm = gdsimp.GDSII_stream(filename, layermap=layermap, layermapmode=layermapmode)
    if scale is None:
        scale = strm.gds_db_unit / gb.gds_db_unit
        # print(f"scale: {os.path.basename(filename)}, {strm.gds_db_unit}, {gb.gds_db_unit}, {scale}")
        if scale != 1.0:
            main_logger(
                f"Imported gds '{filename}' will be scale corrected by factor {scale}.",
                "info"
            )

    if isinstance(filename, str):
        fln_or_buf = filename
    else:
        fln_or_buf = '<buffer>'

    if cellname is not None:
        topcellname = cellname
        topcells = [topcellname]
    else: # find cellname if not provided
        topcells = strm.topcell()
        topcellname = None
        if len(topcells) == 1 and topcells is not None:
            topcellname = next(iter(topcells))
        elif not asdict: # and not native:
            asdict = True
            allcells = set(list(strm.cells.keys()))
            othercells = allcells - topcells
            raise Exception("load_gds: No cellname provided to load from file '{}'."\
                " Tried to find a single topcell, but multiple topcells exist."\
                " To solve this: specify an explicit cellname to load,"\
                " or set asdict=True."\
                " to return a cell dictionary instead of a single cell."
                " Available topcells:\n{}"\
                " \nAvailable non-topcells:\n{}".\
                format(filename, topcells, othercells)) # '\n'.join(list(strm.cells))))

            #print("Error: Could not determine the topcell.")
            #topcellname = None


# =============================================================================
# Build cellmap to avoid name clashes and reload GDS:
# =============================================================================
    if cellmap is None:
        cellmap = dict()
        if newcellname is not None:
            cellmap[topcellname] = newcellname

    # Check for cellname clashes:
    for topcell in topcells:
        names = strm.cell_branch(topcell)
        for name in sorted(names):
            if name not in cellmap:
                cellmap[name] = '{}{}'.format(prefix, name)
            if (cellmap[name] in cfg.cellnames.keys()) and not native: # native solves the conflicts itself
                if cellmap[name] is not name:
                    rename = "is renamed to '{0}' which ".format(cellmap[name])
                else:
                    rename = ''
                if reuse:
                    main_logger(
                        f"use_first=True in load_gds:"
                        f" from file '{name}'",
                        "warning"
                    )
                else:
                    raise Exception(
"""
Error in load_gds: cell name overlap: '{0}'.

Cell name '{3}' in file '{1}' {4}is already in use in the design.

  Five suggestions to solve this issue by giving a <uniquename>:
    In 'load_gds':
    1. Rename the cell (topcell only): load_gds(..., newcellname=<uniquename>)
    2. Apply a cellmap (0 or more cells in a branch): load_gds(...,
       cellmap=<cellmap_dict>)
    3. Apply a prefix (to all cells in a branch): load_gds(..., prefix='new').
       Present prefix is '{2}'
    Other:
    4. Rename cell '{3}' in gds file '{1}'
    5. Rename existing cell name '{0}' in the design
""".format(cellmap[name], fln_or_buf, prefix, name, rename))

        if not native:  # or Cell will see it as existing.
            cfg.cellnames[cellmap[name]] = 'loaded_gds_cell'
    # Reload gds applying cellmaps and layer maps.
    # TODO: change in memory i.o. reload.
    strm = gdsimp.GDSII_stream(
        filename,
        layermap=layermap,
        layermapmode=layermapmode,
        cellmap=cellmap
    )

# =============================================================================
# Create nazca cell to return:
# =============================================================================
    if topcellname is None:
        maptopcellname = None
    else:
        maptopcellname = cellmap[topcellname]

    if native:
        maskcell = _gds2native(
            strm=strm,
            topcellname=maptopcellname,
            cellmap=cellmap,
            bbox=bbox,
            bboxbuf=bboxbuf,
            hull=hull,
            show_hull=show_hull,
            hull_based_bbox=hull_based_bbox,
            scale=scale,
            flat=flat,
            instantiate=instantiate,
            asdict=asdict,
            topcellsonly=topcellsonly,
            select=select,
            reuse=reuse,
            cellsreused=cellsreused,
            cellsnotreused=cellsnotreused)
    else:  # non-native
        with Cell('load_gds', celltype='mask', instantiate=instantiate) as maskcell:
            maskcell._put_gds(
                connect=connect,
                filename=fln_or_buf,
                cellname=topcellname,
                newcellname=maptopcellname,
                layermap=layermap,
                cellmap=cellmap,
                scale=scale,
                strm=strm)
            maskcell.autobbox = bbox
            maskcell.bboxbuf = bboxbuf
            maskcell.hull_based_bbox = hull

    cfg.gdsload = ""
    cfg.gdsloadstore = ""  # clear for generating new warnings in a next gds file.
    return maskcell


def load_gds_raw(filename, cellname=None, newcellname=None, layermap=None,
        cellmap=None, scale=1.0, prefix='', instantiate=False,
        bbox=False, hull=None, hull_based_bbox=None, connect=None, use_first=False,
        layermapmode='all', **kwargs):
    """Load GDS and always force native=False.

    load_gds_raw(...) is the same as load_gds(..., native=False).
    See load_gds() for a detailed method description.

    Returns:
        Cell: Nazca Cell containing the loaded gds cell(tree)
    """
    return load_gds(
        filename=filename,
        cellname=cellname,
        newcellname=newcellname,
        layermap=layermap,
        cellmap=cellmap,
        scale=scale,
        prefix=prefix,
        instantiate=instantiate,
        native=False,
        bbox=bbox,
        hull=hull,
        hull_based_bbox=hull_based_bbox,
        connect=connect,
        use_first=use_first,
        layermapmode=layermapmode)


def show_pin(pin=None, radius=10, width=1, layer='dump'):
    """Draw a ring in the layout to indicate/visualize the pin position.

    Args:
        pin (node | str): pin reference
        radius (float): radius of the ring shaped pin inidcator in um
        width (float): width of the ring-shaped pin indicator in um

    If no pin is provided (None), cp will be use to apply the "show".

    Returns:
        None
    """
    if pin is None:
        pin = cfg.cp
    Polygon(layer=layer, points=ring(radius=radius, width=width)).put(pin)
    # TODO: add direction


def show_cp(radius=10, width=1, layer='dump'):
    """Draw a ring in the layout to indicate/visualize the cp position.

    Args:
        radius (float): radius of the ring shaped pin inidcator in um
        width (float): width of the ring-shaped pin indicator in um

    Returns:
        None
    """
    if isinstance(radius, Node):
        raise Exception(f"To show pin '{radius.name}' use method show_pin() instead of show_cp().")
    show_pin(pin=None, radius=radius, width=width, layer=layer)


def log_pin(pin=None, label="", show=True, layer='dump'):
    """Generate log string where a pin is (cell + location) and lable it.

    Optionally visualize the location of the pin with the label in the layout.

    Args:
        pin (Node): pin to log. default=cp
        label (str): label to identify the log
        show (bool): show the pin in the layout
        layer (str | int | (int, int)): layer indentifier for visualization. default='dump'

    Returns:
        str: message to put in a logfile, e.g. via nd.logging.info(...)
    """
    height = 5
    if pin is None:
        pin = cfg.cp
    try:
        pinname = pin.up.name
    except:
        pinname = pin.name
    logstr = f"log_pin with label '{label}' in cell '{cfg.active_cell().cell_name}' = {pin.fxya()} for pin.name = '{pinname}', pin.id = {pin.id}"
    if show:
        show_pin(pin=pin, layer=layer)
        cp = cfg.cp
        font.text(text=label, height=height, align='cc', layer=layer).put(pin)
        cfg.cp = cp
    try:
        logger.info(logstr)
    except:
        return logstr


# See also show_pin()
wherecnt = count()
def whereami(text='here', size=100, pin=None):
    """Show current pointer position as arrow in the layout.

    Args:
        text (str): annotation text
        size (float): size of the annotation

    Returns:
        None
    """
    layer = 500
    cp_store = cfg.cp
    if pin is None:
        pin = cfg.cp
    points = [
        (0, 0), (-0.5, 0.5), (-0.4, 0.25), (-1, 0.3), (-1, -0.3),
        (-0.4, -0.25), (-0.5, -0.5)]
    points = [(x*size,y*size) for x, y in points]
    with Cell('I am'.format(wherecnt)) as crisis:
        Polygon(layer=layer, points=points).put(0)
        text(text+' ', height=size/4, align='rc', layer=layer).put(0)

    crisis.put(pin)
    cfg.cp = cp_store


tab_elm = '. '
class Netlist:
    """Class for building the netlist."""

    def __init__(self):
        """ Construct a Netlist object. Singleton.
        """
        self.nodes_visited = set()
        self.cnodes_visited = set()


    def solvecell_core(self, node, _cnode0=None):
        """Calculate the geometrical position of nodes in a cell via the netlist.

        Include all nodes in the cell (level 0) and all nodes in instantiated
        cells (level 1). Hence, a node in scope adheres to one of the following
        conditions:

        - the node resides in same cell as the starting node (parent <cnode0>).
        - the node resides in an instance of the cell of the starting node.

        A starting <node> has to be provided from where to solve.
        A good starting point in each cell is the cnode with its pointer
        position set at (x, y, a) = (0, 0, 0).

        Args:
            node (Node): node from where to start solving the cell
            _cnode0 (Node): Internal parameter containing the cnode of the cell
                to solve. It is automatically obtained from <node> in theh
                first pass.

        Returns:
            None
        """

        if _cnode0 is None:
            _cnode0 = node.cnode
        self.nodes_visited.add(node)

        neighbours = node.geo_nb_iter()
        for nextnode, mat in neighbours:
            if nextnode not in self.nodes_visited:
                #print(node)
                #node in main cell or in instance
                if (nextnode.cnode is _cnode0):
                    nextnode.pointer.set_mat(node.pointer.trans(mat))
                    self.solvecell_core(nextnode, _cnode0)
                elif (nextnode.cnode.parent_cnode is _cnode0):
                    nextnode.pointer.set_mat(node.pointer.trans(mat))
                    self.solvecell_core(nextnode, _cnode0)
                else:
                    logger.warning(
                        f"Skipping node {nextnode} in cell '{nextnode.cnode.cell.cell_name}': "
                        f"The node is not in scope of cell '{ _cnode0.cell.cell_name}'.",
                        "warning"
                    )
            # TODO: else: check geo consistency, e.g. floating point delta's.


    def solvecell(self, cell):
        """Solve the cell with. cnode position will be at (x, y, a) = (0, 0, 0).

        Args:
            cell (Cell): cell to solve

        Returns:
            None
        """
        if cfg.solve_direct:
            return None
        cell.cnode.pointer.goto(0, 0, 0)
        self.solvecell_core(cell.cnode)
        return None


# =============================================================================
# Cell-element iterators that return Node and a cell-element object
# =============================================================================
    def polygon_iter(self, cnode, level, infolevel=0):
        """Generator to iterate over all polygons in a cell (cnode).

        Yields:
            (Pointer, Polygom): Next Polygon object and its position
        """
        for node, pgon in cnode.cell.polygons:
            if infolevel > 1:
                print('{}polygon: ML={}, xya={}, bbox={}'.format(
                    '  '*level,
                    pgon.layer,
                    node.pointer.xya(),
                    pgon.bbox
                    )
                )
            yield node.pointer.copy(), pgon


    def polyline_iter(self, cnode, level, infolevel=0):
        """Generator to iterate over all polylines (paths) in a cell (cnode).

        Yields:
            (Pointer, Polyline): Next Polyline object and its position
        """
        for node, pline in cnode.cell.polylines:
            if infolevel > 1:
                print('{}polyline: ML={}, xya={}, bbox={}'.format(
                    '  '*level,
                    pline.layer,
                    node.pointer.xya(),
                    pline.bbox
                    )
                )
            yield node.pointer.copy(), pline


    def annotation_iter(self, cnode, level, infolevel=0):
        """Generator to iterate over all annotations in a cell (cnode).

        Yields:
            (Pointer, annotation): Next Annotation and its position
        """
        for node, anno in cnode.cell.annotations:
            if infolevel > 1:
                print('{}annotation: ML={}, xya={}'.format(
                    '  '*level,
                    anno.layer,
                    node.pointer.xya()
                    )
                )
            yield node.pointer.copy(), anno


    def instance_iter(self, cnode):
        """Generator to iterate over all instances in cell (cnode) one level deep.

        Yields:
            Node: cnode of instance
        """
        for nn in cnode.nb_cnode:
            if nn[2] == 1:
                yield nn[0]


    def gdsfile_iter(self, cnode, level, infolevel=0):
        """Generator to iterate over all GDS files instantiated in a cell (cnode).

        Yields:
            (Pointer, gdsinfo): Next gds file and info to generate it.
                gdsinfo = (filename, cell, newcell, layermap, cellmap, scale)
        """
        gdsfiles = cnode.cell.gdsfiles
        for node, gdsinfo in gdsfiles:
            if infolevel > 1:
                print('{}gds-cell: cell={}, xya={}'.format(
                    '  '*level,
                    gdsinfo.newcell,
                    cnode.pointer.xya()
                    )
                )
            yield node.pointer.copy(), gdsinfo


# =============================================================================
# Cell-element iterators WITH translation and flipping performed
# =============================================================================
    def polygon_transflip_iter(
        self,
        cnode,
        trans,
        flip,
        scale,
        level=0,
        apply=True,
        hierarchy='self',
        infolevel=0,
        hull=False,
    ):
        """Generator to iterate over all polygons in a cell (cnode).

        Args:
            cnode (Node): cnode
            trans (tuple): translation state w.r.t. to (new) parent
            flip (bool): flip state w.r.t. to (new) parent
            scale (float): scale state w.r.t. to (new) parent
            apply (bool): apply translation and flipping (default=True)
            infolevel (int): debug info level to stdout

        Yields:
            (Polygon, list((float, float)), tuple(4*float)) |
            (Polygon, tuple(3*float), tuple(4*float)):

            if apply is True (default)::

                Polygon object, and the polygon points & the bbox w.r.t.
                the (possibly new) parent cell.

            if apply is False::

                Polygon object, its position (x, y, a), 'a' in radians,
                the flip state of the polygon
        """
        # TODO: how to deal with  hierarchy == 'self' and apply
        #if hierarchy == 'self':
        #    apply = False
        for node, pgon in cnode.cell.polygons:
            if infolevel > 1:
                print("{0}polygon: ML='{ld}', "\
                    'xya=({p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}), '\
                    'bbox_wxh=({bx:.3f}, {by:.3f}), '\
                    'N={1}'.format(
                        tab_elm*(level+2),
                        len(pgon.points),
                        ld=pgon.layer,
                        p=node.pointer.xya(),
                        bx=pgon.bbox[2]-pgon.bbox[0],
                        by=pgon.bbox[3]-pgon.bbox[1]
                    )
                )
            org = node.pointer.copy()
            nflip, nflop = 1, 1
            if node.flip:
                nflip = -1
            if node.flop:
                nflop = -1
            s = 1
            if flip:
                s = -1
                org = org.copy()
                org.flip()
            [x, y, a] = org.multiply_ptr(trans).xya()
            a = s*np.radians(a)
            scale_tot = scale * node.scale
            if hull and hasattr(pgon, 'hull'):
                points = pgon.hull
            else:
                points = pgon.points
            if apply:
                if (
                    [x, y, a] != [0, 0, 0]
                    or nflip == -1
                    or nflop == -1
                    or scale_tot != 1.0
                    or s == -1
                ):
                    xy = [(x+scale_tot*(cos(a)*nflip*u-sin(a)*nflop*v),
                           y+scale_tot*s*(sin(a)*nflip*u+cos(a)*nflop*v))
                        for u, v in points]
                else:
                    xy = points
                xyT = list(zip(*xy))
                bbox = (min(xyT[0]), min(xyT[1]), max(xyT[0]), max(xyT[1]))
                yield pgon, xy, bbox
            else:
                yield pgon, pgon.points, pgon.bbox
            # TODO: can this realy be removed?:
            #else:
            #    print(f"other-xya: {x, y, a}")
            #    yield pgon, [x, y, a], s


    def polyline_transflip_iter(
        self,
        cnode,
        trans,
        flip,
        scale,
        level=0,
        apply=True,
        hierarchy='self',
        infolevel=0,
        hull=False,
        polygon=False,
    ):
        """Generator to iterate over all polylines in a cell (cnode).

        Args:
            cnode (Node): cnode
            trans (tuple): translation state w.r.t. to (new) parent
            flip (bool): flip state w.r.t. to (new) parent
            scale (float): scale state w.r.t. to (new) parent
            apply (bool): apply translation and flipping (default=True)
            hull (bool): if True work with hull points
            polygon (bool): if True work with the polygon rather then polyline points
            infolevel (int): debug info level to stdout

        Yields:
            (Polyline, list((float, float)), tuple(4*float)) |
            (Polyline, tuple(3*float), tuple(4*float)):

            if apply is True (default):

                Polyline object, and the polyline points & the bbox w.r.t.
                the (possibly new) parent cell.

            if apply is False:

                Polyline object, its position (x, y, a), 'a' in radians,
                the flip state of the polyline
        """
        #if hierarchy == 'self':
        #    apply = False
        if cfg.export_polyline_as_polygon:
            polygon = True
            hull = False
        for node, pline in cnode.cell.polylines:
            if infolevel > 1:
                print("{0}polyline: ML='{ld}', "\
                    'xya=({p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}), '\
                    'bbox_wxh=({bx:.3f}, {by:.3f}), '\
                    'N={1}'.format(
                    tab_elm*(level+2),
                    len(pline.points),
                    ld=pline.layer,
                    p=node.pointer.xya(),
                    bx=pline.bbox[2]-pline.bbox[0],
                    by=pline.bbox[3]-pline.bbox[1]
                    )
                )
            org = node.pointer.copy()
            nflip, nflop = 1, 1
            if node.flip:
                nflip = -1
            if node.flop:
                nflop = -1
            s = 1
            if flip:
                s = -1
                org = org.copy()
                org.flip()
            [x, y, a] = org.multiply_ptr(trans).xya()
            a = s*np.radians(a)
            scale_tot = scale * node.scale
            if hull and hasattr(pline, 'hull'):
                points = pline.hull
            elif polygon:
                points = pline.polygon
            else:
                points = pline.points
            if apply:
                if (
                    [x, y, a] != [0, 0, 0]
                    or nflip == -1
                    or nflop == -1
                    or scale_tot != 1.0
                    or s == -1
                ):
                    xy = [(x+scale_tot*(cos(a)*nflip*u-sin(a)*nflop*v),
                           y+scale_tot*s*(sin(a)*nflip*u+cos(a)*nflop*v))
                        for u, v in points]
                else:
                    xy = points
                xyT = list(zip(*xy))
                bbox = (min(xyT[0]), min(xyT[1]), max(xyT[0]), max(xyT[1]))
                yield pline, xy, bbox
            else:
                yield pline, pline.points, pline.bbox
            #else:
            #    yield pline, [x, y, a], s


    def annotation_transflip_iter(
        self,
        cnode,
        trans,
        flip,
        scale,
        level=0,
        apply=True,
        hierarchy='self',
        infolevel=0,
    ):
        """Generator to iterate over all annotations in a cell (cnode).

        Args:
            cnode (Node): cnode to iterate
            trans (tuple): translation state of instance
            flip (bool): flip state of instance
            apply (bool): default=True: Apply translation and flipping
            infolevel (int): amount of debug info to display in stdout

        Includes translation and flipping.

        Yields:
            (Annotation, (float, float), int):
                Next Annotation and its position as (x, y),
                and a flip multiplyer -1 (flip) or 1
        """
        for node, anno in cnode.cell.annotations:
            if infolevel > 1:
                print("{0}annotation: ML='{ld}', "\
                    'xy=({p[0]:.3f}, {p[1]:.3f}), '\
                    "text='{t}'".format(
                    tab_elm*(level+2),
                    ld=anno.layer,
                    p=node.pointer.xya(),
                    t=anno.text[:5]
                    )
                )
            org = node.pointer.copy()
            if flip:
                org = org.copy()
                org.flip()
            xy = org.multiply_ptr(trans).xy()
            yield anno, xy


    def instance_transflip_iter(
        self,
        cnode,
        trans,
        flip,
        scale,
        level=0,
        apply=True,
        hierarchy='self',
        infolevel=0,
    ):
        """Generator to iterate over all instances one level deep.

        Includes translation and flipping.

        Args:
            cnode (Node): cnode to iterate.
            trans (tuple): translation state of instance.
            flip (bool): flip state of instance.
            scale (float): scaling factor.
            apply (bool): default=True: Apply translation and flipping.
            infolevel (int): amount of debug info to display in stdout.
            hierarchy (str): setting on how to process the hierarchy.

        Yields:
            cnode, (float, float, float), bool: cnode of instance (inode),
                position (x, y, a), 'a' in degrees, flip state
        """
        if hierarchy == 'flat':
            return None
        for nn in cnode.nb_cnode:
            if nn[2] == 1:  # cell link is down stream
                #if not nn[0].cell.instantiate and hierarchy != 'full':
                #    return None
                org = nn[0].pointer.copy()
                if flip:
                    org = org.copy()
                    org.flip()
                [x, y, a] = org.multiply_ptr(trans).xya()
                #print(f"|instance: {nn[0].cell.cell_name}, instantiate:{nn[0].cell.instantiate}")
                yield nn[0], [x, y, a], flip


    def gdsfile_transflip_iter(
        self,
        cnode,
        trans,
        flip,
        scale,
        apply,
        level=0,
        hierarchy='self',
        infolevel=0,
    ):
        """Generator to iterate over all GDS files instantiated in a cell (cnode).

        Args:
            cnode (Node): cnode to iterate
            trans (tuple): translation state of instance
            flip (bool): flip state of instance
            apply (bool): ignored: always True
            infolevel (int): amount of debug info to display in stdout

        Yields:
             (Polyline, (float, float, float), bool):
                Next gdsfile and its position (x, y, a), 'a' in degrees,
                and the flip state
        """
        gdsfiles = cnode.cell.gdsfiles
        for node, gdsinfo in gdsfiles:
            if infolevel > 1:
                print('{}gds-cell: cell={}, xya={}'.format(
                    tab_elm*(level+2),
                    gdsinfo.newcell,
                    cnode.pointer.xya()
                    )
                )
            org = node.pointer.copy()
            [x, y, a] = org.multiply_ptr(trans).xya()
            yield gdsinfo, [x, y, a], flip


# =============================================================================
# Cell iteration
# =============================================================================
    def celltree_iter_base(
        self,
        icnode,
        hierarchy='full', # ['flat', 'self', 'full', 'apply']
        position=None,
        level=0,
        infolevel=0
    ):
        """Generator to iterate over cell of <cnode>, top-down, go deep first.

        The tree descents from the provided cnode into instantiated cnode(s)
        until the no deeper instances are available.
        For decending into an instance, the inode:
        - look up the cell object the instance represents,
        - take that cell's cnode and use it to seed the next level.

        Args:
            icnode (Node): cnode of cell or inode of instance to iterate into
            level (int): depth in cell tree, (default=0)
            position (Pointer): (default=None)
            flat (bool): flatten array instances if True (default=False)
            infolevel (int): amount of runtime info to stdout (default=0)

        Yields:
            (Node, int, Pointer, bool, scale): (cnode, level, position, flip).
                Yields next cell, its level in the hierarchy and position
        """
        # Get the cnode of the cell cnode or inode
        #   Hence, always start the tree with a cell.
        cnode = icnode.cell.cnode
        if position is None:
            position = Pointer(0, 0, 0)
        if infolevel > 3:
            print('{}: cnode={}, name={}, level={}'.format(
                '  '*infolevel, cnode, cnode.cell.cell_name, level)
            )

        # stack of cnode-tree:
        path = []
        nodes = [(cnode, level, position, False, 1.0, cnode)]
        while nodes:
            level_old = level
            cnode, level, position, flip, scale, next_node = nodes[-1]
            if level > 0 and level <= level_old:
                for _ in range(level_old - level + 1):
                    del path[-1]
            path.append(next_node)
            #P = '/'.join([a.cnode.cell.cell_name for a in path])
            #print(level, len(nodes), P)
            yield cnode, level, position, flip, scale, path
            self.cnodes_visited.add(cnode)
            del nodes[-1]
            next_nodes = list(cnode.instance_iter())
            for next_node in next_nodes:  # append nodes list with all next_nodes
                cnode = next_node.cell.cnode
                flip = next_node.flip
                scale = next_node.scale
                array = next_node.array
                if (array is not None
                    and (not next_node.cell.instantiate or hierarchy=='flat')
                ):
                    nx, (dx1, dy1), ny, (dx2, dy2) = array
                    x0, y0, a0 = next_node.pointer.xya()
                    for xi in range(nx):
                        for yi in range(ny):
                            orgnew = Pointer(x0+xi*dx1+yi*dx2, y0+xi*dy1+yi*dy2, a0)
                            nodes.append((cnode, level+1, orgnew, flip, scale, next_node))
                else:
                    nodes.append((cnode, level+1, next_node.pointer, flip, scale, next_node))


    def celltree_iter(
            self,
            cell,
            position=None,
            hierarchy='self',
            cells_visited=None,
            infolevel=0,
            cellmap=None,
            topdown=False,
            revisit=False
        ):
        """Generator to iterate top-down or bottom-up (default) over the celltree in <cell>.

        Translation, flipping and flattening are applied.
        Go deep first. If a instantiation is found of a cell that already has
        been iterated over (cells_visited), then this cell will not be
        scanned again but reused.

        Iteration over this Generator with infolevel = 1 will print the Nazca
        hierarchy, cell-levels to stdout.

        Args:
            cell (Cell): topcell to iterate into
            position ((float, float float)): xya starting postition of the topcell.
                default=None, which translates into (0, 0, 0)
            hierarchy (str): 'flat', 'self' (default) or full'.
                Determines the handling of the instatiate cell attribute:
                None: follow the cells instaniate setting to flatten cells or not.
                'flat': collapses the hierarcy into a
                single level corresponding to setting instantiate=False in all cells.
                'self' (default): reproduces all cells explicitly while maintaining their instantiate setting.
                'full': instantiates all cells.
                This setting is particular relevant when adding instances to a cell
                as based in this setting. Coordinates returned will be relative to
                the first instantiated parent as set by this keyword.
            cells_visited (set): list of cells_visited. The iterator will skip
                cells in this set and not iterate into them.
            infolevel (int): amount of info printed when iterating
            cellmap (dict): Optional dict to rename cells {oldname: newname}.
                default=None
            topdown (bool): indicate topdown (True) or bottom up yielding iteration (False, default)
            revisit (bool): go into each instance of each cell (default=False)

        Yields:
            named_tuple: Info to (re)build cells in a netlist branch.
        """
        hierarchy_allowed = ['flat', 'apply', 'self', 'full']
        if hierarchy not in hierarchy_allowed:
            raise ValueError("Unkown hierarchy setting '{}', not in {hierarchy_allowed}")

        stack_export_level = {} #  store parent levels that are close after non-instantiated children.
        def up(depth, stop):
            """Move upwards in the celltree.

            Yields:
                named_tuple: Info to (re)build cells in a netlist branch.
            """
            while depth >= stop:
                # get a cellopen value:
                if not topdown:  # bottom-up
                    if export_levels[-1] in stack_instant.keys():
                        # yield cell opening for export level:
                        cellparams = stack_instant.pop(export_levels[-1])
                        stack_export_level[export_levels[-1]] = cellparams
                        yield cellparams

                    # yield a non-instantiated level:
                    if depth not in export_levels:
                        cellparams = stack_noninstant.pop()
                        yield cellparams

                export_flag = depth in export_levels
                if export_flag:
                    cellparams = stack_export_level[depth]

                if translist_loc:  # not []: go level up in location lists
                    translist_loc.pop()
                    fliplist_loc.pop()
                    translist_glob.pop()
                    fliplist_glob.pop()

                if infolevel > 0:
                     print("{}clse cell: '{}', inst={}, export={}".format( # clse i.o. close for flush stdout
                         tab_close*(cellparams.level+1),
                         cellparams.cell.cell_name,
                         cellparams.instantiate,
                         export_flag,
                         )
                     )

                if depth in export_levels:
                    export_levels.pop()
                    CLOSE = CLOSECELL  # if an export level yield cell closure:
                else:  # no export_level
                    scalelist.pop()
                    CLOSE = not CLOSECELL
                yield cellinfo_close(
                    cell_start=None,
                    cell_instantiate=None,
                    cell_close=CLOSE,
                    cell_end=True,
                    cell=cellparams.cell,
                    level=depth,
                    cell_create=False,
                    cell_open=False,
                    iters=iters_close,
                    new_cell_name='',
                    hierarchy=hierarchy,
                    branch=cellparams.branch,
                )
                depth -= 1
            return depth

        cellinfo_open = namedtuple('cellinfo_open',
            ['cell_start',       # Old name for cell_instantiate
             'cell_instantiate', # True for a new nazca cell, instantiated or not
             'cell_create',      # True for a new instantiated cell
             'cell_open',        # = cell_create
             'cell',             # cell object
             'new_cell_name',    # name of the cell
             'level',            # depth od the cell in the hierarchy (top=0)
             'parent_level',     # parent level, taken into account flattend levels
             'transflip_loc',    # position and flip state w.r.t. the parent cell
             'transflip_glob',   # position and flip state w.r.t. the to cell
             'scale',            # scaling factor of the parent cell after export
             'iters',            # list of iterators of the elements in the cell
             'cell_close',       # for compatibility with cellinfo_close tuple only
             'cell_end',         # for compatibility with cellinfo_close tuple only
             'instantiate',      # bool based on hierarchy setting: 'flat', 'self' or 'full'
             'hierarchy',        # how to deal with instantiation
             'branch',
            ]
        )
        cellinfo_close = namedtuple('cellinfo_close',
            ['cell_close',       # True for a cell closure of an instantiated cell
             'cell_end',         # cell ending, instantiated or not
             'level',            # level of the cell in the hierarchy (top=0)
             'cell',             # cell object
             'iters',            # for compatibility with cellinfo_open tuple only
             'cell_start',       # for compatibility with cellinfo_open tuple only
             'cell_instantiate', # for compatibility with cellinfo_open tuple only
             'cell_create',      # for compatibility with cellinfo_open tuple only
             'cell_open',        # for compatibility with cellinfo_open tuple only
             'new_cell_name',    # for compatibility with cellinfo_open tuple only
             'hierarchy',        # how to deal with instantiation
             'branch',
            ]
        )

        tab_open  = '> '
        tab_close = '< '
        tab_reuse = '* '
        tab_flat  = '_ '
        if infolevel > 0:
            print("Start Netlist (infolevel={})\n"
                "{} -> open an instantiated cell\n"
                "{} -> close a cell\n"
                "{} -> reused instantiated cell\n"
                "{} -> non-instantiated cell\n"
                "{} -> elements (polygons, etc.)\n"
                "-------------------------------".format(
                infolevel, tab_open, tab_close, tab_reuse, tab_flat, tab_elm))

        topdown = topdown # return cells in a topdown order if True
        stack_instant = {} # TODO: clearer description of relation with 'export_levels'.
        stack_noninstant = []
        if cellmap is None:
            cellmap = {}
        NEWCELL = True
        CLOSECELL= True
        iters_close = {
            'polyline': [],
            'polygon': [],
            'annotation': [],
            'instance': [],
            'gdsfile': []}
        export_levels = []  # stack for exported netlist levels.
        translist_loc = []  # stack for local translation
        fliplist_loc = []   # stack for local flip state
        translist_glob = [] # stack for global translation
        fliplist_glob = []  # stack for global flip state
        scalelist = [1]     # stack for scaling of flattened instances
        level_last = 0      # top level is 0
        level_reuse = 100   # level

        if cells_visited is None:
            self.cells_visited = set()
        else:
            self.cells_visited = cells_visited
        # TODO: do not reuse a top cell in a flattend design
        #       it may become missing in a gds libary as a separate top cell.

        cell_iter = self.celltree_iter_base(
            icnode=cell.cnode,
            level=0,
            position=position,
            hierarchy=hierarchy,
            infolevel=infolevel,
        )

        for cnode, level, org, flip, scale, branch in cell_iter:
            if level > level_reuse:
                # skip cells that are part of a reused cell.
                continue
            else:
                level_reuse = 100

            # close cells:
            depth, stop = level_last, max(1, level)
            depth = yield from up(depth, stop)

            cellname = cnode.cell.cell_name

            if infolevel > 0:
                info = ("open cell: '{}' @ ({:.3f}, {:.3f}, {:.3f}), flip={}, inst={}".\
                    format(cellname, org.x, org.y, org.a,
                        flip, cnode.cell.instantiate))

            if hierarchy == 'self':
                inst = cnode.cell.instantiate
                makecell = True
            elif hierarchy == 'apply':
                inst = cnode.cell.instantiate
                makecell = cnode.cell.instantiate
            elif hierarchy == 'flat':
                if level == 0:
                    inst = True
                    makecell = True
                else:
                    inst = False
                    makecell = False
            elif hierarchy == 'full' :
                inst = True
                makecell = True
            else:
                raise ValueError(f"hierarchy value '{hierarchy}' unknown")

            visited = cellname in self.cells_visited

            if not makecell or not visited or revisit: # create new cell or recreate flat cell:
                if infolevel > 0:
                    text = 'process ' if infolevel > 2 else ''
                    if inst:
                        print('{}{}{}'.format(tab_open*(level+1), text, info))
                    elif infolevel > 0:
                        print('{}{}{}'.format(tab_flat*(level+1), text, info))

                self.cells_visited.add(cellname)

                # global transflipping:
                if level == 0:
                    translist_glob.append(Pointer(0, 0, 0))
                    fliplist_glob.append(False)
                else:
                    translate = translist_glob[-1].copy()
                    orgtrans = org.copy()
                    if fliplist_glob[-1]:
                        orgtrans.flip()
                    fliplast = fliplist_glob[-1] ^ flip
                    fliplist_glob.append(fliplast)
                    translate.move_ptr(orgtrans)
                    translist_glob.append(translate)

                # local transflipping:
                if makecell or level == 0: # create new export level
                    createcell = True
                    export_levels.append(level)
                    fliplast = flip
                    translist_loc.append(Pointer(0, 0, 0))
                    fliplist_loc.append(False)
                else: # flattened cell -> content translates into a parent cell.
                    createcell = False
                    scalelist.append(scalelist[-1]*scale)
                    translate = translist_loc[-1].copy()
                    orgtrans = org.copy()
                    if fliplist_loc[-1]:
                        orgtrans.flip()
                    fliplast = fliplist_loc[-1] ^ flip
                    fliplist_loc.append(fliplast)
                    translate.move_ptr(orgtrans)
                    translist_loc.append(translate)

                transflip_loc = translist_loc[-1], fliplist_loc[-1]
                transflip_glob = translist_glob[-1], fliplist_glob[-1]
                apply_transflip = True # TODO: connect apply to hierarchy?
                internal_params = (
                    cnode,
                    translist_loc[-1],
                    fliplist_loc[-1],
                    scalelist[-1],
                    level,
                    apply_transflip,
                    hierarchy,
                    infolevel
                )
                element_iters = {
                    'polyline':   self.polyline_transflip_iter(*internal_params),
                    'polygon':    self.polygon_transflip_iter(*internal_params),
                    'annotation': self.annotation_transflip_iter(*internal_params),
                    'instance':   self.instance_transflip_iter(*internal_params),
                    'gdsfile':    self.gdsfile_transflip_iter(*internal_params),
                }
                parent_level = export_levels[-1]
                new_cell_name = cellmap.get(cellname, None)

                iterparams = cellinfo_open(
                    cell_start=NEWCELL,
                    cell_instantiate=NEWCELL,
                    cell=cnode.cell,
                    cell_create=createcell,
                    cell_open=createcell,
                    cell_close=False,
                    cell_end=False,
                    new_cell_name=new_cell_name,
                    level=level,
                    parent_level=parent_level,
                    transflip_loc=transflip_loc,
                    transflip_glob=transflip_glob,
                    scale=scalelist[-1],
                    iters=element_iters,
                    instantiate=inst,
                    hierarchy=hierarchy,
                    branch=branch.copy(),
                )
                if topdown:
                    yield iterparams
                else:
                    if createcell:
                        stack_instant[level] = iterparams
                    else:
                        stack_noninstant.append(iterparams)
                level_last = level

            elif makecell: # reuse previously processed cell:
                if infolevel > 0:
                    text = 'reuse ' if infolevel > 1 else ''
                    print('{}{}{}'.format(tab_reuse*(level+1), text, info))
                level_reuse = level
                level_last = level-1
        # end cell_iter

        # close cells:
        if topdown:
            depth = level
        else:
            depth = level_last
        yield from up(depth, 0)
        if infolevel > 0:
            print("End celltree iteration")


def cell_iter(
        topcell,
        cellmap=None,
        flat=False,
        hierarchy='self',
        infolevel=0,
        topdown=False,
        revisit=False,
):
    """Get a cell iterator fo rebuilding cells.

    Args:
        topcell (Cell): cell to iterate into.
        cellmap (dict): Optional dict to rename cells {oldname: newname}.
            default=None.
        flat (bool): for backward compatibility, same as hierarchy='flat' if True
            and hierarchy not provided, but ignored if hierarchy is provided explicitly.
        hierarchy (str): 'flat', 'self' (default), 'full', or 'apply'.
        topdown (bool): indicate topdown (True) or bottom up iteration (False, default)
        revisit (bool): go into each instance of each cell (default=False)

    Returns:
        Iterator: celltree_iter()
    """
    if flat:
        hierarchy = 'flat'

    return Netlist().celltree_iter(
        topcell,
        hierarchy=hierarchy,
        cellmap=cellmap,
        topdown=topdown,
        revisit=revisit,
        infolevel=infolevel,
    )


#little bit ugly, but needed for now to make the Cell.rebuild method
import nazca.layout as layout
#-----------------------------------------------------------------------------
if __name__ == '__main__':
    pass
