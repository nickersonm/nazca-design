#!/usr/bin/env python3
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
# -*- coding: utf-8 -*-
# 2017 (c)  Ronald Broeke


"""
The 'cp' module provides short syntax methods to operate on current pointer cp.
The cp holds position of the standard output of the last element
in a mask layout. This module should be rarely needed be designers as
translations of positions can typically be best obtained by applying
the translation methods of the Pin object.
"""

from . import cfg
from .netlist import show_cp
from .logging import logger

stack = []

def show(radius=10, width=1):
    """Show location of cp by drawing a circle around it.

    Args:
        radius (float): radius of the circle
        width (float): width of the circle line

    Returns
        Node: cp
    """
    show_cp(radius, width)
    return cfg.cp

def here():
    """Get cp.

    Returns:
        Node: cp
    """
    return cfg.cp

def push():
    """Push cp on the FILO stack.

    Returns:
        Node: pushed cp
    """
    stack.append(cfg.cp)
    return cfg.cp

def pop():
    """Pop cp from the FILO stack.

    Returns:
        Node: popped cp
    """
    cfg.cp = stack.pop()
    return cfg.cp

def goto(x=0, y=0, a=0):
    """Goto absolute position with cp.

    Args:
        x (float): x-position
        y (float): y-position
        a (float): a-position

    Returns:
        Node: cp at goto position
    """
    cfg.cp.goto(x, y, a)
    return cfg.cp

def goto_pin(p):
    """Goto absolute position with cp.

    Args:
        p (Node): (x, y, a) position

    Returns:
        Node: cp at goto position
    """
    cfg.cp = p
    return cfg.cp

def move(x=0, y=0, a=0):
    """Move relative with respect cp.

    Args:
        x (float): x-position
        y (float): y-position
        a (float): a-position

    Returns:
        Node: cp at goto position
    """
    cfg.cp = cfg.cp.move(x, y, a)
    return cfg.cp

def shift(x=0, y=0):
    """Move relative with respect cp by x, y.

    Args:
        x (float): delta x
        y (float): delta y

    Returns:
        Node: cp at shifted position
    """

    cfg.cp = cfg.cp.move(x, y, 0)
    return cfg.cp

#def move2(x=0, y=0, a=0):
#    cfg.cp = cfg.cp.move2(x, y, a)
#    return cfg.cp

def rotate(a=0):
    """Rotate the cp.

    Args:
        a (float): a-position

    Returns:
        Node: cp at rotated position
    """
    cfg.cp = cfg.cp.rotate(a)
    return cfg.cp
rot = rotate

def skip(x=0):
    """Move in the direction of cp.

    Args:
        x (float): delta x

    Returns:
        Node: cp at skip position
    """
    cfg.cp = cfg.cp.skip(x)
    return cfg.cp

def offset(y=0):
    """Move perpendicular with respect to cp.

    Args:
        y (float): delta y

    Returns:
        Node: cp at offset position
    """
    cfg.cp = cfg.cp.offset(y)
    return cfg.cp
os = offset

#TODO: flip and flop
#TODO: this give the coordinates local to the defined cell.
#  Could be better to make is with respect to the active cell scope.
def get_xya():
    """Get cp coordinates x, y, a.

    Returns:
        (float, float, float): x, y, a
    """
    return cfg.cp.pointer.get_xya()
xya = get_xya

def get_xy():
    """Get cp coordinates x, y.

    Returns:
        (float, float): x, y
    """
    return cfg.cp.pointer.get_xy()
xy = get_xy

def get_a():
    """Get cp coordinates a.

    Returns:
        float: a
    """
    return cfg.cp.pointer.a

def x():
    """Get cp coordinates x.

    Returns:
        float: x
    """
    return cfg.cp.pointer.x

def y():
    """Get cp coordinates y.

    Returns:
        float: y
    """
    return cfg.cp.pointer.y

a = get_a


