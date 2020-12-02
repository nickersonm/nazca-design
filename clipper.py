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
#
# Grow and merge polygons.
# Uses the pyclipper module
#
# (c) 2017 Xaveer Leijtens, Ronald Broeke 2018 (c)
#
"""
Module with grow and merge polygon functions.
"""
import nazca.cfg as cfg
from nazca.logging import logger

try:
    import pyclipper
    cfg.PYCLIPPER = True
except Exception as e:
    cfg.PYCLIPPER = False
    #print('Warning: Could not import pyclipper:', e)
if cfg.PYCLIPPER:
    st = pyclipper.scale_to_clipper
    sf = pyclipper.scale_from_clipper


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
            logger.error("Could not load module 'pyclipper'."\
                " Skipping requested pyclipper functionality.")
            clipper_check = True
        return False


def merge_polygons(paths, accuracy=0.001):
    """Merge polygons and return a list with the resulting polygon(s).

    Args:
        paths (list): list of polygons. A polygon is list of coordinates (x, y).
        accuracy (float): accuracy [µm] of the location of the intermediate points.
            There is no good reason to deviate from the default (default = 0.001 µm).

    Returns:
        list: list of polygons  [ [(x1, y1), (x2, y2), ...], [...], ... ]
    """
    if not _has_pyclipper():
        return paths
    sc = 1 / accuracy
    pc = pyclipper.Pyclipper()
    pc.AddPaths(st(paths, sc), pyclipper.PT_SUBJECT, True)
    mp = pc.Execute(pyclipper.CT_UNION, pyclipper.PFT_NONZERO,
        pyclipper.PFT_NONZERO)
    return sf(mp, sc)


def grow_polygons(paths, grow=5, accuracy=0.1, jointype='round'):
    """Grow polygons and return the grown structures.

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
    if grow == 0: # 0-grow would affect polygons because of the grid.
        return paths
    sc = 1 / accuracy
    # Offset is the term for growing.
    pco = pyclipper.PyclipperOffset()
    # Join type
    jt = {'round': pyclipper.JT_ROUND, 'square': pyclipper.JT_SQUARE,
        'miter': pyclipper.JT_MITER}
    if jointype not in jt:
        print("jointype '{}' unknown.".format(jointype))
        print("jointype should be one of 'round', 'square', 'miter'.")
        print("Using default ('round')")
        jointype = 'round'
    pco.AddPaths(st(paths, sc), jt[jointype], pyclipper.ET_CLOSEDPOLYGON)
    return sf(pco.Execute(grow*sc), sc)


def diff_polygons(paths_A, paths_B, accuracy=0.001):
    """Subtract polygons and return a list with the resulting polygon(s).

    Args:
        paths_A (list): list of polygons that each have a list of (x,y) coordinates.
        paths_B (list): list of polygons that each have a list of (x,y) coordinates.
        accuracy (float): accuracy [µm] of the location of the intermediate points.
        There is no good reason to deviate from the default (default = 0.001 µm).

    Returns:
        list: List of polygons that result from subtraction: paths_A - paths_B.
    """
    if not _has_pyclipper():
        return paths_A + paths_B
    sc = 1 / accuracy
    pc = pyclipper.Pyclipper()
    pc.AddPaths(st(paths_A, sc), pyclipper.PT_SUBJECT, True)
    pc.AddPaths(st(paths_B, sc), pyclipper.PT_CLIP, True)
    sp = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO,
            pyclipper.PFT_NONZERO)
    return sf(sp, sc)


# TODO: add clipping

def clip_polygons(paths, clip, accuracy=0.001):
    if not _has_pyclipper():
        return paths
    sc = 1 / accuracy
    pc = pyclipper.Pyclipper()
    pc.AddPath(st(clip, sc), pyclipper.PT_CLIP, True)
    pc.AddPaths(st(paths, sc), pyclipper.PT_SUBJECT, True)
    sp = pc.Execute(
        pyclipper.CT_INTERSECTION,
        pyclipper.PFT_EVENODD,
        pyclipper.PFT_EVENODD)
    return  sf(sp, sc)
    # solution (a list of paths): [[[240, 200], [190, 200], [190, 150], [240, 150]], [[200, 190], [230, 190], [215, 160]]]
