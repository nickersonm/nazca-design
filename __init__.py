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
# Ronald Broeke (c) 2016-2019
# Main initialization


"""Nazca Design module for photonic integrated circuit design."""

#__all__ = ['Cell', 'Pin', 'Polygon', 'Polyline', 'Annotation', 'load_gds',
#    'export_gds', 'export_plt', 'add_layer', 'get_layer', 'show_layers',
#    'show_xsections', 'layout_open', 'cell_open', 'cell_close', 'layout_close',
#    'gds_filter', 'replaceCells', 'print_structure', 'image',
#    'strt', 'bend', 'taper', 'ptaper', 'cp', 'add_xs', 'text', 'linelength',
#    'textheight', '__version__', 'modulepath']

import pandas as pd

from .version import __version__
from . import cfg
from .simglobal import sim
#from .slabmode import *
#from .slabsolver import *
from .xsection import *
from .plotting import *
from .gds import *
from .util import *
from .extras import *
from .convert_AWGtxt2nazca import *
from .gds_import import *
from .gds_base import *
from .colormap import *

from .netlist import *
from .layout import *
from .pathfinder import *
from .index import *

from . import geometries as geom
from . import icons as icons
from .mask_layers import *
from .mask_elements import *
from .gds_cellreplace import *
from .font import Font, text, linelength, default_font

from . import interconnects as ic
from . import cp
from .logging import *
from .angle_drc import *

from .pdk_template import *
from .polysplitter import *

# solvers
from .solver_slab2020 import *
from .solver_multilayerc import *



# Nazca module path
import os
modulepath = os.path.dirname(os.path.abspath(__file__))
