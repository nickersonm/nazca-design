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
#
# 2017-2022 (c) Ronald Broeke
#-----------------------------------------------------------------------

"""
This module is for configuring Nazca settings and sharing variables between
Nazca modules.
"""


import matplotlib.pyplot as plt


#==============================================================================
# Matplotlib layout settings
#==============================================================================
plt_cmap_name = 'Set2'  # colormap used in mask layout in Matplotlib
plt_alpha = 0.3  # tranparency
plt_figsize = 8  # size of layout
plt_fontsize = 10  # fontsize of annotations in the layout
plt_background_inside = '#FFFFFF'  # background color inside the axes
plt_background_outside = '#EEEEEE'  # background color outside the axes
plt_cmap = [] # list of the colors in the cmap
# Flag for redirecting matplotlib export to user defined Figure and Axes.
# To use, set plt_fig = (<figure>, <axes>), where <figure> and <axes> are one Figure and one Axes object from matplotlib
# e.g. <figure> and <axes> may come from functions as matplotlib.pyplot.subplots or  matplotlib.pyplot.subplot2grid.
# Default is None (required object will be created at runtime, and result displayed interactively).
plt_fig = None


# font setting for matplotlib generated mask output
matplotlib_font = {
    #'family': 'normal',
    'style' : 'normal',
    'weight': 'normal', #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
    'size'  : plt_fontsize}


# figsize for mode field plots
modeplotsize = (8, 5)


#==============================================================================
# Layout Settings
#==============================================================================
# The name of the cell that is generated automatically as the first topcell.
defaultcellname = 'nazca'


# add stubs to the cell bounding box
bbox_stubs = True


# Store pins inside the cell as annotation. One annotation per pin.
# Layer for output is 'bb_pin_attr'.
# This is a global setting applied to all pins created.
# default = False
store_pin_attr = False


# Store pins inside the cell annotation in a single annotation.
# Layer for output is default_layers['pin_text']
# This will be the default setting that can be overruled per cell:
# default = False
store_pins = False


# Store the cell name as cell annotation
store_bbname = True


# Redirect an unknown layer to the 'dump' layer.
# Set True for PDKs where only known layers should be allowed.
# Default = False works best for tutorials and examples
# and create layers by using them.
redirect_unknown_layers = False


# Log if Interconnect use non-user defined values for radius, width or euler
# This is useful in PDKs, where it should be set to True
# default = False
log_interconnect_defaults = False


# Allow for xsection without mask layers if True.
# A warning will be issued if False and a layer will be created automatically.
# Default = False
allow_empty_xsections = False


# Dictionary for default layers for special purpose objects in Nazca:
# They are defined in a dict so they can be adjusted during runtime
# and are easily accessible by 'key'.
# These layers will automatically be created when called for in a layer setting
# if not yet defined.
default_layers = {
    'bb_pin':             (1001, 0),  # layer for the pin symbol (arrow)
    'bb_pin_text':        (1002, 0),  # layer for pin name annotation
    'bb_pin_attr':        (1002, 1),  # layer for pin attributes annotations
    'bb_parameter_text':  (1003, 0),  # layer to store pcell parameters
    'bb_name':            (1004, 0),  # layer to store pdk-version and owner info: see bb_version_layer

    'package_pin':        (1011, 0),  # layer for package pin symbol
    'package_pin_text':   (1012, 0),  # layer for package pin annotation

    'tube_pin':           (1013, 0),  # layer for tube internal pins
    'tube_pin_text':      (1013, 1),  # layer for tube internal pin annotation
    'ribbon':             (1013, 2),  # layer for indication tube left and tube right side

    'bbox':               (1021, 0),  # layer for Nazca bounding box
    'bbox_name':          (1022, 0),  # layer for the name in the bbox
    'bbox_pin':           (1023, 0),  # layer for bounding box pins
    'bbox_pin_text':      (1024, 0),  # layer for bounding box pin annotations
    'userbbox':           (1021, 1),  # layer to (optionally) show the convex hull in.
    'hull':               (1021, 2),  # layer to (optionally) show the convex hull in.

    'pin_text':           (1005, 0),  # layer to store pin attribute information
    'docu_pin':           (6000, 0),  # Layer to store pin shapes for documentation purposes
    'dump':               (1111, 0),  # redirect content for unknown layers here
    'error':              (1111, 1),  # layer to display errors, e.g. to a bend with too small radius

    'bb_flag':            (1051, 0),  # layer to indicate a black box
    'validation_stub':    (1052, 0),  # layer for black box BB validation

    'drc_xs':             (1500, 0),  # display pin2pin drc on xsection
    'drc_angle':          (1500, 1),  # display pin2pin drc on angle
    'drc_width':          (1500, 2),  # display pin2pin drc on width
    'drc_radius':         (1500, 4),  # display pin2pin drc on radius
    'drc_sym':            (1500, 3),  # display pin2pin drc on width
    'drc_instance_angle': (1510, 0),  # display instantiation drc on angle
}


# Layer name to scan for BB versioning:
bb_version_layer = 'bb_name'

# default xs to create when an unknown xsection is set:
default_xs_name = 'nazca'
default_xs_layer = {
    'name': [default_xs_name],
    'layer': [default_layers['dump'][0]],
    'datatype': [default_layers['dump'][1]],
    'accuracy': [0.001],
    'growx': [0],
    'growy': [0],
    'origin': ['nazca']
}
default_xs_width = 1.0
default_xs_angle = 90.0
default_xs_radius = 20.0


# default error layer to place errorneous structures, e.g. impossible interconnects:
default_xserror_name = 'error'
default_xserror_width = 1.0
default_xserror_radius = 20.0

# Default xsections which will be automatically generated when a section name
# in the list (key) is used. The attributes (values) will be set for that xs.
default_xs_list = {
    default_xs_name: {
        'layer': 'nazca',
        'pinstyle': 'default',
    },
    'tube': {
        'layer': "tube_pin",
        'pinstyle': 'tube',
    },
    default_xserror_name: {
        'layer': 'error,',
        'pinstyle': 'default',
    },
}
# default technology
default_tech = None


# instantiate pin-shapes as separate cells.
# default=False
instantiate_pin = False


# instantiate stub cell as separate cells.
# default=False
instantiate_stub = False


# instantiate mask_elements as separate cells.
# default=False
instantiate_mask_element = False


# instantiate tube in interconnects as separate cells.
# default=None; leave it to internal interconnect setting
instantiate_tube = None


# instantiate ribbons in interconnects as separate cells.
# default=False
instantiate_ribbon = False


# Always use a ribbon, also for N=1 cases.
# default=False
always_use_ribbon = False


def instantiate_all(val=True):
    """Set all global intantiation flags at once."""
    global instantiate_pin, instantiate_stub, instantiate_mask_element
    global instantiate_tube, instantiate_ribbon
    instantiate_pin = val
    instantiate_stub = val
    instantiate_mask_element = val
    instantiate_tube = val
    instantiate_ribbon = val



# Add linenumbers to instance objects and gds
_linenumbers = False
def add_linenumbers():
    """Set linenumbers and instantiate mask element so they ill be tracked."""
    global _linenumbers, instantiate_mask_element, instantiate_ribbon, always_use_ribbon
    _linenumbers = True
    always_use_ribbon = True  # make sure all maskelements are in a ribbon in the end.
    instantiate_mask_element = False  # reduce the number of lower level instances for speed.
    instantiate_ribbon = True


#hull_mask_element = False


# use high-level layermapping in the lower-level gds routines.
# default = True
gds_over_getlayer = True


# Font default
# default = 'nazca'
default_font = 'nazca'


LDT2layername = {}


# maximum height of the cell name in gds in um
cellname_max_height = 50


# size of the cell_name in the cell w.r.t. the size of the cell
# default = 0.5
cellname_scaling = 0.5


# perform validate basenames on substrings if True.
validate_basename = False


# Connect pins based on proximity when placed.
# This enables pin2pin DRC and full circuit connectivity.
# It is possible to switch the group_connect on/off flag during a layout process,
# which can be useful to speed up the layout generation when needed.
# default = True
group_connect = True


# group_connect all pins based on proximity only, apply no filters.
# default = False
group_connect_all = False


# Maximum distance for autoconnecting pins when groupconnect = True
autoconnectdistance = 0.0014


# Angle in degrees from ideal chain-connect (180) that still will connect
# in group_connect. This helps to avoid "reflecting" paths on "splitting nodes".
# default = 175
group_connect_acceptance_angle = 175


# Connection logging for debugging only
# default = False
connections_logging = False


# Log file name for connection logging.
connections_logfile = 'connections.log'

# Log file handle for connection logging.
connections_logfile_handle = None


def log_connections():
    global connections_logging, connections_logfile, connections_logfile_handle
    connections_logging = True
    connections_logfile_handle = open(connections_logfile, 'w')


# Clear the mask after an export, or not.
# Examples:
# 1. Export multiple independent layouts from a cell -> export_clear = True
# 2. Incremental layout export, e.g. in a Jupyter Notebook -> export_clear = False
# default = True
export_clear = True


# Default setting for cells to add convex hull polygon to the mask export.
# This only has effect if True AND the cell's hull setting is hull is True.
# default = False
show_hull = False


# Default hull setting for all Cells
# Set True to use a convex hull of instances in the cell
# for bounding box calculations.
# A convex hull gives an instance a "rotation proof" bbox.
# The hull may cause the the layout generation to take several times longer.
# Note that with use_hull = False for flattened instances the bbox of
# the parent cell may look "too large" as it is based on
# the corners of rotated bbox if the instance that is not visualized.
# default = False
use_hull = False


# Use hull of instances (if they have use_hull is True) as input to
# the bbox or hull calculation.
# default = False
hull_based_bbox = False


# Generate polygons and polylines in mask elements
# Default is True, do not change.
generate_shapes = True


# logging level to include source coe linenumbers in log files.
# 10 DEBUF, 20 INFO, 30 WARNING, 40 ERROR, 50 CRITICAL
# default = 30 for WARNING and up.
logging_level_linenumber = 30

# Flag to raise an exception on netlist connection errors
# This helps identifying where the exception occurs in the code.
# default = False
netlist_raise = False


# Flag to raise an exception on interconnect warnings
# This helps identifying where the exception occurs in the code.
# default = False
interconnect_raise = False

# Flag to raise an exception on the cummulative error counter
# This helps identifying where the exception occurs in the code.
# default = False
total_raise = False


# Switch rangechecking on/off (True/False)
# Switching it off can be convenient to bypass a "wrong" blocking range exception.
# Do not forget to set it to True again after the bypass.
rangecheck = True


# Show unknown GDS records in logging when loading a GDS in stdout.
# default = False
show_unknown_GDS_records = False


# Style of pcell parameter annotation in ip-blocks/cells
# options ["default", 'hhi', 'yaml']
parameter_style = 'default'


# Set tubepin naming
tubepin_prefix = "tube"


# Convert polylines into polygons in gds export
# Note this conversion happens automatically for formats without polylines, e.g. png.
# default = False
export_polyline_as_polygon = False
# GDS standard defines 600 points as maximum. This seems unnecessarily few points.
maxpolygonpoints = 4000  # Safe maximum for most software.
# GDS standard defines 200 points as maximum. This seems unnecessarily few points.
maxpolylinepoints = 4000  # Safe maximum for most software.


#==============================================================================
# DRC visualisation related variables
#==============================================================================
# pin2pin DRC layers:
drc_layer_xs     = 'drc_xs'
drc_layer_angle  = 'drc_angle'
drc_layer_width  = 'drc_width'
drc_layer_radius = 'drc_radius'
drc_layer_sym    = 'drc_sym'

drc_layer_instance_angle = 'drc_instance_angle'

# DRC gds ring dimensions:
drc_ring_xs    = (10, 1) # symbol displayed on xsection drc location
drc_ring_angle = (9, 1)
drc_ring_width = (8, 1)
drc_ring_width = (7, 1)
drc_ring_sym   = (6, 1)

drc_rule_xs = {}  # dict map which xsections are allowed to be connected
drc_raise = False  # raise a drc error for trace-back the error if True
drc = False  # perform drc on pin connections if True

drc_instance_angle = {'angle': {}}


# In case of a geometrically discontinous pin2pin connection
# (as determined by drc_max_distance) the xs and width of a new node will
# be cleared (set to None) if this parameters is True.
# Clearing will supress pin2pin drc on 'width' and 'xs'.
# Setting this parameter to true requires explicit xs and width settings in
# move method to clear DRC errors.
# default = False
clear_pin_on_move = True

# Maximum distance between node to be considered a smooth/continuous node connection
drc_max_distance = 1e-4

# Maximun angle to be considered smooth (in degrees)
drc_max_angle = 5e-6


# =============================================================================
# uPDK
# =============================================================================
allowed_vartypes = ['int', 'str', 'float', 'bool', None]
allowed_units = ['um', '']

# =============================================================================
# initialize global variables into existence
# =============================================================================
cp = None
self = None  # active cell reference
cells = []  # store open cells
cellnames = dict()  # list of all cells {cellname: Cell}
basenames = dict()  # {basenames" function_id}
topcells = []  # Cells that have been collected to be parsed for GDS export
xsall = dict()  # foundry layer table
xs_layers = dict()  # foundry layer table
#layerdict = dict()
xsmap = None  # map scriptnames of xs to technology names of xs.
share = set()  # set of cellname not to prefix.
stubmap = {}
XSdict = dict()
patchcell = False  # Keep set to False at all times.
reuse_pathlayers = set() # Allow the pathfinder method to reuse an existing path name by storing names allowed for reuse.


# =============================================================================
# bb_util.hashme
# =============================================================================
hashme = False  # for checking @hashme status
# Default parameter representation settings:
hashme_config = {
    "prefix":  "_",    # Prefix for variable
    "format":  ".5g",  # Format of variable (None redirects to cfg.hashme_default_format)
    "func":    None,   # option to provide a custom formatting function
    "suffix":  "",     # default value for suffix
    "addhash": True,   # default value for addhash
}
# Check if the function id using a certain cell base name is unique.
# Check happens in hashme.
# default = True
check_basename = True


# =============================================================================
# define pin shape polygons in a dictionary. Idealy shapes are normalized to 1.
# =============================================================================
pinshapes = {
    'arrow_full_asym': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.5, 0.0), (-0.7, 0.0), (-0.7, 0.25), (-0.5, 0.25), (-0.5, 0.5)],
    'arrow_full':      [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25), (-0.7, 0.25), (-0.5, 0.25), (-0.5, 0.5)],
    'arrow_left':      [(0, 0), (-0.5, 0.5), (-0.5, 0.25), (-0.7, 0.25), (-0.7, 0)],
    'arrow_right':     [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25), (-0.7, 0)],
    'circle':          [(0.5, 0), (0.353, 0.353), (0, 0.5), (-0.353, 0.353), (-0.5, 0), (-0.353, -0.353), (0, -0.5), (0.353, -0.353)],
    'inv_pointer':     [(0, 0), (-0.25, 0.5), (-0.5, 0.5), (-0.5,-0.5), (-0.25, -0.5)],
    'arrow_head':      [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.25, 0.0), (-0.5, 0.25), (-0.5, 0.5)],
    'box1':            [(0, 0.25), (-0.5, 0.25), (-0.5, -0.25), (0, -0.25)],
    'box3':            [(0, -0.5), (0, 0.5), (-0.5, 0.5), (-0.5, -0.5), (0, -0.5), (-1, -0.5), (-1, 0.5), (-0.5, 0.5), (-1.5, 0.5), (-1.5, -0.5), (-1, -0.5)],
}


# list of recognized pinstyles to validate new pinstyle requests against
# Note that a pinstyle includes stub parameters and other info too.
pinstylelabels = [
    # pinshape cell settings:
    'shape',  # "pinshapes" label to access pinshape polygon
    'size',  # size of the pin, default = 1.0
    'layer',  # string name of layer to place the pinshape polygon in

    # annotation setttings:
    'annotation_layer',  # string name of layer to place the pin name annotation in
    'annotation_move',  # tuple (x, y) distance of pin annotation from actual pin location

    # stub display settings:
    'stub_length',  # length of the stub in um, set to 0 to not show any stub layers.

    # documentation settings:
    'scale',  # scaling of the pin for png export
    'pin_ignore',  # list pin names to not visualize in png export.
    'edgecolor',  #  pinshape edge color in png export
    'fillcolor',  #  pinshape fill color in png export

]

def reset_pinstyles():
    """Define default pinstyle dicts."""
    global pinstyle, pinstyles, default_pinstyle

    # pin connections:
    pinstyle = {
        'shape': 'arrow_full',
        'size': 1.0,
        'layer': 'bb_pin',
        'annotation_layer': 'bb_pin_text',
        'annotation_move': (-0.3, 0),
        'stub_length': 2.0,
        'scale': 1/15,
        'pin_ignore': []
    }

    # bbox pins (abutment pins)
    pinstyle_bbox_left = {
        'shape': 'arrow_left',
        'size': 0.5,
        'layer': 'bbox_pin',
        'annotation_layer': 'bbox_pin_text',
        'annotation_move': (-0.3, 0),
        'stub_length': None
    }
    pinstyle_bbox_right = {
        'shape': 'arrow_right',
        'size': 0.5,
        'layer': 'bbox_pin',
        'annotation_layer': 'bbox_pin_text',
        'annotation_move': (-0.3, 0),
        'stub_length': None
    }
    pinstyle_bbox_center = {
        'shape': 'arrow_full',
        'size': 0.5,
        'layer': 'bbox_pin',
        'annotation_layer': 'bbox_pin_text',
        'annotation_move': (-0.3, 0),
        'stub_length': None
    }

    # documenatation:
    pinstyle_docu = {
        'shape': 'arrow_full',
        'scale': 1/15,
        'pin_ignore': [],
        'edgecolor': "#000000",
        'fillcolor': "#FFFFFF",
    }

    # tube connectivity support:
    pinstyle_tube = {
        'shape': 'arrow_head',
        'size': 0.7,
        'layer': 'tube_pin',
        'annotation_layer': 'tube_pin_text', # use string name
        'annotation_move': (-0.3, 0),
        'stub_length': 0.0,
    }

    # ribbon order support:
    pinstyle_ribbon0 = {
        'shape': 'box1',
        'size': 1.0,
        'layer': 'ribbon',
        'stub_length': 0.0,
    }
    pinstyle_ribbonN = {
        'shape': 'box3',
        'size': 1.0,
        'layer': 'ribbon',
        'stub_length': 0.0,
    }

    # add pin styles to a global container to be referenced by string name:
    pinstyles = {
        'default':      pinstyle.copy(),
        'bbox_left':    pinstyle_bbox_left,
        'bbox_right':   pinstyle_bbox_right,
        'bbox_center':  pinstyle_bbox_center,
        'docu_default': pinstyle_docu,
        'tube':         pinstyle_tube,
        'ribbon0':      pinstyle_ribbon0,
        'ribbonN':      pinstyle_ribbonN,
    }

    default_pinstyle = 'default'


default_pinstyle_name = 'default'
reset_pinstyles()
pinstyles_overrule = {}  # set before loading a foundry to update pin representation
documentation_pins = False

# For pytesting it may be usefull to not have the version annotation change
# in this case force another dict content in your code, e.g. {}
# default = None  # no version forcing
force_pdk_version = None
force_nazca_version = None


# For black box export with simplified layers and pins for validation
# default = None  # add xs name string to overrule the existing stub
validation_stub_xs = None

# For black box export with simplified layers and pins for validation.
# Use pinstyle name string to overrule the existing pin style
# default = None
validation_pinstyle = None

# For black box export with simplified layers and pins for validation
# default = "none"  # no layer filtered
validation_layermapmode = None

# default = None  # No change w.r.t. layermapmode setting
validation_layermap = None


def active_cell():
    """Get the active cell.

    The 'active cell' is the last opened Cell object that has not been closed yet.

    Returns:
        Cell: active cell
    """
    if cells:
        return cells[-1]
    else:
        return None

def gds_cell_validate(cell):
    """Placeholder function that operates on a cell at Cell.close.

        Args:
            cell (cell): The cell which will be operated on.

        Returns:
            bool: None
        """
    return None

def gds_cellname_cleanup(name):
    """Placeholder function that can be overruled in a foundry PDK.

    This function should make sure only valid cellnames are used, i.e.
    formats excepted by the foundry.

    Args:
        name (str): cellname

    Returns:
        str: valid cellname
    """
    return name


def formatplot():
     """update matplotlib plotting style for graphs.

     Returns:
         None
     """
     lines = {'linewidth':3}
     plt.rc('lines', **lines)
     font = {'size': plt_fontsize}
     plt.rc('font', **font)


def save_mode():
    """Switch on all relevant cfg flags to save mode.

    Use for PDKs where you do not want unnoticed definitions or default settings.

    Returns:
        None
    """
    global redirect_unknown_layers, log_interconnect_defaults
    redirect_unknown_layers = True
    log_interconnect_defaults = True
