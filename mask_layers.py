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
#
# @author: Ronald Broeke (c) 2016-2019
# @email: ronald.broeke@brightphotonics.eu
#-----------------------------------------------------------------------
# -*- coding: utf-8 -*-
"""
This module defines
1- layers in DataFrame 'layer_table'
2- xsections in DataFrame 'xsection_table.
3- xsection layers in DataFrame 'xsection_layer table', adding layers to each 'xs'
4- The 'layer_table' and 'xsection_layer table' are joined using the 'layer2xsection table'.

The result is stored in a dictionary containing the xsection-layer information
per 'xs' needed in layout: 'xsection_layers'

Layers and xsection can be added in two ways
1 - Loading from csv file: load_layers, load_xs
2 - Functions add_layer, add_xsection
When adding layers to a 'xs' the xs_map is built automatically.

It is also possible to first load layers and xs, and then extend them via
'add_layer' and/or 'add_layer2xsection'.
"""

from collections import defaultdict, namedtuple
import pandas as pd
import numpy as np
from pprint import pprint
from math import isnan
import matplotlib as mpl
from matplotlib.colors import rgb2hex

from . import cfg
from . import xsection
from nazca.logging import logger
import nazca as nd
#from nazca.mask_elements import _handle_missing_xs


# keep track of layers used by the designer that have not been defined:
unknown_layers = set()
unknown_xsections = set()
#cfg.layerset = set() #set of all (L, D) layers defined

cfg.layername2LDT = {}
cfg.layername2LDT_cnt = defaultdict(int)

cfg.xsection_table = pd.DataFrame()
xsection_table_attr = ['xsection', 'xsection_foundry', 'origin', 'stub', 'description']


#layer-table:
cfg.layer_table = pd.DataFrame() # layer_name -> index
layer_table_attr_csv  = ['layer', 'datatype', 'tech', 'layer_name_foundry',
    'accuracy', 'origin',       'remark', 'description']
#rename attributes for the xs-layer join to avoid name clashes:
layer_table_attr_join = ['layer', 'datatype', 'tech', 'layer_name_foundry',
    'accuracy', 'origin_layer', 'remark', 'description']


#layercolor table:
cfg.colors = pd.DataFrame()
cfg.xs_list = []
colors_attr_csv = ['name', 'layer', 'datatype', 'depth', 'fill_color', 'frame_color', 'frame_brightness', 'fill_brightness', 'dither_pattern', 'valid', 'visible', 'transparent', 'width', 'marked', 'animation']


#xsection-layer-map table, mapping layers to xsections:
cfg.xsection_layer_map = pd.DataFrame()

xsection_layer_map_attr_csv  = ['layer_name', 'xsection', 'growx', 'leftedgefactor', 'leftedgeoffset', 'rightedgefactor', 'rightedgeoffset', 'growy', 'growy1', 'growy2', 'polyline' , 'origin']
#rename attributes for the xs-layer join to avoid name clashes:
xsection_layer_map_attr_join = ['layer_name', 'xsection', 'growx', 'leftedgefactor', 'leftedgeoffset', 'rightedgefactor', 'rightedgeoffset', 'growy', 'growy1', 'growy2', 'polyline', 'origin_xs']


#final set of export attributes for joined xsection-layers table:
#mask_layers_attr = ['layer_name', 'xsection', 'layer', 'datatype', 'tech', 'growx', 'growy', 'accuracy']
mask_layers_attr = ['layer_name', 'xsection', 'layer', 'datatype', 'tech', 'growx', 'leftedgefactor', 'leftedgeoffset', 'rightedgefactor', 'rightedgeoffset', 'growy', 'growy1', 'growy2', 'accuracy', 'polyline']


layer_debug = False # flag for debug output to stdout.
block_load_colors = False
None_layerflag = False  # register if layer=None message has been issued


def _check_duplicates(df, tablename):
    """Find duplicated layers entries in a layer table.

    Obsolete

    Duplicates looks for ['layer', 'datatype', 'tech'].

    Return:
        bool: True if duplicate layers are found.
    """
    ID = ['layer', 'datatype', 'tech']
    G = df.duplicated(subset=ID, keep=False)
    if G.any():
        if layer_debug:
            print("\nDuplicate layers found in '{}':".format(tablename))
            print("\nSame ['layer', 'datatype', 'tech'], but different layer_name.")
            print(df[['layer_name']+ID].sort_values(by=ID)[G])
        return True
    return False


def add_layer(
    name=None,
    layer=None,
    tech=None,
    accuracy=None,
    fab_name=None,
    origin=None,
    remark=None,
    description='',
    overwrite=False,

    # color info:
    frame_color=None,
    fill_color=None,
    frame_brightness=None,
    fill_brightness=None,
    dither_pattern=None,
    valid=None,
    visible=None,
    transparent=None,
    width=None,
    marked=None,
    animation=None,
    alpha=None,

    # flow control:
    unknown=False,
    merge=True,
):
    """Create a new mask layer.

    Args:
        name (str): layer name
        layer (int | tuple): layer number or (layer, datatype)
        tech (str): technology ID
        accuracy (float): mask resolution of the layer in um
        fab_name (str): optional name used by the foundry for this layer
        origin (str): who created the layer
        remark (str): extra info for the layer
        overwrite (bool): overwrite an existing layer if True (default=False)
        unknown (bool): True if redirected from get_layer (default = False)
        merge (bool): if True merge xscetion and layers (default is True)
        ...: color attributes
        skip_color: skip setting of color attributes if True

    Returns:
        DataFrame: Table of mask layers
    """
    global None_layerflag

    isDefault = False
    if unknown:
        # Redirected from get_layer -> layer unknown.
        #   - If the layer corresponds to a default_layer in cfg then
        #     create it here as an actual layer.
        #   - If a layer is redirected to 'dump' then raise a warning (first time it occurs)
        LDT_or_name = cfg.default_layers.get(name, False)
        if LDT_or_name:
            if isinstance(LDT_or_name, str): # for ID (str) mapped to ID (str)
                _layer, _datatype, _tech = cfg.layername2LDT[LDT_or_name]
            elif len(LDT_or_name) == 2:
                _layer, _datatype = LDT_or_name
                tech = cfg.default_tech
            else:
                _layer, _datatype, _tech = LDT_or_name
            isDefault = True
        elif layer is None:
            isDefault = False
        else:  # map number to layerID
            _layer, _datatype, tech = layer
            newML = (_layer, _datatype)
            if newML in cfg.default_layers.values():
                for MLname, ML in cfg.default_layers.items():
                    if ML == newML:
                        name = MLname
                        isDefault = 'by_value'
                        break

        if isDefault is False:
            # create a new layer or redirect to 'dump'
            if cfg.redirect_unknown_layers: # do not auto-create unknown layers.
                if layer not in unknown_layers:
                    # TODO: user option to raise exception to correct the error.
                    msg = f"Unknown layer {layer} in get_layer()."\
                        f" Moving layout content to the 'dump 'layer 'dump' instead,"\
                        f" because redirect_unknown_layers = True."
                    if isinstance(cfg.gdsload, str):
                        msg += f" Layer found while loading '{cfg.gdsload}'."
                    nd.main_logger(msg, "warning")
                    unknown_layers.add(layer)
                return get_layer('dump')
            elif layer is None:
                if not None_layerflag:
                    nd.main_logger(f"Using layer=None. Setting it to 'dump' layer instead. "
                        f"Alternatively create this layer explicitly using:\n"
                        f"add_layer(name='{name}', layer=(<layerinfo>))\n"
                        f"This message is only shown once.",
                        "warning"
                    )
                    None_layerflag = True
                return get_layer('dump')

    else:  # not unknown (not a redirect from get_layer())
        if tech is None:
            tech = cfg.default_tech
        if isinstance(layer, tuple):
            if len(layer) == 2:
                _layer, _datatype = layer[0], layer[1]
            elif len(layer) == 3:
                _layer, _datatype, tech = layer[0], layer[1], layer[2]
            else:
                mes = "Invalid layer value: {}".format(layer)
                logger.exception(mes)
                raise Exception(mes)
        else:
            _layer = layer
            _datatype = 0

    # generic layer handling:
    if _layer is None:
        mes = "Error: no layer number provided in call to add_layer()."
        logger.exception(mes)
        raise Exception(mes)
    elif _layer != int(_layer):
        mes = "layer value '{}' should be of type integer".format(_layer)
        logger.exception(mes)
        raise Exception(mes)
    if _datatype != int(_datatype):
        mes = "datatype value '{}' should be of type integer".format(_datatype)
        logger.exception(mes)
        raise Exception(mes)

    if accuracy is None:
        accuracy = cfg.default_xs['accuracy'][0]
    if origin is None:
        origin = cfg.default_xs['origin'][0]
    if name is None:
        name = '{}/{}/{}'.format(_layer, _datatype, tech)

    if name in cfg.layername2LDT.keys() and not overwrite: # existing layer_name (ID):
        if isDefault == 'by_value':
            if (_layer, _datatype, tech) != cfg.layername2LDT[name]:
                mes = f"Provided layer ({_layer}, {_datatype}) maps to an already used layer ID '{name}'"
                logger.warning(mes)
        else:
            mes = "Reusing an existing layer_name '{}'."\
            " Use a different layer name or set overwrite=True.".format(name)
            logger.warning(mes)
    else:
        LDT = (_layer, _datatype, tech)
        cfg.layername2LDT[name] = LDT
        cfg.layername2LDT_cnt[LDT] += 1
        if cfg.layername2LDT_cnt[LDT] == 1:
            cfg.LDT2layername[LDT] = name

        layer_name = str(name)
        newlayer = {
            'layer_name': layer_name,
            'layer': [_layer],
            'datatype': [_datatype],
            'tech': tech,
            'accuracy': [accuracy],
            'layer_name_foundry': fab_name,
            'origin': [origin],
            'remark': remark,
            'description': description,
            'source': 'add_layer'}

        df = pd.DataFrame(newlayer)
        df.set_index('layer_name', inplace=True)
        try: #pandas in python 3.7
            cfg.layer_table = pd.concat([cfg.layer_table, df], sort=True)
        except:  #pandas in python <= 3.6
            cfg.layer_table = pd.concat([cfg.layer_table, df])

        set_layercolor(layer_name,
            frame_color, fill_color, frame_brightness, fill_brightness,
            dither_pattern, valid, visible, transparent, width, marked, animation,
            alpha)

        merge_xsection_layers_with_layers()
    return name


def add_XSdict(name, name_foundry=None, origin='add_xsection', stub=None, description=''):
    """Create a new Xsection object named <name>.

    If a Xsection with <name> already exists, the existing Xsection is returned.

    Args:
        name (str): xsection name.
        name_foundry (str): optional xsection name used by the foundry
        origin (str): source/creator of the layer, e.g. 'user', 'foundry' or any string.
        stub: xsection name to be used for the stub of this xsection

    Returns:
        Xsection: Xsection object with name <name>
    """
    if name not in cfg.XSdict.keys():
        XS = xsection.Xsection(name=name)
        cfg.XSdict[name] = XS
        XS.description = description
        XS.origin = origin
        XS.stub = stub

        # add new info to xsection_layers:
        D = {
            'xsection': name,
            'xsection_foundry': name_foundry,
            'origin': origin,
            'description': description,
        }

        if cfg.xsection_table.empty:
            cfg.xsection_table = cfg.xsection_table.append(D, ignore_index=True)
        elif cfg.xsection_table[cfg.xsection_table['xsection'] == name].empty:
            cfg.xsection_table = cfg.xsection_table.append(D, ignore_index=True)

        #update mask_layers:
        merge_xsection_layers_with_layers()
    else:
        XS = get_xsection(name)
        if description != '':
            XS.description = description
        if stub is not None:
            XS.stub = stub
        if origin != '':
            XS.origin = origin

    return cfg.XSdict[name]

add_xsection = add_XSdict


def parse_grow(edge, side, default):
    """Parse growth parameter into a tuple.

    Returns:
        tuple: (a, b) grow information as a*width+b
    """
    if edge is None:
        return (side, default)
    if isinstance(edge, (float, int)):
        return (side, float(edge))
    if not isinstance(edge, tuple):
        logger.error("Growth value for layer in xsection should be a tuple like (width_factor, shift)."\
            " Using value (1.0, 0.0) instead of given {}".format(edge))
        return (1.0, 0.0)
    else:
        return edge


def add_layer2xsection(xsection='name', layer=None,
        growx=None, leftedge=None, rightedge=None,
        growy=None, growy1=None, growy2=None,
        accuracy=None, fab_name=None, origin=None,
        overwrite=False, polyline=False, remark=None):
    """Add a layer to the Xsection object with name <xsection>.

    If <xsection> does not exist yet it will be created.
    The layer entry is described by (layer, datatype).
    If the (layer, datatype) does not exist yet it will be created.
    The keyword 'accuracy' is only processed in case of a new layer.

    The "leftedge" and "rightedge" keywords describe the the waveguide sidewalls
    by the distance of the guide sidewalls. The waveguide in this concept is oriented
    on the direction of the positive x-axis and the sidewalls are expressed as
    the y-coordinate as a*width+b, or in short (a, b), where width is the width
    of the guide. For a standard waveguide we have leftedge=(0.5, 0.0)
    and rightedge = (-0.5, 0.0). Swapping the values of leftedge and rightedge
    generate the same polygon, so the naming is indicative only.

    NOTE: growy (or growy1 and growy2 for that matter) will only have an
    effect on the straight guide and the linear taper.

    Args:
        xsection (str): xsection name for use in nazca
        layer (int | tuple): layer number or (layer, datatype)
        growx (float): growth in x-direction of polygon in the (layer, datatype)
            (default = 0.0)
        growy (float): growth in y-direction of polygon in the (layer, datatype)
            (default = 0.0)
        accurary (float): accuracy of the grid in um of the (layer, datatype)
            (default = 0.001)
        fab_name (str): technology xsection name, e.g. the foundry name for the layer
        origin (str): source/creator of the layer, e.g. 'user', 'foundry'
        overwrite (bool): Overwrite a layer if it already exists (default=False)
        polylilne (bool); Render the structures a polyline rather than polygon (defaul=False)
        remark (str): extra info about the layer
        leftedge (float | (float, float)): b or (a, b) for an edge at a*width+b, a=0.5 by default
        rightedge (float | (float, float)): b or (a, b) for an edge at a*width+b, a=-0.5 by default

    Returns:
        Xsection: Xsection object having name <xsection>
    """
    layer = get_layer(layer)
    if growx is None:
        growx = cfg.default_xs['growx'][0]
    if growy is None:
        growy = cfg.default_xs['growy'][0]
    leftedgefactor, leftedgeoffset = parse_grow(leftedge, 0.5, growx)
    rightedgefactor, rightedgeoffset = parse_grow(rightedge, -0.5, -growx)
    if growy1 is None:
        growy1 = growy
    if growy2 is None:
        growy2 = growy

    if origin is None:
        origin = cfg.default_xs['origin'][0]
    if accuracy is None:
        accuracy = cfg.layer_table.loc[layer]['accuracy']
    if xsection not in cfg.XSdict.keys():
        add_xsection(name=xsection, origin='add_layer2xsection')

    if not cfg.xsection_layer_map.empty:
        select = (
            (cfg.xsection_layer_map['xsection']==xsection) &
            (cfg.xsection_layer_map['layer_name']==layer)
        )
        LayerRow = cfg.xsection_layer_map[select]
        if not LayerRow.empty and overwrite: # overwrite layer entry
            cfg.xsection_layer_map.loc[
                select,
                ['growx',
                 'leftedgefactor',
                 'leftedgeoffset',
                 'rightedgefactor',
                 'rightedgeoffset',
                 'growy',
                 'growy1',
                 'growy2',
                 'polyline']
            ] = (growx,
                leftedgefactor,
                leftedgeoffset,
                rightedgefactor,
                rightedgeoffset,
                growy,
                growy1,
                growy2,
                polyline
            )
            cfg.xs_list = cfg.xsection_layer_map['xsection'].unique()
            merge_xsection_layers_with_layers()
            return cfg.XSdict[xsection]
    else:
        pass # first xs is a but random choice and confusing when using nd.strt() type of elements.
        # cfg.default_xs_name = xsection # set default xs to first user defined xs.

    #define new layer entry
    D = {
        'layer_name': str(layer),
        'xsection': xsection,
        'xsection_foundry': fab_name,
        'growx': growx,
        'leftedgefactor': leftedgefactor,
        'leftedgeoffset': leftedgeoffset,
        'rightedgefactor': rightedgefactor,
        'rightedgeoffset': rightedgeoffset,
        'growy': growy,
        'growy1': growy1,
        'growy2': growy2,
        'accuracy': accuracy,
        'origin': origin,
        'remark': remark,
        'polyline': polyline,
        'source': 'add_layer'
    }
    cfg.xsection_layer_map = cfg.xsection_layer_map.append(D, ignore_index=True)

    merge_xsection_layers_with_layers()
    return cfg.XSdict[xsection]

add_xsection_layer = add_layer2xsection


def load_layers(filename, tech=None, clear=False, autocolor=False):
    """Load layer definitions from csv file into the Nazca layer table.

    The load can replace existing layer definitions or extend them (default).

    Expected csv file format (no spaces after thhe comma!):

        layer_name, layer_name_foundry, layer, datatype, accuracy, origin, remark
        name1,name_alias,2,10,0.001,nazca,this is a remark

    Columns may be omitted but at the least use

        layer_name,layer,datatype,accuracy
        name1,2,10,0.001

    layer_name: layer name used in Nazca
    layer_name_foundry: layer name as used by the foundry, if applicable
    layer: layer number
    datatype: datatype number
    accuracy: mask layer resolution in um
    origin: indicates who defines the layer (foundry, nazca, you?)
    remark: optional information on the layer

    Args:
        filename (str): name of layer definition file in csv format
        tech (str): technology name
        clear (bool): clear existing layers (default = False)
        autocolor (bool): Create automated color settings per layer (default = False)

    Returns:
        DataFrame: table with layers
    """
    if tech is None:
        tech = cfg.default_tech
    load_table = pd.read_csv(filename, delimiter=',')
    if 'tech' not in load_table:
        load_table['tech'] = tech
        if layer_debug:
            print("Column 'tech' missing in file '{}'."\
                " Column has been added after load."\
                " For explicitly use of 'tech' add it to the table".format(filename))
    if 'origin' not in load_table:
        load_table['origin'] = ''
    if 'remark' not in load_table:
        load_table['remark'] = ''
    if 'description' not in load_table:
        load_table['description'] = ''
    else:
        load_table.description = load_table.description.fillna('')

    load_table.dropna(subset=['layer_name'], inplace=True)
    load_table['layer'] = load_table['layer'].astype(int) # note this gets you dtype numpy.int64
    load_table['datatype'] = load_table['datatype'].astype(int)
    load_table['source'] = filename

    if clear:
        clear_layers()
    for i, row in load_table.iterrows():
        add_layer(
            name=row.layer_name,
            layer=(row.layer, row.datatype),
            accuracy=row.accuracy,
            origin=row.origin,
            remark=row.remark,
            description=row.description,
            merge=False)
            # merge only after adding multiple layers to save time
            # add_layer() check if the layername is unique
    merge_xsection_layers_with_layers()

    if autocolor:
        for i, (name, L, D, T) in load_table[['layer_name', 'layer', 'datatype', 'tech']].iterrows():
            set_layercolor(name, (L, D, T),
                frame_color=None, fill_color=None, frame_brightness=None,
                fill_brightness=None, dither_pattern=None, valid=None,
                visible=None, transparent=None, width=None, marked=None,
                animation=None, alpha=None)

    return cfg.layer_table


def load_xsection_layer_map(filename, tech=None):
    """Load the assignment of layers to xsections from <filename>.

    Layers or xsections in the file that do not yet exist are created automatically.

    Args:
        filename (str): xsection file in csv format
        tech (str): technology ID.

    Returns:
        DataFrame: merged xsections and layers.
    """
    if tech is None:
        tech = cfg.default_tech
    cfg.xsection_layer_map = pd.read_csv(filename, delimiter=',')

    cfg.xsection_layer_map = cfg.xsection_layer_map[pd.notnull(cfg.xsection_layer_map['xsection'])]
    cfg.xsection_layer_map = cfg.xsection_layer_map.where((pd.notnull(cfg.xsection_layer_map)), None)
    # cfg.xsection_layer_map.dropna(inplace=True)

    if 'tech' not in cfg.xsection_layer_map.columns:
        cfg.xsection_layer_map['tech'] = tech
        if layer_debug:
            logger.warning("Column 'tech' missing in file '{}'."\
                "Column has been added after load"\
                " with value {}".format(filename, tech))
    if 'leftedgefactor' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['leftedgefactor'] = None
    if 'leftedgeoffset' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['leftedgeoffset'] = None
    if 'rightedgefactor' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['rightedgefactor'] = None
    if 'rightedgeoffset' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['rightedgeoffset'] = None
    if 'growy1' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['growy1'] = None
    if 'growy2' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['growy2'] = None
    if 'polyline' not in cfg.xsection_layer_map.columns:
       cfg.xsection_layer_map['polyline'] = False

    for i, row in cfg.xsection_layer_map.iterrows():
        if row.leftedgefactor is None:
            cfg.xsection_layer_map.at[i, 'leftedgefactor'] = parse_grow(None, 0.5, row.growx)[0]
        if row.leftedgeoffset is None:
            cfg.xsection_layer_map.at[i, 'leftedgeoffset'] = parse_grow(None, 0.5, row.growx)[1]
        if row.rightedgefactor is None:
            cfg.xsection_layer_map.at[i, 'rightedgefactor'] = parse_grow(None, -0.5, -row.growx)[0]
        if row.rightedgeoffset is None:
            cfg.xsection_layer_map.at[i, 'rightedgeoffset'] = parse_grow(None, -0.5, -row.growx)[1]
        if row.growy1 is None:
            cfg.xsection_layer_map.loc[i, 'growy1'] = row.growy
        if row.growy2 is None:
            cfg.xsection_layer_map.loc[i, 'growy2'] = row.growy

    # add xsections in the map to the xsection_table if missing:
    for name in cfg.xsection_layer_map['xsection']:
        if name not in cfg.XSdict.keys():
            add_XSdict(name, origin='xs2layer')

    return merge_xsection_layers_with_layers()


def load_xsections(filename):
    """Load list of xsection from file with stub mapping.

    In addition create a stub map dictionary of a xsection to its stub-xsection.

    Args:
        filename (str): xsection map filename in csv format

    Returns:
        DataFrame: table with loaded xsections
    """
    cfg.xsection_table = pd.read_csv(filename, delimiter=',')
    cfg.xsection_table.dropna(how='all', inplace=True)

    #temptable = cfg.xsection_table.set_index('xsection')
    if 'description' not in cfg.xsection_table:
        cfg.xsection_table['description'] = ''
    cfg.xsection_table['description'] = cfg.xsection_table['description'].fillna('')
    cfg.stubmap = cfg.xsection_table.set_index('xsection')['stub'].dropna().to_dict()

    for i, row in cfg.xsection_table.iterrows():
        add_XSdict(name=row.xsection, description=row.description,
            origin=row.origin, stub=row.stub)

    #TODO: check naming consistency of stubs
    #TODO: check for unique names
    #TODO: fill missing nazca_name with fab_name?
    #TODO: check if all xs are present in the xs definition.

    #update mask_layers:
    merge_xsection_layers_with_layers()
    return cfg.xsection_table


def add_xsection_stub(xsection, stub):
    """Assign a stub xsection to a xsection.

    The default stub is the xsection itself.
    """
    cfg.stubmap[xsection] = stub


def merge_xsection_layers_with_layers():
    """Create a dictionary containing tables of xsections joined with layers.

    The xsection names are the dictionary keys.
    The dictionary contains all attributes for mask  rightedgefactor, rightedgeoffsetexport.

    Left-join xsection table with layer table on (layer, datatype) and
    chop it up into a dictionay with xsection as key and the join result as
    values.

    Returns:
        None
    """
    if cfg.layer_table.empty or\
        cfg.xsection_layer_map.empty or\
        cfg.xsection_table.empty:
        return None

    #TODO: check if column names exists in cfg.xsection_layer_map.empty
    #TODO: deal with NaN entries, i.e. no (layer, datatype) match.

    # rename to internal column names and create layer-xsection merge.
    layers = cfg.layer_table[layer_table_attr_csv].rename(
        columns=dict(zip(
            layer_table_attr_csv,
            layer_table_attr_join))
        )
    layers2 = layers.reset_index()
    xsections_layers = cfg.xsection_layer_map[
        xsection_layer_map_attr_csv].rename(
            columns=dict(zip(
                xsection_layer_map_attr_csv,
                xsection_layer_map_attr_join))
            )


    Merged = pd.merge(
        left=xsections_layers,
        right=layers2,
        how='left',
        on=['layer_name']
    )

    # check if all layers in the xsection_layer_map are defined
    isna = Merged['layer'].isna()
    if isna.any():
        unknown_layers = Merged['layer_name'][isna].values
        raise Exception(f"Refering to undefined layer(s) in xsection_layer_map: {unknown_layers}. "
            "Correct the spelling or add the layer(s) explicitly.")

    # reset the xsextion's layer-table when merge is called for a fresh build.
    cfg.xs_list = cfg.xsection_table['xsection'].unique()
    for xs in cfg.xs_list:
        rows = Merged['xsection'] == xs
        xs_layer_table = Merged[rows][mask_layers_attr]
        add_xsection(xs).mask_layers = xs_layer_table.set_index('layer_name')
    return None


def load_masklayers(layer_file=None, xsection_layer_file=None):
    """Load layer and xsection files and merge them.

    This function combines
    1. load_layers()
    2. load_xsection_layer_map()
    3. merge()

    Args:
        layer_file (str): layer file name of csv file
        xsection_layer_file (str): xsection file name of csv file

    Returns:
        dict: {xsection name: DataFrame with mask layer}
    """
    load_layers(filename=layer_file)
    load_xsection_layer_map(filename=xsection_layer_file)
    return merge_xsection_layers_with_layers()


def get_xsection(name):
    """Return the Xsection object corresponding to <name>.

    Args:
        name (str | Node): xsection name or a pin object

    Returns:
        Xsection: Xsection object with name <name>
    """
    if hasattr(name, 'xs'): # if name is a pin (Node) get its xsection name
        name = name.xs
    if not name in cfg.XSdict.keys():
        if name == cfg.default_xs_name:
            #add_xsection(cfg.default_xs_name)
            nd.handle_missing_xs(name)
        else:
            msg = "\nNo xsection object existing under xsection name '{}'.".format(name)
            msg += ' Available xsections are:'
            for xname in cfg.XSdict.keys():
                msg += "\n  '{}'".format(xname)
            msg += "\nAlternatively, add a new xsection as:"
            msg += "\n  add_xsection(name='{}')".format(name)
            raise Exception(msg)
    return cfg.XSdict[name]


def get_layer(layer, aslist=False):
    """Get the layer ID for a specific layer reference.

    If the <layername> does not exist a default layer ID is returned.

    Args:
        layer (str | tuple | int): layer reference by name | (layer, datatype) | layer
        aslist (bool): return layer IDs as a list, which may include multiple entries.
            default=False.

    Returns:
        str: layer_ID (layer_name) (or a list of matching layer_names is aslist=True)
    """
    if isinstance(layer, str) or layer is None:
        if not layer in cfg.layername2LDT.keys():
            layer = add_layer(name=layer, layer=None, unknown=True)
        if not aslist:
            return layer
        else:
            return [layer]

    # parse LDT:
    elif isinstance(layer, tuple):
        n = len(layer)
        if n == 3:
            LDT = layer
        elif n == 2:
            LDT = (layer[0], layer[1], cfg.default_tech)
        elif n ==1:
            LDT = (layer[0], 0, cfg.default_tech)
        else:
            mes = "Error: invalid tuple format for layer {}. Tuple length must be <= 3.".\
                format(layer)
            logger.exception(mes)
            raise Exception(mes)

    elif isinstance(layer, int):
        LDT = (layer, 0, cfg.default_tech)

    status = cfg.layername2LDT_cnt.get(LDT, 0)
    if status == 1:  # unique mapping possible
        if not aslist:
            return cfg.LDT2layername[LDT]
        else:
            return [cfg.LDT2layername[LDT]]
    if status > 1:
        # if LDT not in cfg.LDT2layername.keys():
        names = [name for name, ldt in cfg.layername2LDT.items() if ldt == LDT]
        if not aslist:
            cfg.LDT2layername[LDT] = names[0]
            if cfg.gdsload != cfg.gdsloadstore:  # avoid repeating messages during single gds load
                msg = "Non-unique layer identifier used: {0}. "\
                    "Matching layer_names for {0} are {1}. "\
                    "Continuing here with layer_name = '{2}'.".\
                    format(layer, names, cfg.LDT2layername[LDT])
                if isinstance(cfg.gdsload, str):
                    cfg.gdsloadstore = cfg.gdsload
                    msg +=  f" It occured when loading gds file '{cfg.gdsload}'. "\
                        f"This message can be removed by mapping the layer explicitly "\
                        f"in the load using the layermap keyword, e.g. layermap={{{layer}: '{cfg.LDT2layername[LDT]}'}}"
                nd.main_logger(msg, "warning")
                # Note that a non-unique (L, D) becomes tricky when importing GDS;
                # The GDS format only uses (L, D) as the layer identifier,
                # and this has to be mapped to a unique layer_name.
                # The designer has to assign that in a layermapping:
                # (L, D) -> layer_name (str) to resolve conflicts.
                # As a default the first (L, D) defined can be used.
                # Issue this warning only once in stdout, or move to a log-file.
            return cfg.LDT2layername[LDT]
        else:
            return names

    # if we get here, then the layer has not been not found yet:
    layer = add_layer(name=None, layer=LDT, unknown=True)
    if not aslist:
        return layer
    else:
        return [layer]


layer_tuple = namedtuple('layer_info',
    ['layer',
     'datatype',
     'tech'
    ]
)

def get_layer_tuple(layer):
    """Get layer information as tuple components (layer, datatype, technology).

    Returns:
        layer_tuple: layer info (L, D, T) as a named tuple
    """
    lay = get_layer(layer)
    return layer_tuple (
        layer=cfg.layername2LDT[lay][0],
        datatype=cfg.layername2LDT[lay][1],
        tech=cfg.layername2LDT[lay][2]
    )


def clear_layers():
    """Clear all objects with layer information.

    Returns:
        None
    """
    global unknown_layers
    unknown_layers = set()
    cfg.layername2LDT = {}
    cfg.layername2LDT_cnt = defaultdict(int)
    cfg.LDT2layername = {}

    cfg.layer_table.drop(cfg.layer_table.index, inplace=True)
    cfg.xsection_layer_map.drop(cfg.xsection_layer_map.index, inplace=True)
    cfg.colors.drop(cfg.colors.index, inplace=True)
    for xs in cfg.xs_list:
        get_xsection(xs).mask_layers = pd.DataFrame()
    cfg.xs_list = []


def clear_xsection_layers():
    """Drop all rows from layers and xs tables.

    Returns:
        None
    """
    cfg.xsection_layer_map.drop(cfg.xsection_layer_map.index, inplace=True)
    for xs in cfg.xs_list:
        get_xsection(xs).mask_layers = pd.DataFrame()
    cfg.xs_list = []


def clear_xsections():
    """Clear xsection dict"""
    cfg.xs_list = []
    cfg.XSdict = {}


def clear_all():
    """Clear Nazca layer and xsection structures.

    For now only call certain funtions.

    Returns:
        None
    """
    clear_xsections()
    nd.cfg.stubmap = {}
    clear_layers()


def empty_xsections():
    """Delete info on all xs. Keep layer info.

    Returns:
        None
    """
    cfg.xsection_layer_map = pd.DataFrame()
    merge_xsection_layers_with_layers()


def load_parameters(filename):
    """Obsolete"""
    cfg.dfp = pd.read_csv(filename, delimiter=',')


def get_parameter(name):
    """Obsolete.

    Set parameter value for specific name.
    """
    return float(cfg.dfp[cfg.dfp['name'] == name]['value'])


#==============================================================================
# plt
#==============================================================================
def set_plt_properties(figsize=14, cmap=None, N=32, alpha=0.3):
    """Set the default colormap to use in Matplotlib output for mask layout.

    Args:
        figsize (float): size of the matplotlib figure (default = 14)
        cmap (str): name of the colormap, (default = None uses cfg.plt_cmap_name)
        N (int): the number of colors in the map in case of a 'linearSegmentedColormap' default = 32)
        alpha (float): transparency of the cmap to see through layers (default = 0.3).

    Returns:
        None
    """
    if cmap is None:
        cmap = cfg.plt_cmap_name

    mp = mpl.cm.get_cmap(cmap)
    if isinstance(mp, mpl.colors.LinearSegmentedColormap):
        cfg.plt_cmap = []
        for n in range(N):
            cfg.plt_cmap.append(mp(n/N))
    else:
        cfg.plt_cmap = mp.colors

    cfg.plt_cmap_size = len(cfg.plt_cmap)
    cfg.plt_figsize = figsize
    cfg.plt_alpha = alpha

#make sure colors are set:
set_plt_properties()


#==============================================================================
# colors
#==============================================================================
color_defaults = {
    'fill_color': 0,
    'frame_color': 0,
    'transparent': True,    # Klayout transparancy flag
    'frame_brightness':32,  # Klayout transparency level
    'fill_brightness': 32,  # Klayout transparency level
    'dither_pattern': 'I0', # KLayout dither
    'marked': False,
    'animation': 1,         # blinking if True
    'valid': True,          # KLayout valid flag (editable)
    'visible': True,        # layer visible if True
    'width': 0,             # frame width
    'alpha' : 0.3           # Matplotlib transparency level in
}


def load_layercolors(filename):
    """Read colormap from a .csv Nazca color table file.

    Args:
        filename (str): filename

    Returns:
        None
    """
    block_load_colors = True
    # For Klayout tabs with groups the layer column contains a '*'
    # and Pandas will create type->str
    # For tabs without groups the layer is read as type->int
    df1 = pd.read_csv(filename)

    # filter out expected columns
    df1 = df1[colors_attr_csv]

    df1['layer'] = df1['layer'].astype(str)
    #df1.loc[df1['width'] == 'na', 'width'] = 0

    df1.width.replace('na', 0, inplace=True)
    #df1.replace('na', '', inplace=True)

    #TODO: creates a brand new table and forgets present settings
    #  needs an update to handle mutliple layer/color sets.
    cfg.colors = df1[df1['layer'] != '*'].dropna(subset=['layer'])
    cfg.colors.reset_index(inplace=True)
    try:
        cfg.colors['layer'] = cfg.colors['layer'].astype(int)
        cfg.colors['datatype'] = cfg.colors['datatype'].astype(int)
    except:
        logger.exception('Check load_layercolors.')
        #print('Group in load_layercolors:', cfg.colors['layer'], '\n')
    #TODO: the try will trigger on any layer entry with text.
    # as a resuls all the layers and datatype remain a string.

    # Using a color_name key rather than layername, because mask viewer use the
    #  (L, D) and in Nazca these have to map to the all layer_names
    #  having these (L, D). Hence the color_name based on "L/D"
    cfg.colors['color_name'] = cfg.colors.apply(
         lambda row: "{}/{}".format(row.layer, row.datatype), axis=1)
    cfg.colors['alpha'] = 0.3
    return None


def set_layercolor(
    layer=None,
    frame_color=None,
    fill_color=None,
    frame_brightness=None,
    fill_brightness=None,
    dither_pattern=None,
    valid=None,
    visible=None,
    transparent=None,
    width=None,
    marked=None,
    animation=None,
    alpha=None):
    """Set layer color information.

    Generate a tabel (DataFrame) with layer color information.
    For missing attribute values for fill_color and frame_color,
    the default_colors as specified in colormap cfg.plt_cmap are be applied.
    Note that the default colormap can be adjusted by the user.

    Args:
        layer (int | tuple | str): mask layer
        tech (str) : technology (default = None)
        ...: attributes

    Returns:
        None
    """
    if layer is None:
        raise ValueError('Need to provide a layer number or name to set layer colors.')
    layer_name = get_layer(layer)
    layer, datatype, tech = cfg.layername2LDT[layer_name]
    color_name = "{}/{}".format(layer, datatype)
    kltab_name = "{}/{} {}".format(layer, datatype, layer_name)

    colors = {
        'depth': 1, # group depth
        'layer_name': layer_name,
        'name': kltab_name,
        'color_name': color_name,
        'layer': layer,
        'datatype': datatype,
        'fill_color': fill_color,
        'frame_color': frame_color,
        'frame_brightness': frame_brightness,
        'fill_brightness': fill_brightness,
        'dither_pattern': dither_pattern,
        'valid': valid,
        'visible': visible,
        'transparent': transparent,
        'width': width,
        'marked': marked,
        'animation': animation,
        'alpha' : alpha
    }

    assert isinstance(layer, int)
    assert isinstance(datatype, int)
    if not cfg.colors.empty:
        mask = cfg.colors['color_name'] == color_name
        color_duplicate = cfg.colors.loc[mask]
        coloritems = len(color_duplicate)
    else:
        coloritems = 0

    df = None
    if coloritems == 0:
        #TODO: not here or the colors can't be changed later on.
        #col = mpl.colors.to_hex(cfg.plt_cmap[layer % cfg.plt_cmap_size])
        col = rgb2hex(cfg.plt_cmap[int(layer % cfg.plt_cmap_size)])
        if fill_color is None:
            colors['fill_color'] = col
        if frame_color is None:
            colors['frame_color'] = col
        for attr in colors.keys():
            if colors[attr] is None:
                if attr in color_defaults.keys():
                    colors[attr] = color_defaults[attr]

        df = pd.DataFrame(colors, index=[0])

        # note the df will get a dtype numpy int, not native Python int
        try: #pandas in python 3.7
            cfg.colors = pd.concat([cfg.colors, df], ignore_index=True, sort=True)
        except:  #pandas in python <= 3.6
            cfg.colors = pd.concat([cfg.colors, df], ignore_index=True)

    elif coloritems == 1:
        for attr in cfg.colors.columns: # use existing values for undefined attributes
            if attr in colors.keys():
                if colors[attr] is None:
                    #TODO: Better use a try to catch color_defaults are missing?:
                    if attr in color_defaults.keys():
                        colors[attr] = color_duplicate[attr].iloc[0]
            else:
                colors[attr] = color_duplicate[attr].iloc[0]
        df = pd.DataFrame(colors, index=color_duplicate.index) #.set_index(idx)

        # TODO: this could be a new (extra) layername with an old layer number
        # now it's always overwritten when add_layer() is called:
        cfg.colors.loc[mask] = df
    else:
        logger.warning('Multiple color entries for ML({}, {}). Skipping...'.format(layer, datatype))

    return df


def save_colormap(filename=None):
    """Save a colormap to csv file.

    Args:
        filename (str): output filename

    Returns:
        None
    """
    if filename is None:
        filename = "colors_saved.csv"
    cfg.colors.reset_index()
    bools = ['valid', 'visible', 'transparent', 'marked']
    colors = cfg.colors.copy()
    for colm in bools:
        mask = colors.loc[:, colm] == True
        colors.loc[mask, colm] = 'true'
        colors.loc[~mask, colm] = 'false'
    colors[colors_attr_csv].to_csv(path_or_buf=filename, sep=',')
    return None


#==============================================================================
# Show or get tables
#==============================================================================
def get_layers():
    """Get predefined selection of columns from layer_table DataFrame.

    Returns:
        DataFrame: layer info ['layer_name', 'layer', 'datatype', 'tech', 'accuracy']
    """
    df = cfg.layer_table[['layer', 'datatype', 'tech', 'accuracy']]
    return df


def show_xsections():
    """Print the xsections table.

    Returns:
        None
    """
    #print('-------------------------------')
    print('xsections:')
    if not cfg.xsection_table.empty:
        pprint(cfg.xsection_table[xsection_table_attr])
    else:
        print('No xsections defined.')


def show_layers():
    """Print the layer table."""
    #print('-------------------------------')
    #print("DataFrame 'cfg.layer_table' with layer definitions:")
    print("layers:")
    try:
        #df = cfg.layer_table[['layer_name', 'layer', 'datatype', 'tech', 'accuracy']]
        df = cfg.layer_table
        pprint(df)
    except:
        print('  ...no layers defined.')


def show_layercolors():
    """Print the layercolor table."""
    #print('-------------------------------')
    #print("DataFrame 'cfg.layer_table' with layer definitions:")
    print("layercolors:")
    df = cfg.colors[['fill_color', 'frame_color', 'width']]
    pprint(df)


def show_xsection_layer_map():
    """Print the xsection_layer_map table to stdout.

    Returns:
        None
    """
    #print('-------------------------------')
    print("xsection-layer_map:")
    if not cfg.xsection_layer_map.empty:
        df = cfg.xsection_layer_map[['xsection', 'layer_name', 'growx', 'growy']]
        pprint(df)
        #return df
    else:
        print('No xs defined.')


def show_mask_layers():
    """Print mask_layers dictionary with layer export definitions to stdout.

    Returns:
        None
    """
    #print('-------------------------------')
    #print("Dictionary of DataFrames 'cfg.xs_layers' with the description of each xs:")
    columns = ['layer_name', 'layer', 'datatype', 'tech', 'growx', 'growy', 'accuracy']
    print("mask_layers:")
    for xs in cfg.xs_list:
        #for sorted(cfg.mask_layers.keys()):
        under = '_'*(29-len(xs))
        print("xsection: '{}' {}".format(xs, under))
        pprint(get_xsection(xs).mask_layers.reset_index()[columns])
        print('')




