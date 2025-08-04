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

# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
#==============================================================================
# (c) 2016-2022 Ronald Broeke, Katarzyna Lawniczuk, Xaveer Leijtens
#==============================================================================

"""Module with a set of Buiding Block templates for faciliating PDK creation."""

import os
from collections import OrderedDict
import inspect
from functools import wraps
from ast import literal_eval
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import yaml
from pprint import pprint
from PIL import Image

import hashlib
import nazca as nd
from nazca import cfg
from nazca.logging import logger
from nazca.clipper import merge_polygons

import nazca.geometries as geom



#==============================================================================
# hash
#==============================================================================
_hash_id = None
_hash_params = {}
_hash_name = 'NONAME'
cfg._hashme = []
HASH_LENGTH = 4  # number of hash characters to add to the cellname


def hashme(name="", parreprs=[], suffix=None, addhash=None, config=None):
    """Decorator to reuse building blocks (cells) based on the same parameters.

    The @hashme decorator can be used for parametrized cells.
    It checks if a function call (that returns a Cell) occured before
    with the same parameters. If that is the case @hashme returns a reference
    to the previously created cell with these parameters, rather than creating
    a full copy with a new cell name. The hashme ensures a unique cellname by
    adding an md5 hash at the end of the cellname. Note that in gds cell the cellname
    string is the unique identifier for referencing cells. In order to add
    parameters to the cellname each specific parameter name must be added to "parrepr".

    Args:
        name (str): Cell name.
        parreprs (list of dict | list of str): list of parrepr 'parameter representations'
            in the cellname. In case of proving a list-of-str, only the parameter names
            shall be provided, which will automatically be translated in a dict
            internally. A parameter entry in parreprs adds the parameter
            to the cellname in order of appearence in parreprs.
            Each parrepr dict describes not only if, but also how to format
            cell parameter in the cellname. See parrepr formatting rules below.
        config (dict): if set, overrules the default setting in cfg.hashme_config
        suffix (str): cellname suffix. Overules 'config'.
        addhash (bool): flag to show a hash in the cellname.
            Overules 'addhash' setting in 'config'.
            There maybe foundry rules that ban the use of hashme. Note that in such
            cases the parreprs need to be used to uniquely define cell names.

    The default settings for the hashme cellname formatting are set in the
    cfg.hashme_config dict. These setting can be overruled with the "config" dict
    in the hashme call, which has the same structure. Using the "config" keyword
    also facilitate reuse of settings, and simplifies the "parreprs" entry.
    Settings in parrepr overrule settings in "config".
    Note that there can be multiple parameters per cellname, but only one "addhash"
    and one "suffix". The config and parrepr dicts are defined as follows:

        config = {
           "prefix": <parameter specific prefix>,  # Generic prefix for the parameter value, e.g. "", "_", "_L", "_w"
           "format": <format string>,  # parameter format, e.g. ".5f", "2.1g"
           "func": <function(value)>,  # custom formatting function for callback, returns a string
           "suffix": <suffix string>,  # can be overruled by "suffix" keyword in hashme.
           "addhash": <addhash True|False>  # can be overruled by "addhash" keyword in hashme.
        }

        parrepr = {
            "parameter": <parameter name>,  # mandatory parameter name as used in the cell call, e.g. "L", "width"
            "prefix": <parameter specific prefix>,  # optional prefix for the parameter value, e.g. "", "_", "_L", "_w"
            "format": <format string>,  # optional parameter format, e.g. ".5f", "2.1g"
            "func": <function(value)>,  # optional custom formatting function for callback, returns a string
        }

    Example:

        Usage of the hashme decorator::

            @hashme('mycellname')
            def parametrized_cell(W, L):
                with nd.Cell() as C:
                    # cell stuff
                    nd.strt(length=L, width=W).put()
                return C

            parametrized_cell(W=3.0, L=10.0)
            # cell_name ='mycellname'

        Usage of the hashme with adding keyword values to the cellname (default representation)::

            @hashme('mycellname', ['L', 'W'])
            def parametrized_cell(W, L):
                with nd.Cell() as C:
                    # cell stuff
                    nd.strt(length=L, width=W).put()
                return C

            parametrized_cell(W=3.0, L=10.0)
            # cell_name ='mycellname_10.0_3.0'

        Usage of hashme with adding keyword values to the cellname (custom representation)::

            parreprs = [
                {"parameter": "L", "prefix": "_L"},  # Default float format: '.5g'
                {"parameter": "W", "format": ".5g", "min_decimals": 1},  # Force a minimum number of decimals
            ]
            @hashme('mycellname', addhash=False, parrepr=parreprs, suffix="_end")
            def parametrized_cell(W, L):
                with nd.Cell() as C:
                    # cell stuff
                    nd.strt(length=L, width=W).put()
                return C

            parametrized_cell(W=3.0, L=10.0)
            # cell_name = 'mycellname_L10_3.0_end'

        Usage of hashme with config::

            config = {
                "prefix": "/",
                "format": ".2g,"
                "addhash: False",
                "suffix": "_the_end",
            }
            @hashme(
                'name=mycellname',
                 parrepr=[
                     {"parameter": "L", "prefix": "_L"},
                     {"parameter": "W", "prefix": "_W"},
                 ],
                 config=config,
            )
            def parametrized_cell(W, L):
                with nd.Cell() as C:
                    # cell stuff
                    nd.strt(length=L, width=W).put()
                return C

            parametrized_cell(W=3.0, L=10.0)
            # cell_name = 'mycellname_L10_3.0_end'


    Note: do *not* use any global parameters in a parameterized Cell function
    that is decorated with hashme, i.e. put all parameters that may change the
    cell in the function call. The reason is that a change of a global parameter
    could possibly cause a different result of the generated cell without
    being noticed by the function introspection. For this reason function
    in classes should not be decorated with @hashme.

    Returns:
        decorator
    """

    hash_class_warning =\
"""Warning: @hashme decorator on Class method:

By default Nazca will guarantee a unique cell name each time a cell is created.
It does so by suffixing cell names with an ordinal counter.
This may result in multiple copies of the same cell where only one would do.

Decorating the cell generating function with @hashme avoids cell copies of
identical cells by returning an existing cell if the function has been
called with the same parameters earlier.

The @hashme 'state checker' can only be used safely if the state of the
function it decorates is *not* dependent on a state *external* to that
function, i.e. only if the variations of the cell it generates *only* depends
on the function's parameters. In contrast, Class methods typically depend
on attributes stored at class level, hence, by their very nature Class methods
do not fit well with the @hashme concept to guarantee unique cells for
unique names.

Example of a cell function where @hashme should NOT be used, i.e. in a class method:

    class myElement():
        # stuff
        def make_cell(a, b, c):
            # stuff
            return cell_object

    elm = myElement()

Now compare two ways to place a cell multiple times.

* method 1:
    elm.make_cell(a, b, c).put()
    elm.make_cell(a, b, c).put()

a) creates two copies of the cell and suffixes the name with an ordinal counter.
b) if @hashme is used it returns an existing cell if the call profile
   of make_cell has been used earlier.

* method 2:
A secure way in all cases to reuse cells is to call a cell generation function
only once and assign its returns value, i.e. the Cell object, to a variable.
Use this variable to put instances of the cells:

    E = elm.cell(a, b, c)
    E.put()
    E.put()
"""
    #print('params:' , params)
    if config is None:
        config = {}
    if name == "":
        name = _hash_name
        print(f"warning: no cellname given to @hashme, calling it '{name}'.")

    if len(parreprs) > 0:
        if isinstance(parreprs, str):  # single string
            parreprs = [{"parameter": parreprs}]
        elif isinstance(parreprs, list) and isinstance(parreprs[0], str):  # list of str
            parreprs = [{"parameter": par} for par in parreprs]
            parreprs = parreprs
        # TODO, check beyond item [0] may be needed:
        elif isinstance(parreprs, list) and isinstance(parreprs[0], dict):  # list of dict, desired format
            pass
        else:
            raise Exception(
                "Hashme input variable 'parreprs' (parameter representations) should be a list-of-str or list-of-dict. "\
                f"Current value: {parreprs}."
            )

        # Fill all values for parameter representation in the cellname:
        for parrepr in parreprs:
            for key in ["prefix", "format", "func"]:
                if key not in parrepr.keys():
                     parrepr[key] = config.get(key, cfg.hashme_config[key])

    if addhash is None:
        addhash = config.get("addhash", cfg.hashme_config['addhash'])
    if suffix is None:
        suffix = config.get("suffix", cfg.hashme_config['suffix'])
    assert isinstance(addhash, bool), \
        f"Hashme parameter 'addhash' should be a boolean. Current value: {addhash}."

    def decorator(cellfunc):
        """Decorator for grabbing and hashing a function's full parameter list.

        While executing the wrapper a global hashid and parameter list are set.
        Before leaving the wrapper the hashid is reset to None.

        Returns:
            function: wrapper function
        """
        @wraps(cellfunc)
        def wrapper(*args, **kwargs):
            global _hash_id
            global _hash_params
            global _hash_name
            nonlocal name

            if cfg.hashme:
                raise Exception("Can not call an @hashme-decorated function before the Cell() call "\
                    "inside another @hashme decorated function. Move the function call(s) after the 'with Cell()' context. "\
                    f"Cellname: '{name}'"
                )
            name_long = name
            getargs = inspect.getfullargspec(cellfunc)
            funcargs = getargs.args
            funcdefaults = getargs.defaults
            if funcdefaults is not None:
                delta = len(funcargs) - len(funcdefaults)
            else:
                delta = len(funcargs)
            hashstr = ''
            _hash_params = OrderedDict()

            allow_kwargs = getargs.varkw is not None  # Check if kwargs can be passed to the function
            unknown_params = list(set(kwargs.keys()) - set(funcargs))  # Check if non-defined function arguments are passed
            if unknown_params and not allow_kwargs:
                raise TypeError(f"{cellfunc.__name__}() got unexpected keyword argument(s) {unknown_params}")

            if funcargs:
                if funcargs[0] == 'self':  # Assuming consistent use of string 'self' for class refs.
                    funcargs = funcargs[1:]
                    print(hash_class_warning)

                for i, a in enumerate(funcargs):
                    if i < len(args):
                        default = args[i]
                    elif i - delta >= 0:
                        default = funcdefaults[i - delta]
                    else:
                        default = None
                    _hash_params[a] = kwargs.get(a, default)

                for parrepr in parreprs:
                    val = parrepr["parameter"]
                    parfunc = parrepr["func"]
                    if val in _hash_params.keys():
                        pvalue = _hash_params[val]
                        try:
                            pstr = format(pvalue, parrepr["format"])
                        except:
                            pstr = pvalue
                        if parfunc is not None:
                            pstr = parfunc(pstr) # call formatting function
                        name_long += f"{parrepr['prefix']}{pstr}"
                    else:
                        raise ValueError(
                            f"Parameter '{val}' not found in function parameters {list(_hash_params.keys())}."
                        )

                for p in sorted(_hash_params.keys()):
                    hashstr += '{}_{}'.format(p, _hash_params[p])
                _hash_id = hashlib.md5(hashstr.encode()).hexdigest()[:HASH_LENGTH]

            else:
                _hash_id = ''

            if suffix != '':
                name_long = f"{name_long}{suffix}"

            if (_hash_id != '') and addhash:
                _hash_name = f"{name_long}_${_hash_id}"
            else:
                _hash_name = name_long
            _hash_name = cfg.gds_cellname_cleanup(_hash_name)
            # add cfg to avoid explicit import of this module for hashme attributes
            cfg.hashme = True
            cfg.hash_id = _hash_id
            cfg.hash_cellnameparams = tuple(parrepr['parameter'] for parrepr in parreprs)
            cfg.hash_params = _hash_params
            basename = cfg.gds_cellname_cleanup(name)
            cfg.hash_basename = basename                               # base
            cfg.hash_paramsname = cfg.gds_cellname_cleanup(name_long)  # base + params
            cfg.hash_name = _hash_name    # base + params + hash
            cfg.hash_func_id = id(cellfunc)

            # check if basename comes from a unique function:
            funcid = cfg.basenames.get(basename, False)
            if funcid:
                if funcid != cfg.hash_func_id and cfg.check_basename:
                    logger.warning(f"Reusing a basename across function IDs: '{basename}'")
            else:
                logger.debug(f"Adding basename '{basename}' to cfg.basenames")
                if cfg.validate_basename:
                    nd.validate_basename(basename)
                # add basename:
                cfg.basenames[basename] = cfg.hash_func_id

            cell_exists = _hash_name in cfg.cellnames.keys()

            if cell_exists:
                cfg.hashme = False
                return cfg.cellnames[_hash_name]
            else:
                cell_new = cellfunc(*args, **kwargs)

            #reset:
            _hash_id = None
            _hash_params = {}
            _hash_name = ''
            cfg.hashme = False
            cfg.hash_id  = ''
            cfg.hash_params = ''
            cfg.hash_cellnameparams = {}
            cfg.hash_basename = ''
            cfg.hash_paramsname = ''
            cfg.hash_name = ''
            cfg.hash_func_id = None
            return cell_new
        return wrapper
    return decorator


def trial_cell(func):
    """Decorator to delete all Nazca cells generated in the decorated functions.

    It will reset the cell counter to before the trial

    Returns:
        callable: wrapped function
    """
    def wrapper(*args, **kwargs):
        """Wrapper function that deletes new Cells generated in its scope.

        Create a trial cell to extract cell properties .

        Returns:
            Cell | float: result of wrapped function.
        """
        cell = kwargs.pop('cell', False)
        if not cell:
            S1 = set(nd.cfg.cellnames.keys())
            num = nd.netlist.num
        out = func(*args, cell=cell, **kwargs)
        if not cell:
            nd.netlist.num = num
            S2 = set(nd.cfg.cellnames.keys())
            for s in S2 - S1:
                nd.cfg.cellnames.pop(s)
        return out
    return wrapper


# create layers_ignore list to exclude layer from being part of bbox calculations:
cfg.bbox_layers_ignore = []
def bbox_layers_ignore(layers=None, reset=False):
    """Set the layers that should be ignored in bbox calculations.

    Args:
        layers (layer): list of layers
        reset ((bool): If True, clear existing ignore list and only set <layers>, default=False.

    Returns:
        None
    """
    if layers is None:
        layers = []
    layerIDs = [nd.get_layer(ml) for ml in layers]
    if reset:
        cfg.bbox_layers_ignore = layerIDs
    else:
        cfg.bbox_layers_ignore.extend(layerIDs)


def get_Cell(name):
    """Get Cell object by providing cell name.

    Args:
        name (str): cell name

    Returns:
        Cell: cell having cellname <name>
    """
    if name in cfg.cellnames.keys():
        return cfg.cellnames[name]
    else:
        mes = f"Error: Requested cell with name '{name}' not existing. "\
            f"Availabe are the following {len(cfg.cellnames)} cells: {sorted(cfg.cellnames.keys())}."
        logger.exception(mes)
        raise Exception(mes)


def rangecheck(allowed_values):
    """Check if parameter values are in range.

    A dictionary with the allowed values for parameters is sent to rangecheck.
    If any parameter is out of range, a ValueError will be raised and relevant
    error info is provided to solve the issue.

    Args:
        allowed_values (dict): {'<var_name>': (low, <var_name>, high)},
        where the key is the <varname> in a str type,
        and low and high are the min and max values

    Raises:
        ValueError if out of range.

    Returns:
        None

    Example:

        Add a check in a function on 0 <= a <= 10 and -5 <= b <= 5.
        The call will raise a ValueError::

            def func(a, b):
                nd.rangecheck({'a': (0, a, 10), 'b': (-5, b, 5)})

            func(a=100, b=100)

            # output:
            # ValueError: a=100 too large. Allowed values: 0<=a<=10
            # b=100 too large. Allowed values: -5<=b<=5
    """
    err = ''
    for varname, values in allowed_values.items():
        if isinstance(values, tuple): # for backward compatibility
            #logger.warning("obsolete, use ...")
            if len(values) != 3:
                raise ValueError(f"Need to provide 3 values, got {len(values)}")
            else:
                values= {
                    'min': values[0],
                    'default': values[1],
                    'max': values[2],
                    'type': None,
                    'doc': '',
                    'unit': '',
                }
                if values['default'] is None:
                    raise Exception(f"Value None is not allowed in rangecheck() in cell '{nd.cfg.cells[-1].cell_name}'")
        elif not isinstance(values, dict):
            raise Exception (f"Expected type dict, not type '{type(values)}'.")

        if values['type'] not in cfg.allowed_vartypes:
            err += f"vartype '{values['type']}' unknown. Allowed vartypes: '{cfg.allowed_vartypes}'"
        if values['min'] is None:
            values['min'] = float('-inf')
        if values['max'] is None:
            values['max'] = float('inf')
        #if low <= value <= high:
        #    continue
        if values['default'] < values['min']:
            err += f"{varname}={values['default']} too small. Allowed values: {values['min']}<={varname}<={values['max']}"
        elif values['default'] > values['max']:
            err += f"{varname}={values['default']} too large. Allowed values: {values['min']}<={varname}<={values['max']}"

        # store variable info in cell for later use:
        cell = nd.cfg.cells[-1]
        cell.ranges[varname] = values

    if err != '':
        if cfg.rangecheck:
            raise ValueError(err)
        else:
            nd.main_logger(f"{err}, continuing anyway (rangecheck is set to False).", "error")

    return None


#==============================================================================
#
#==============================================================================
stubs = dict() # dict of all stubs. {stub_name: stub_cell_obj}
stubmap = dict() # dict of xsection_name to its stub's xsection_name. {xs_name: stub_xs_name}

class Functional_group():
    """Class to group building blocks syntactically together.

    Example:
        Create bb_splitters group with various splitter components
        such that the mmi can be found under bb_splitters.<tab>::

            bb_splitters = Functional_group()
            bb_splitters.mmi1x2 = mmi1x2
            bb_splitters.mmi2x2 = mmi2x2
            bb_splitters.mmi3x3 = mmi3x3
    """

    def __init__(self, name=''):
        """Create a Functional group object."""
        self.name = name
BB_Group = Functional_group


# Put in this module to avoid loading utils.py:
def parameters_to_string(param):
    """Create a string from a parameter dictionary.

    Args:
        param (dict): (parameter_name, value)

    Note that the output format, i.e. inside the annotation object)
    may be foundry specific depending on the setting of
    cfg.parameter_style

    'default' output format:

    Parameters:
    <parameter> = <value>
    <parameter> = <value>
    ...

    'yaml' output format:

    <parameter>: <value>
    <parameter>: <value>
    ...

    Returns:
        str: parameters as a string
    """
    if cfg.parameter_style == 'yaml':
        return yaml.dump(dict(param)) #  change from ordered dict
    elif cfg.parameter_style == 'hhi':
        plist = []
        for key, value in param.items():
            plist.append("{}:{}".format(key, value))
        return '\n'.join(plist)
    elif  cfg.parameter_style == 'default':
        plist = ['Parameters:']
        for key, value in param.items():
            plist.append("{} = {}".format(key, value))
        return '\n'.join(plist)
    else:
        nd.logger.error(f"cfg.parameter_style = '{cfg.parameter_style}' not recognized.")


def string_to_parameters(string):
    """Convert a string to a parameter dictionary.

    The returned parameter values are represented as type str.

    Expected format of <string>:

    "parameters:
    <parameter> = <value>
    <parameter> = <value>
    ..."

    Header 'parameters:' is case incensitive and spaces will be stripped.

    Args:
        string (str): parameters

    Returns:
        OrderedDict: {<parameter_name>: <parameter_value>}
    """
    if cfg.parameter_style == 'default':
        lines = string.split('\n')
        p = OrderedDict()
        if (lines[0].lower() == 'parameters:'):
            for line in lines[1:]:
                param = line.split('=', 1)
                if len(param) == 2:
                    p[param[0].strip()] = param[1].strip()
                else:
                    logger.error(
                        f"string_to_parameter: Expected one keyword "
                        f"and one value,\nbut found instead: {param}. "
                        f"Provided string instead: '{string}'"
                    )
        else:
            logger.exception(
                f"string_to_parameter: "
                f"Expected header 'parameters:', because cfg.parameter_style = '{cfg.parameter_style }'. "
                f"Provided string instead: '{string}'"
            )
        return p

    elif cfg.parameter_style == 'yaml':
        p = yaml.load(string, Loader=yaml.Loader)
        return p

    else:
        logger.error("string_to_parameter: "
            "Expected header 'parameters:'\n"
        )
        return {}


def put_parameters(parameters=None, pin=None):
    """Put a parameter list as annotation in a building block.

    Args:
        parameters (dict): {parameter_name: value}

    Returns:
        Annotation: Annotation object with the parameterlist
    """
    cell = cfg.cells[-1]
    if isinstance(parameters, dict):
        cell.parameters = parameters
        text = parameters_to_string(cell.parameters)
    elif cell.hashme:
        text = parameters_to_string(cell.parameters)
    else:
        text = ''

    if cfg.parameter_style == 'hhi' :
        for i, plist in enumerate(text.split('\n')):
            anno = nd.Annotation(text=plist, layer='bb_parameter_text')
            if pin is None:
                anno.put(0, 10*i)
            else:
                anno.put(pin.move(0, 5*i))
    else: # standard
       anno = nd.Annotation(text=text, layer='bb_parameter_text')
       if pin is None:
           anno.put(0)
       else:
           anno.put(pin)
    return text


def cellname(cellname=None, length=0, width=None, align='lc', name_layer=None):
    """Create the cellname as a text cell.

    Args:
        cellname (str): name of the cell
        length (float): length available for the BB name in um
        width (float):
        align (str): text alignment (see nazca.text, default = 'lc')
        name_layer (str): name of the layer to be used for bbox name

    Returns:
        Cell: text with cellname
    """
    cell = cfg.cells[-1]
    if cellname is None:
       cellname = cell.cell_paramsname

    if width is not None:
        maxheight = min(width*0.8, cfg.cellname_max_height)
    else:
        maxheight = cfg.cellname_max_height

    if name_layer is None:
        name_layer = 'bbox_name'

    length = length * cfg.cellname_scaling
    texth = min(maxheight, abs(length) / nd.linelength(cellname, 1))
    txt = nd.text(cellname, texth, layer=name_layer, align=align)
    return txt


#==============================================================================
# stubs and pins
#==============================================================================
def add_pinstyle(name, styledict):
    """Set a custom pin style for a technology.

    Start with the style in <name> if existing, otherwise with the
    default cfg.pinstyle settings and overrule these
    for the keywords provided in the <styledict>.

    Each pinstyle name will point to a unique dictionary internally.

    Args:
        name (str): pinstyle name
        styledict (dict): pinstyle

    Returns:
        None
    """
    present = nd.cfg.pinstyles.get(name, None)
    if present is None:
        present = nd.cfg.pinstyle
    newstyle = present.copy()

    #keys = nd.cfg.pinstyle.keys()
    if styledict is not None:
        for item, setting in styledict.items():
            if item in cfg.pinstylelabels:
                newstyle[item] = setting
            else:
                nd.main_logger(
                    f"No valid pinstyle key '{item}' for pinstyle. Available keys: {cfg.pinstylelabels}",
                    "warning",
                )
    nd.cfg.pinstyles[name] = newstyle


# for backward compatibility:
raise_set_pinstyle = False
set_pinstyle_flag = False
def set_pinstyle(name, styledict):
    global set_pinstyle_flag
    if not set_pinstyle_flag:
       print("Method set_pinstyle() is obsolete. Use add_pinstyle() instead."\
            " Find the source of this issue by adding nd.raise_set_pinstyle=True after import nazca." )
       if nd.raise_set_pinstyle:
           raise Exception("set_pinstyle() is obsolete.")
       add_pinstyle(name, styledict)

       set_pinstyle_flag = True
    add_pinstyle(name, styledict)


def add_pinshape(name, shape, overwrite=True):
    """Add a pin shape to the pinshapes dict.

    Args:
        name (str): reference name of the shape
        shape (list of (float, float)): pin shape [(x1, y1), ...]

    Returns:
        None
    """
    if name in cfg.pinshapes.keys() and not overwrite:
        print(f"pinshape name '{name}' already exists. Use an other name of set overwrite=True.")
    else:
        cfg.pinshapes[name] = shape


def set_default_pinstyle(stylename):
    """Provide the name to use for the standard pinstyle.

    The default standard style utilized in 'default'.

    Returns:
        None
    """
    cfg.default_pinstyle = stylename


#@hashme('arrow', 'layer', 'pinshape', 'pinshape')
# Do NOT @hashme this function,
# because parameter values are finalized inside the function.
def make_pincell(layer=None, shape=None, size=None, style=None):
    """Create a cell to indicate a pin position.

    The cell contains a shape, e.g. an arrow, to point out a location in the layout.
    Available pin shapes are in a dictionary in cfg.pinshapes: {name: polygon}.
    The predefined shapes have been normalized to unit size.

    Args:
        layer (float): layer number to place the pin symbol/shape
        shape (str): name (dict key) of the pin symbol/shape
        size (float): scaling factor of the pin symbol/shape
        style (str): set a new pinstyle, to be overruled by other keywords, e.g. size, if provided.

    Returns:
        Cell: cell with pin symbol
    """
    if style is not None:
        pinstyle = cfg.pinstyles.get(style, None)
        if pinstyle is None:
            print("Error: pinstyle '{}' unknown. Known styles: {}".\
                format(style, list(cfg.pinstyles.keys())))
            pinstyle = cfg.pinstyles['default']
    else:
        pinstyle = cfg.pinstyles[cfg.default_pinstyle]

    if shape is None:
        shape_use = pinstyle['shape']
    else:
        shape_use = shape

    if layer is None:
        layer = pinstyle['layer']
    layer = nd.get_layer(layer)

    if size is None:
        size = pinstyle['size']

    #if cfg.black_pinstyle is not None:
    #    pinstyle = cfg.pinstyles[cfg.black_pinstyle]
    #    shape_use = cfg.pinstyles[cfg.black_pinstyle]['shape']

    if shape_use not in cfg.pinshapes.keys():
        mes = "pinshape indentifier '{}' not recognized.".format(shape)
        mes += "Available options are: {}.".format(cfg.pinshapes.keys())
        mes += "Fall back type '{}' will be used.".format(shape_use)
        logger.warning(mes)
    shape = shape_use

    name = f"arrow_{layer}_{shape}_{size}"
    name = cfg.gds_cellname_cleanup(name)
    if name in cfg.cellnames.keys():
        return cfg.cellnames[name]
    with nd.Cell(name, instantiate=cfg.instantiate_pin, store_pins=False) as arrow:
        arrow.auxiliary = True
        outline = [(x*size, y*size) for x, y in cfg.pinshapes[shape]]
        nd.Polygon(layer=layer, points=outline).put(0)
    return arrow


def stubname(xs, width, thick, stubshape=None, pinshape=None, pinsize=None, pinlayer=None):
    """Construct a stub name.

    Args:
        xs (str); xsection name
        width (float): stub width
        thick (float): thickness of stub into cell (length)

    Returns:
        str: stub name
    """
    if width is None:
        return 'stub_{}'.format(xs)
    else:
        name = 'stub_{}_w{}_t{}_{}_{}_s{}_l{}'.\
            format(xs, width, thick, stubshape, pinshape, pinsize, pinlayer)
        return cfg.gds_cellname_cleanup(name)


missing_xs = []  # keep a list of encountered misssing xs to flag them only once.
def _makestub(
    xs_guide=None,
    width=0,
    length=2.0,
    shape=None,
    pinshape=None,
    pinsize=None,
    pinlayer=None,
    cell=None,
    pinstyle=None,
):
    """Create a stub cell (in the logical layers) and add it to the stub dict.

    A stub is the stub of a xsection shape around a pin to visualize a connection.
    A pincell is added to the stub to indicate the pin position inside the stub.
    The new stub is added to the stubs dictionary: {name: stubcell}

    Args:
        xs_guide (str): name of xsection
        width (float): stub width
        thick (float): thickness of stub into cell (length)
        shape (str): shape of the stub: 'box' | 'circ' (default = 'box')
        pinshape (string): pinshape used in the stub
        pinsize (float): scaling factor of the pinshape (default = 1)
        cell (Cell): use the provided cell as stub instead of creating a new stub cell
        pinstyle (dict):

    Returns:
        str: key/name of the stub
    """
    if cfg.validation_stub_xs is not None:
        xs_logic = cfg.validation_stub_xs
    else:
        xs_logic = cfg.stubmap.get(xs_guide, None)
        if xs_logic is None:  # if no stub defined, use the xs as its own stub.
            if xs_guide is None:  # no stub defined
                arrow = make_pincell(
                    layer=pinlayer,
                    shape=pinshape,
                    size=pinsize,
                    style=pinstyle,
                )
                return arrow
            xs_logic = xs_guide

        if xs_guide not in cfg.XSdict.keys():
            if xs_guide in nd.cfg.default_xs_list:
                nd.add_xsection(name=xs_guide, pinstyle=nd.cfg.default_xs_list[xs_guide]['pinstyle'])
            elif xs_guide not in missing_xs:
                missing_xs.append(xs_guide)
                if xs_guide != cfg.default_xs_name:
                    logger.error("Can not make a stub in cell '{3}' in undefined xsection '{0}'.\n"\
                       "  Possible causes: '{0}' is misspelled or not yet defined.\n"\
                       "  Will use xsection '{2}' instead and continue.\n"
                       "  To define a new xsection:\n"\
                       "      add_xsection(name='{0}')\n"\
                       "  or with layers info and adding a custom stub:\n"\
                       "      add_xsection(name='{0}', layer=1)\n"\
                       "      add_xsection(name='{1}', layer=2)\n"\
                       "      add_stub(xsection='{0}', stub='{1}')".\
                           format(xs_guide, 'stubname', cfg.default_xs_name, cfg.cells[-1].cell_name))
                xs_guide = cfg.default_xs_name

            cfg.stubmap[xs_guide] = xs_guide
            xs_logic = xs_guide

    name = stubname(xs_guide, width, length, shape, pinshape, pinsize, pinlayer)
    if name in stubs.keys():
        return name

    # make a new stub:
    stubshapes = ['box', 'circle']
    arrow = make_pincell(
        layer=pinlayer,
        shape=pinshape,
        size=pinsize,
        style=pinstyle,
    )

    if width is None:
        width = 0

    with nd.Cell(name, instantiate=cfg.instantiate_stub, store_pins=False) as C:
        C.auxilliary = True
        C.default_pins('a0', 'a0')
        nd.Pin('a0', chain=0).put(0)
        arrow.put(0)

        if cell is not None:
            cell.put()
        else:
            if length > 0:
                for lay, grow, acc, polyline in nd.layeriter(xs_logic):
                    (a1, b1), (a2, b2), c1, c2 = grow
                    if shape == 'circle':
                        outline = geom.circle(radius=0.5*length, N=32)
                    else:
                        #outline = geom.box(width=width+(b1-b2), length=length)
                        outline = [(0, a2*width+b2), (0, a1*width+b1), (length, a1*width+b1), (length, a2*width+b2)]
                    if abs((a2*width+b2) - (a1*width+b1)) > 0.0001:
                        nd.Polygon(layer=lay, points=outline).put(0, 0, 180)

                    if shape not in stubshapes:
                        logger.warning("stub shape '{}' not recognized, possible options are {}.".
                            format(shape, stubshapes))
    stubs[name] = C
    return name  # name can have change to default_xs


def put_stub(
    pinname=None,
    length=None,
    shape='box',
    pinshape=None,
    pinsize=None,
    pinlayer=None,
    pinshow=True,
    annotation_layer=None,
    pinstyle=None,
    cell=None,
):
    """Add a xsection stub to one or more pins.

    The pinstyle is a key poining to a dict which contains the setting for the stub.
    These setting are overruled by the explicit keyword in the call, e.g. shape.

    The method put_stub() places a stub cell (pinshape cell + stub layers) on a Node if
    stub information (see below) has been defined for the xsection assigend to the pin.
    If no xsection has been defined for the pin, only the pincell will the created and used.
    If an annotation_layer is defined, the annotation is added to the parent cell.

    The xsection and width of a stub are obtained from its pin attributes.
    If no attributes are set, the stub layers are empty and the stubcell
    only contains the pincell. It is possible to supply a list of pins.

    The pinstyle is determined by the following sequence:
    1. use the pinstyle of the xsection, if set
    2. use the pinstyle 'default'
    3. update the pinstyle with explicit settings via the keyword parameters

    The pinshape is one of the following:

    1. A reference to a key in the cfg.pinshapes dictionary which contains
        the raisepin polygon description of which a pin cell is created.

    2. A cell.

    For a pin having attribute xs=None, only the pin is drawn.

    Args:
        pinname (str | list of str | None): name(s) of the pin(s) (default = None)
            For pinname=None (default) all chain-pins will obtain a stub.
        length (float): length of the stub (thickness into cell)
        shape (string | Cell): shape of the stub 'box' | 'cirlce' (default = 'box')
        pinshape (string): pinshape used in the stub
        pinsize (float): scaling factor of the pinshape (default = 1)
        layer (str | int | (int, int)): pin layer
        annotation_layer (str | int | (int, int)): annotation layer
        pinstyle (str): pin style name used to look up a dict with pinstyle settings.
            Default it will look in the xsection pinstyle setting. If that is None
            the pinstyle 'default' is used.
        pinshow (bool): set the pin show attribute (default=True)
        cell (str | Cell): reference to a cell to place the stub in.
            Default=None, which adds the stub to the active cell.

    Returns:
        None
    """
    if cell is not None:
        # add user defined cell as temporary active cell:
        if isinstance(cell, str):
            if cell in cfg.cellnames.keys():
                cfg.cells.append(cfg.cellnames[cell])
            else:
                nd.main_logger(f"Trying to put a stub in unknown cell '{cell}'.", "error")
        else:
            cfg.cells.append(cell)
        nd.cfg.patchcell = True
    _cell = cfg.cells[-1]  # active cell

    # prepare pin:
    if pinname is None or pinname == []:
        pinname = [name for name, pin in _cell.pin.items() if pin.chain == 1]
    elif isinstance(pinname, str):
        pinname = [pinname]
    elif isinstance(pinname, nd.Node):
        pinname = [pinname.name]

    # prepare stub
    if isinstance(shape, nd.Cell):
        stubcell = shape
        shape = shape.cell_name
    else:
        stubcell = None

    # prepare pinstyles
    stylestr = None  #'default'
    if cfg.validation_pinstyle is not None:
        pinstyle = cfg.validation_pinstyle

    pinstylenames = cfg.pinstyles.keys()
    if pinstyle in pinstylenames:
        stylestr = pinstyle
    elif pinstyle is not None:
        logger.warning("pinstyle '{}' unknown. Using default instead. " \
              "Known pinstyles: {}".format(pinstyle, pinstylenames))

    #prepare annotation layer
    _annotation_layer = annotation_layer

    for p in pinname:
        _stylestr = stylestr
        try:
            node = _cell.pin[p]
            node.show = pinshow # display pin name in layout and fingerprint
            if stubcell is not None: # use existing Nazca cell for stub+pin
                name = "stub_{}".format(shape)
                if name not in stubs.keys():
                    stubs[name] = stubcell
                stubs[name].put(node)
            else:
                if node.xs is None: # pin cell only, no stub
                    xs = None
                    width = 0
                else: # define stub + pin from style
                    xs = node.xs
                    if cfg.validation_pinstyle is None and _stylestr is None:
                        xs_pinstyle = nd.get_xsection(xs).pinstyle
                        if xs_pinstyle in pinstylenames and pinstyle is None:
                            _stylestr = xs_pinstyle
                        else:
                            _stylestr = 'default'
                    if node.width is None:
                        width = 0
                    #TODO: why allow str widths here?
                    elif isinstance(node.width, str):
                        width = literal_eval(node.width)
                    else:
                        width = node.width

                # make and place cell
                if xs is None: # pin only
                    arrow = make_pincell(
                        layer=pinlayer,
                        shape=pinshape,
                        size=pinsize,
                        style=_stylestr,
                    )
                    inst = arrow.put(node)
                else: # stub
                    if length is not None:
                        if isinstance(length, (float, int)):
                            stub_length = float(length)
                    else:
                        stubstyle = _stylestr
                        if stubstyle is None:
                            stubstyle = 'default'
                        stub_length = nd.cfg.pinstyles[stubstyle].get('stub_length', 0)
                        if stub_length is None:
                            stub_length = 0.0
                    name = _makestub(
                        xs_guide=xs,
                        width=width,
                        length=stub_length,
                        shape=shape,
                        pinshape=pinshape,
                        pinsize=pinsize,
                        pinlayer=pinlayer,
                        pinstyle=_stylestr,
                    )

                    # place stubs without drc check
                    drc = cfg.drc
                    cfg.drc = False
                    if node.io & 1:
                        inst = stubs[name].put(node.rot(180), flip=True)
                    else:
                        inst = stubs[name].put(node.rot(180))
                    inst.sourcepin = node  # for accessing pin properties for cell visualisation
                    cfg.drc = drc


            # add annotation (keep outside of stub and pin cell, as these are reused for different pin names)
            if _stylestr is None:
                _stylestr = 'default'
            if _annotation_layer is None:
                annotation_layer = nd.cfg.pinstyles[_stylestr].get('annotation_layer', False)
            else:
                annotation_layer = _annotation_layer
            # Move the annotation location for readability
            drc = cfg.drc  # switch off drc check for annotations
            cfg.drc = False
            if annotation_layer:
                anno = nd.Annotation(layer=annotation_layer, text=p)
                anno.put(node.move(*nd.cfg.pinstyles[_stylestr]['annotation_move']))
                anno.sourcepin = node  # for accessing pin properties for cell visualisation
            if cfg.store_pin_attr:  # TODO: check how this overlaps with yaml pinstyle for validation.
                if node.type != 'bbox':
                    text = f"name: {node.name}\nwidth: {node.width}\nxs: {node.xs}\nxya: {node.fxya(digits=6)}\nflip: {node.io}\nremark: {node.remark}"
                    nd.Annotation(layer='bb_pin_attr', text=text).put(node)
            cfg.drc = drc

        except Exception as E:
            mes = "Can't add stub in cell '{}' on pin '{}'.".\
                format(_cell.cell_name, p)
            logger.exception(mes)
            raise Exception(E)
    if cell is not None:
        nd.cfg.patchcell = False
        cfg.cells.pop(-1)


def validation_style(
    on=True,
    stub_layer=None,
    pin_style=None,
    validation_layermapmode="none",
    validation_layermap= None,
):
    """Set layers and pinstyle and layermap for BB validation stubs.

    Args:
        on (bool): default=True, turn on validation style
        stub_layer (int, (int, int), str): layer of the stub.
        pin_style (dict): optional pinstyle dict for the arrow and stub.
        validation_layermapmode (str):
        validation_layermap (dict): {in_layer: out_layer}

    Returns:
        None
    """
    if not on:  # reset validation settings
        cfg.validation_pinstyle = None
        cfg.validation_stub_xs = None
        cfg.validation_layermapmode = "none"
        cfg.validation_layermap = None
        return None

    if stub_layer is None:
        stub_layer = cfg.default_layers['validation_stub']

    nd.add_xsection(name="validation_stub")
    nd.add_layer(name="validation_stub", layer=stub_layer)
    nd.add_layer2xsection(xsection="validation_stub", layer=stub_layer)

    nd.cfg.store_pin_attr = True

    if pin_style is None:
        validation_pinstyle = {
            'shape': 'arrow_full',
            'size': 1.0,
            'layer': "bb_pin",  # use string name
            'annotation_layer': 'bbox_pin_text',  # use string name
            'annotation_move': (0, 0),
            'stub_length': 1.0,
        }
        nd.add_pinstyle('validation', validation_pinstyle)
    else:
        pin_style['annotation_move'] = (0, 0)  # always force to pin location
        nd.add_pinstyle('validation', pin_style)

    # set global parameters that for overrideing the stub xs and the pinstyle:
    cfg.validation_pinstyle = "validation"
    cfg.validation_stub_xs = "validation_stub"

    cfg.validation_layermapmode = validation_layermapmode
    cfg.validation_layermap = validation_layermap


def export_bb2png(cellcalls, path='', close=True, bbox_pins=True):
    """Create png image for all cells in <cellcalls>.

    Only the Cell objects are used from <cellcalls>, but it does accept as well
    the dict as needed in method 'bb_fingerprint'.

    cellcalls (dict | list of Cells): All the cells {function_call_str: cell_from_call} | list of Cells
    close (bool): close cell immediately after drawing (default = True)
    path (str): path for saving png files
    bbox_pins (bool): add pins to the bbox (default = True)

    Returns:
        None
    """
    export = True

    if isinstance(cellcalls, dict):
        cells = cellcalls.values()
    else:
        cells = cellcalls

    N = len(cells)

    for i, cell in enumerate(cells):
        print('{}/{}'.format(i+1, N))
        bbox_stubs_store = cfg.bbox_stubs
        cfg.bbox_stubs = bbox_pins
        if export:
            output = 'file'
        else:
            output = None
        nd.export_plt(cell, output=output, path=path,
            title=cell.cell_paramsname)
        cfg.bbox_stubs = bbox_stubs_store
        plt.close("all")
    return None


def bb_fingerprint(cellcalls, save=False, filename='fingerprint.json', infolevel=0):
    """Generate a dict with parameter list and pin list of a list off BBs.

    This generates a blueprint (aka fingerprint).

    Args:
        cellcalls (dict): All the cells {function_call_str: cell_from_call}

    Returns:
        dict: Dictionary with building block info, Structure is as follows::

            {'<cellname>':
                {'pins':
                    { '<pinname>':
                        {'width': <value>,
                         'xs'   : <value>
                         'xya'  : (<x>, <y>, <a>)
                         'show' : <values>
                        }
                        , ...
                    },
                 'parameters':
                     {'<keyword>': <value>}
                     , ...
                 'paramsname'   : <value>,
                 'basename'     : <value>,
                 'groupname'    : <value>,
                 'module'       : <value>,
                 'function'     : <value>,
                 'length'       : <value>,
                 'width'        : <value>
                }
            }
    """
    fingerprint = {}
    cellcallsorted = sorted(cellcalls.items()) # not effective for Python <3.6 with scrambled dict.
    for callstr, cell in cellcallsorted:
        if hasattr(cell, 'group'):
            group = cell.group
        else:
            group = ''

        module = cell.module

        pinsetting = {}
        if infolevel > 0:
            logger.info(cell.cell_basename)
        for name, pin in sorted(cell.pin.items()):
            if name != 'org':
                pinsetting[name] = {
                    'width': pin.width,
                    'xsection': pin.xs,
                    'pointer': [float('{:.5f}'.format(val+1e-10)) for val in pin.xya()],
                    'show': pin.show,
                    'type': pin.type,
                    'remark': pin.remark
                }
                if infolevel > 0:
                    logger.info(name,  pinsetting[name])

        paramsetting = {}
        try: # 'parameters' may be empty
            for p in cell.parameters:
                paramsetting[p] = cell.parameters[p]
        except:
            logger.debug("Empty parameters in {cell.cell_name}")

        #bbox
        try:
            length = cell.length
        except:
            logger.debug("Empty length in {cell.cell_name} setting it to None.")
            length = 'None'
        try:
            width = cell.width
        except:
            logger.debug("Empty width in {cell.cell_name} setting it to None.")
            width = 'None'

        cellnameparameters = cell.properties.get('cellnameparameters', None)
        fingerprint[cell.cell_name] = {
            'pins': pinsetting,
            'parameters': paramsetting,
            'paramsname': cell.cell_paramsname,
            'cellnameparameters': cellnameparameters,
            'basename': cell.cell_basename,
            'group': group,
            'modulename': module,
            'function': callstr,
            'length': length,
            'width': width,
            'connect': (cell.default_in, cell.default_out)}

    if save:
        try:
            with open(filename, 'w') as fp:
                json.dump(fingerprint, fp, indent=2, sort_keys=True)
            print("Saved fingerprint to {}".format(filename))
        except Exception as E:
            print(f"ERROR in json.dump of fingeprint: {E}")
            nd.main_logger(f"ERROR in json.dump of fingeprint: {E}", "error")
    if infolevel > 1:
        pprint(fingerprint)
    return fingerprint


def validate_black_to_white_mapping(black2whiteMap, allBBcells, infolevel=0):
    """Validate if all white and black boxes are mapped.

    By increasing the infolevel more information will be displayed on where
    black <-> white mappings are missing.

    Args:
        black2whiteMap (dict): {blackbox-basename: whitebox-function}
        allBBcells (list of Cell): all black box cells.
        infolevel: amount of feedback to stdout (default = 0)

    Returns:
        bool: True if mapping is successful.
    """

    if infolevel > 0:
        print('Validate black <-> white box mappings:')

    okay = True
    mesg = []
    #print("\nStatus of white-box implementation:")
    mappings = ""
    header = "   {0:30}{1}".format("BASENAME", "FUNCTION")
    black_in_map = list(black2whiteMap.keys())
    basenames = {cell.cell_basename: cell for cell in allBBcells}
    black_count = len(basenames)
    white_found = 0
    for basename, cell in sorted(basenames.items()):
        func = ""
        defined = '?'
        if basename in black_in_map:
            defined = '*'
            white_found += 1
            func = black2whiteMap[basename]
        mappings += "{0}  {1:30}{2}\n".format(defined, basename, func)

    black_orphan_count = black_count - white_found
    text = ("Found {} black orphans".format(black_orphan_count))
    logger.info("{}".format(text))
    if black_orphan_count > 0:
        okay = False
        mesg.append(text)
        if infolevel > 0:
            print(header)
            print(mappings)

    white_orphans = set(black2whiteMap) - set(basenames.keys())
    white_orphan_count = len(white_orphans)
    text = "Found {} white orphans".format(white_orphan_count)
    logger.warning("{}".format(text))
    if white_orphan_count > 0:
        okay = False
        mesg.append("{}{}".format(text, '.'))
        if infolevel > 0:
            logger.info(header)
            for name in white_orphans:
                logger.info("?  {:30}{}".format(name, black2whiteMap[name]))

    if not okay:
        logger.error("Incorrect white <-> black box mappings.")
        for m in mesg:
            logger.error("- {}{}".format(m, '.'))
    else:
        logger.info("\nSuccesful black <-> white mapping.")


bbox_pinnames = [
    'lb', 'lc', 'lt',
    'tl', 'tc', 'tr',
    'rt', 'rc', 'rb',
    'br', 'bc', 'bl',
    'cc']

def put_boundingbox(pinname, length, width, raise_pins=True, outline=True,
        align='lc', name=True, params=True, move=(0, 0, 0), name_layer=None):
    """Create bounding box (bbox) cell inside the active cell.

    This function places a bbox cell and raises by default the bbox pins into
    the active cell. The bbox displays a bbox outline (can be switched off).
    By default it also adds the active cellname and parameters.

    Args:
        pin (str): pin to place bbox on (center left of bbox)
        length (float): length of the bbox
        with (float): width of the bbox
        raise_pins (bool): raise bbox pins into active cell (default = True)
        outline (bool): draw bbox outline (default = True)
        align (str): align the bbox on the specified bbox pin <pinname> (default = 'lc')
        name (bool): display the (active) cell name in the bbox. (default = True)
        params (bool): add parameter annotation to the bbox
        move (tuple): move the bbox placement by (float, float, float)
        name_layer (str): name of the layer to be used for bbox name

    Returns:
        None
    """
    cell = cfg.cells[-1]
    _paramsname = cell.cell_paramsname
    _parameters = cell.parameters
    with nd.Cell(name='bbox', instantiate=False) as C:
        outline = geom.box(length, width)
        nd.Polygon(layer='bbox', points=outline).put(0)

        stbs = cfg.bbox_stubs and outline
        nd.Pin('lb', type='bbox', show=stbs).put(0, -0.5*width, 180)
        nd.Pin('lc', type='bbox', show=stbs).put(0, 0, 180)
        nd.Pin('lt', type='bbox', show=stbs).put(0, 0.5*width, 180)

        nd.Pin('tl', type='bbox', show=stbs).put(0, 0.5*width, 90)
        nd.Pin('tc', type='bbox', show=stbs).put(0.5*length, 0.5*width, 90)
        nd.Pin('tr', type='bbox', show=stbs).put(length, 0.5*width, 90)

        nd.Pin('rt', type='bbox', show=stbs).put(length, 0.5*width, 0)
        nd.Pin('rc', type='bbox', show=stbs).put(length, 0, 0)
        nd.Pin('rb', type='bbox', show=stbs).put(length, -0.5*width, 0)

        nd.Pin('br', type='bbox', show=stbs).put(length, -0.5*width, -90)
        nd.Pin('bc', type='bbox', show=stbs).put(0.5*length, -0.5*width, -90)
        nd.Pin('bl', type='bbox', show=stbs).put(0, -0.5*width, -90)

        nd.Pin('cc', type='bbox', show=cfg.bbox_stubs).put(0.5*length, 0, 0)

        if stbs:
            pin_placement= {
                'bbox_left': ['lt', 'tr', 'rb', 'bl'],
                'bbox_right': ['lb', 'tl', 'rt', 'br'],
                'bbox_center': ['lc', 'tc', 'rc', 'bc']
                }
            for style, pins in pin_placement.items():
                nd.put_stub(
                    pins,
                    pinshape=cfg.pinstyles[style]['shape'],
                    pinsize=cfg.pinstyles[style]['size'],
                    pinlayer=cfg.pinstyles[style]['layer'],
                    annotation_layer=cfg.pinstyles[style]['annotation_layer']
                )

        if name:
            cellname(
                cellname=_paramsname,
                length=length,
                width=width,
                align='cc',
                name_layer=name_layer,
            ).put(C.pin['cc'])
        if params:
            put_parameters(parameters=_parameters, pin=C.pin['cc'])

    # options to align the bbox in this cell w.r.t. to its pins:
    align_shift = {
        'lb': (0, 0.5*width),
        'lc': (0, 0),
        'lt': (0, -0.5*width),
        'rb': (-length, 0.5*width),
        'rc': (-length, 0),
        'rt': (-length, -0.5*width),
        'bc': (-0.5*length, 0.5*width),
        'cc': (-0.5*length, 0),
        'tc': (-0.5*length, -0.5*width)
    }
    dx = align_shift[align][0]
    dy = align_shift[align][1]

    box = C.put(cfg.cells[-1].pin[pinname].move(dx, dy).move(*move))
    cell.length = length
    cell.width = width
    cell.add_property({
        'dimension': {
            'length': length,
            'width': width,
            'area': length*width}
        })

    if raise_pins:
        box.raise_pins(bbox_pinnames, show=False)

    return None


def load_gdsBB(
    gdsfilename,
    cellname,
    pinfilename=None,
    newcellname=None,
    layermap=None,
    cellmap=None,
    flip=False,
    flop=False,
    scale=None,
    stubs=True,
    native=True,
    bbox=False,
    bboxbuf=0,
    hull=None,
    prefix='',
    suffix='_stubs',
    instantiate=None,
    flat=False,
    layermapmode=None,
):
    """Load a gds cell and the corresponding pin info from file.

    This function uses method 'load_gds' for loading the gds file
    and adds pins and stubs to the cell according to the pinfile.
    This creates building block cell with connectivity.
    If there is not pinfile available, then use method 'load_gds" instead.

    Format of the pin file, including the single line header:

    port, x, y, a, xs, w
    a0, 100, 0, 90, nazca, 2.0
    ...

    If stubs == True then both the stubs and the loaded GDS are instantiated
    inside the cell returned by this function.
    This to ensure that for this cell that the gds and stubs remain aligned
    under flipping and flopping.

    Args:
        gdsfilename (str): gds filename to import from
        cellname (str): cellname to import
        pinfilename (str): optional name of csv file containing pin info for stubs.
        newcellname (str): new cell name
        cellmap (dict): mapping of cellnames {cellname_in: celllname_out}
        layermap (dict): mapping of layer {number_in: number_out}
        flip (bool): mirror the gds without changing pin direction (default = False)
        flop (bool): mirror the gds with reversing pin direction (default = False)
        scale (float): scaling factor (default = 1.0). Use with great care
            as scaled building blocks will not make much sense from a functional
            prespective.
            pin positions will be scaled along with the gds, however,
            xs information, like 'width' will *not* be scaled
        stubs (bool): place stubs in the cell (default = True).
            The stub value also set the instantiation value of the returned cell
        native (bool): load gds into a Nazca native format (default = True)
        bbox (bool): add a bounding box if True (default = False)
        bboxbuf (float): add a buffer layer around the structure and the bounding box.
        prefix (str): prefix to add the cell name to avoid cellname clashes.
            See help(nazca.load_gds)) more info on cellname resolution
        suffix (str): suffix to indicate it is cell wuth stubs added (default = '_stubs')
        instantiate (bool): intantiate the returned cell

    Returns:
        Cell: Nazca Cell with the loaded gds and (if provided) pins and stubs.
    """
    if stubs:
        _instantiate = True
        bboxgds = False
        topbbox = bbox
    else:
        _instantiate = False
        bboxgds = bbox
        topbbox = False

    # instantiate override:
    if instantiate is not None:
        instantiate = instantiate
    else:
        instantiate = _instantiate

    with nd.Cell(prefix+cellname+suffix, instantiate=instantiate, autobbox=topbbox) as C:
        #TODO: if instantiate is true the instance and gds file may get the same name
        #      causing a topology error in GDS.
        nd.load_gds(
            filename=gdsfilename,
            cellname=cellname,
            newcellname=newcellname,
            layermap=layermap,
            layermapmode=layermapmode,
            cellmap=cellmap,
            scale=scale,
            native=native,
            bbox=bboxgds,
            bboxbuf=bboxbuf,
            hull=hull,
            prefix=prefix,
            instantiate=instantiate,
            flat=flat,
        ). put(0, flip=flip, flop=flop)
        if flip:
            s = -1.0
        else:
            s = 1.0
        if pinfilename is not None:
            if not os.path.exists(pinfilename):
                nd.main_logger(f"Pin file not found '{pinfilename}'. No pins and stubs placed.", "error")
            else:
                df = pd.read_csv(pinfilename, delimiter= ',', skiprows=1, header=None,
                    names = ['io', 'x', 'y', 'a', 'xs', 'w'])
                df.fillna(0.0, inplace=True)
                for row in df.itertuples():
                    i, io, x, y, a, xs, w = row
                    try:
                        xs = xs.strip()
                        if xs.upper() == 'NONE':
                            xs = None
                    except:
                        logger.exception('xsection not parseble.')
                        xs = None
                    try:
                        w = float(w)
                    except:
                        if w is not None and w != 'None':
                            logger.warning("load_gdsBB: Could not cast width to float on pin '{}': {}".\
                                format(io, w))
                        w = 0.0
                    if scale is None:
                        scale = 1.0
                    nd.Pin(name=io, xs=xs, width=w).put(scale*x, s*scale*y, a)
                    if stubs: #and xs is not None :
                        put_stub(io)
        if instantiate and bbox:
            C.autobbox = True
            C.bboxbuf = bboxbuf
    return C


def _PIL2array(img):
    return np.array(img.getdata(), bool).reshape(img.size[1], img.size[0])


def image(
    name,
    layer=1,
    size=256,
    pixelsize=1,
    threshold=0.5,
    cellname=None,
    invert=False,
    align="cc",
    box_layer=None,
    box_buf=0,
    merge=True,
    bg_white=True
):
    """Read an image file and return a nazca cell with the image.

    Image format can be png, jpg, gif, bpm, eps and others, as supported by Pillow.
    Note that the output resolution (size) does not exceed the image resolution.
    Increase <pixelsize> in this case to obtain a larger image in gds.

    A rectangular box can be added around the image by providing a <box_layer>.
    This box can be enlarged beyond the original image size by setting <box_buf> > 0.

    Args:
        name (str): name of the image file
        layer (int): layer number that the image will be written to (default 1)
        size (int): maximum bounding box size in pixels (default 256). Setting this to
            zero keeps the natural size of the image. This may be slow and can result
            in a large GDS file, depending on the image size.
        pixelsize (float): pixel size in micron (default 1)
        threshold (float): black/white threshold (default 0.5)
        cellname (str): Nazca cell name (default image filename)
        invert (bool): flag to invert black & white (default False)
        align (str): two character string for image alignment (default 'cc')
            allowed:
            lt, ct, rt,
            lc, cc, rc,
            lb, cb, rb

        box_layer (str | int | tuple): layer reference to generate a rectangular
            box behind the image, e.g. for tiling exclusion areas (NOFILL)
            (default = None)
        box_buf (float): extra buffer for the box_layer in um
        merge (bool): merge polygons after generating image (merge is True) or leave
            many small polygons (merge is False), which is faster but results in a
            somewhat larger file size.
        bg_white (bool): transparent background (if present) is interpreted as
            white (True, default) or black (False).

    Returns:
        Cell: cell with image

    Examples:
        Load an image in a cell and put and/or export it to gds::

            import nazca as nd

            logo = nd.image('mylogo.png', align='lb') # left/bottom alignment
            logo.put(0)
            nd.export_gds()

        or::

            import nazca as nd

            nd.export_gds(logo, filename='mylogo.gds')
    """
    if cellname is None:
        cellname = os.path.basename(name)
    p = pixelsize
    threshold = int(threshold * 255)
    a = {"lb", "cb", "rb", "lc", "cc", "rc", "lt", "ct", "rt"}
    if align not in a:
        mes = "Invalid alignment specification '{}' for image '{}'.".format(align, name)
        mes += "Allowed values are {}.".format(a)
        mes += "Using default value 'cc'"
        logger.warning(mes)
        align = "cc"
    halign = {"l": 0, "c": -0.5, "r": -1}
    valign = {"b": 0, "c": -0.5, "t": -1}

    # Treat transparent background as white (usually good) or black.
    im = Image.open(name).convert("RGBA")
    new = Image.new("RGBA", im.size, "WHITE" if bg_white else "BLACK")
    new.paste(im, mask=im)
    gray = new.convert("L")
    # resize keep aspect, only if smaller
    gray.thumbnail((size, size), Image.LANCZOS)
    bw = gray.point(lambda x: 0 if x < threshold else 255, "1")
    pix = _PIL2array(bw)
    width, height = bw.size
    width_tot = width * p + 2 * box_buf
    height_tot = height * p + 2 * box_buf
    logger.info(
        "Generating {}x{} pixels image of {:.0f}x{:.0f} um2, edge is {} um.".format(
            width, height, width * p, height * p, box_buf
        )
    )
    h0 = halign[align[0]] * width_tot
    v0 = valign[align[1]] * height_tot
    # Ensure pixel alignment on the minimum GDS grid size (default 1 nm)
    h0 = nd.gds_base.round_to_db_unit(h0) * nd.gds_base.gds_db_user
    v0 = nd.gds_base.round_to_db_unit(v0) * nd.gds_base.gds_db_user
    box_buf = nd.gds_base.round_to_db_unit(box_buf) * nd.gds_base.gds_db_user

    polygons = []

    with nd.Cell(cellname) as C:
        for line in range(height):
            x1 = x0 = h0 + box_buf
            lb = lw = 0
            y0 = (height - line) * p + v0 + box_buf
            for pixel in pix[line]:
                if pixel == invert:
                    lb += 1
                    if lw > 0:
                        x0 += lw * p
                        lw = 0
                else:
                    lw += 1
                    if lb > 0:
                        x1 = x0 + lb * p
                        xy = [(x0, y0), (x1, y0), (x1, y0 - p), (x0, y0 - p)]
                        polygons.append(xy)
                        x0 = x1
                        lb = 0
            if lb > 0:
                x1 = x0 + lb * p
                xy = [(x0, y0), (x1, y0), (x1, y0 - p), (x0, y0 - p)]
                polygons.append(xy)
        if merge:
            if polygons:
                polygons = merge_polygons(polygons)
            else:
                logger.warning(f"Resulting image is empty for '{name}'.")
        for pol in polygons:
            nd.Polygon(points=pol, layer=layer).put()

        if box_layer is not None:
            nd.Polygon(
                layer=box_layer,
                points=[
                    (0, 0),
                    (0, height_tot),
                    (width_tot, height_tot),
                    (width_tot, 0),
                ],
            ).put(h0, v0, 0)
    return C


def vector_image(
    name,
    layer=1,
    width=None,
    height=None,
    threshold=0.5,
    cellname=None,
    invert=False,
    align="cc",
    box_layer=None,
    box_buf=0,
    bg_white=True,
):
    """Read an image file and return a nazca cell with the vector image.

    Image format can be png, jpg, gif, bpm, eps and others, as supported by Pillow.

    A rectangular box can be added around the logo by providing a box_layer.
    This box can be enlarged beyond the original image size by setting box_buf > 0.

    Args:
        name (str): name of the image file
        layer (int): layer number that the image will be written to (default 1)
        width (float): width of image. If height is not specified, the aspect ratio is
            maintained.
        height (float): height of image. If width is not specified, the aspect ratio is
            maintained.
        threshold (float): black/white threshold (default 0.5)
        cellname (str): Nazca cell name (default image filename)
        invert (bool): flag to invert black & white (default False)
        align (str): two character string for image alignment (default 'cc')
            allowed:
            lt, ct, rt,
            lc, cc, rc,
            lb, cb, rb

        box_layer (str | int | tuple): layer reference to generate a rectangular
            box behind the image, e.g. for tiling exclusion areas (NOFILL)
            (default = None)
        box_buf (float): extra buffer for the box_layer in um
        bg_white (bool): transparent background (if present) is interpreted as
            white (True, default) or black (False).

    Returns:
        Cell: cell with image

    Examples:
        Load an image (typically a logo) in a cell and put and/or export it to gds::

            import nazca as nd

            logo = nd.vector_image('mylogo.png', align='lb') # left/bottom alignment
            logo.put(0)
            nd.export_gds()

        or::

            import nazca as nd

            nd.export_gds(logo, filename='mylogo.gds')
    """
    if cellname is None:
        cellname = os.path.basename(name)
    try:
        import potrace
    except ModuleNotFoundError:
        logger.warning(
            "Can't load module 'potrace'. Function vector_image() not available.\n"
            "Install 'pypotrace' or use (bitmap) nd.image()."
        )
        # Return empty cell.
        with nd.Cell(f"{cellname}_NoVector") as C:
            pass
        return C
    threshold = int(threshold * 255)
    a = {"lb", "cb", "rb", "lc", "cc", "rc", "lt", "ct", "rt"}
    if align not in a:
        mes = f"Invalid alignment specification '{align}' for image '{name}'."
        mes += f"  Allowed values are {a}."
        mes += "  Using default value 'cc'"
        logger.warning("\n".join(mes))
        align = "cc"
    halign = {"l": 0, "c": -0.5, "r": -1}
    valign = {"b": 0, "c": -0.5, "t": -1}

    # Treat transparent background as white (usually good) or black.
    im = Image.open(name).convert("RGBA")
    new = Image.new("RGBA", im.size, "WHITE" if bg_white else "BLACK")
    new.paste(im, mask=im)
    gray = new.convert("L")
    if invert:
        bw = gray.point(lambda x: 0 if x < threshold else 255, "1")
    else:
        bw = gray.point(lambda x: 0 if x >= threshold else 255, "1")
    imwidth, imheight = bw.size

    # Scalings
    scalew, scaleh = None, None
    if width:
        scalew = width / imwidth
    if height:
        scaleh = height / imheight
    if not scalew:
        scalew = scaleh
    if not scaleh:
        scaleh = scalew
    if not scaleh and not scalew:
        scalew, scaleh = 1, 1
    width = imwidth * scalew
    height = imheight * scaleh
    width_tot = width + 2 * box_buf
    height_tot = height + 2 * box_buf
    logger.info(
        f"Generating image of {width:.0f}x{height:.0f} um2, edge is {box_buf} um."
    )
    h0 = halign[align[0]] * width_tot
    v0 = valign[align[1]] * height_tot

    pix = _PIL2array(bw)  # 2D numpy array
    bmp = potrace.Bitmap(np.flipud(pix))
    path = bmp.trace()
    pols = []
    sc = 1e8  # Scale to prevent rounding
    for curve in path:
        # Convert bezier to points, scale and move.
        pols.append(
            [
                [(x * scalew) * sc, (y * scaleh) * sc]
                for x, y in curve.tesselate(curve.adaptive)
            ]
        )
    # Handle polygons with holes.
    pols = nd.clipper._clipper2GDS(pols)  # doesn't need clipper in spite of the name
    # Write out the polygons
    with nd.Cell(cellname) as C:
        for pol in pols:
            nd.Polygon(points=[[x / sc, y / sc] for x, y in pol], layer=layer).put(
                h0 + box_buf, v0 + box_buf, 0
            )
        if box_layer is not None:
            nd.Polygon(
                layer=box_layer,
                points=[
                    (0, 0),
                    (0, height_tot),
                    (width_tot, height_tot),
                    (width_tot, 0),
                ],
            ).put(h0, v0, 0)
    return C


# =============================================================================
# predefined example cell for documentation examples
# =============================================================================
def example_cell():
    """Create a Nazca example cell.

    Returns:
        cell
    """
    with nd.Cell('example_cell') as example:
        s = nd.strt(length=20).put(0)
        b = nd.bend(angle=45).put()
        nd.Pin('a0', pin=s.pin['a0']).put()
        nd.Pin('b0', pin=b.pin['b0']).put()
    return example

dir_path = os.path.dirname(os.path.abspath(__file__))

def nazca_logo(layers={'ring':1, 'box':2, 'bird':3}, cellname='nazca_logo',
        scale=1.0, bbox=True):
    """Return the nazca logo as a cell.

    Args:
        layers (dict): gds layers of the logo. default = {'ring':1, 'box':2, 'bird':3}
        cell_name (str): default = 'nazca_logo'
        scale (float): scaling factor (default = 1)
    """
    return nd.load_gds(
	    filename=os.path.join(dir_path, 'nazca_logo.gds'),
        instantiate=True,
        cellname='nazca',
		newcellname=cellname,
        scale=scale,
        bbox=bbox,
	    layermap={1:layers['ring'], 2:layers['box'], 3:layers['bird']})

