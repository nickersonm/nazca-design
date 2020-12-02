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
# Test replacement of gds cells from one file with those from another.
# This will be used for replacing black-box building blocks with their real
# implementation.
#
# (c) 2016-2020 Katarzyna Lawniczuk, Ronald Broeke
#==============================================================================

"""
Module for replacing static cells and parametric cells.
"""

import os
from ast import literal_eval
import time
from collections import Counter
from pprint import pprint
import nazca as nd
from nazca.logging import logger


def _exclude_substrings(allBlackNames):
    """Check for blackbox names with overlapping substrings at start.

    Args:
        allBlackNames (list of str): list of all blackbox names in the libraries

    Returns:
        list of list: list per blackbox name containing which extended blackbox
            names to exclude in cellmap match
    """
    b = sorted(allBlackNames)
    excludes = {}
    Lb = len(b)
    for i in range(Lb):
        exclude = []
        for j in range(i, Lb-1):
            if b[j+1].startswith(b[j]):
                exclude.append(b[j+1])
            else:
                break
        excludes[b[i]] = exclude
    return excludes


def _findBlackCells(blackBasename, cells=[], excludes=None, infolevel=0):
    """Find all "black" cells starting with <blackBasename> in list <cells>.

    Args:
        blackBasename (str): string matching first part (basename) of black cells
        cells (list of str): list of cell names
        exclude (list): list of names to exclude due to matching substrings

    Returns:
        list of str: list of cell names in <cells> starting with <blackBasename>
    """

    matchingcells = []
    for name in cells:
        if name.startswith(blackBasename):
            exclude = False
            for e in excludes[blackBasename]:
                if name.startswith(e):
                    exclude = True
                    break
            if not exclude:
                matchingcells.append(name)

    if infolevel>0:
        logger.info("* basename={}, excludes={}\n  matching-cells:{}".\
            format(blackBasename, excludes[blackBasename], matchingcells))
    return matchingcells


def _readPcellParams(
    blackBasename,
    gdsin,
    allInputCellNames,
    excludes=None,
    infolevel=0,
):
    """Create a dict off cells each containing a dict of all Pcell parameter value.

    Args:
        blackBasename (str): pcell base name (without parameter string suffix)
        allInputCellNames (list of str): list of cell names
        gdsin (str): gds input file name
        infolevel (int): amount of debug info that is printed, 0 is minimum.

    Returns:
        dict: {cellname: {<var1>: <value1>, <var2>: <value>2}},
            cellnames and their parameters dict
    """
    allCellsOf1type = _findBlackCells(
        blackBasename=blackBasename,
        cells=allInputCellNames,
        excludes=excludes,
        infolevel=infolevel
    )
    if infolevel > 1:
        logger.info("blackbox '{}'\ninstances:".format(blackBasename))
        pprint(allCellsOf1type)
    allPcellParams = {} #  dictionary of {cell: PcellStrParams}
    for cellname in allCellsOf1type:
        PcellStrParams = {}
        PcellValueParams = {}

        for lay, pos, text in nd.get_cell_annotation(gdsin.cells[cellname]):
            if 'Parameter' in text or nd.get_layer(lay) == 'bb_parameter_text':
                PcellStrParams = dict(**PcellStrParams, **nd.string_to_parameters(text))

        for varname, valuestr in PcellStrParams.items():
            valuestr = str(valuestr)  # needed for 'yaml' type
            value_lit = literal_eval(valuestr)
            if valuestr.upper() == 'TRUE':
                value = True
            elif valuestr.upper() == 'FALSE':
                value = False
            elif isinstance(value_lit, (int, float, tuple)):
                value = value_lit
            else:
                logger.info('ERROR: parameter value not recognized: {}'.format(valuestr))

            PcellValueParams[varname] = value

        allPcellParams[cellname] = PcellValueParams
    if infolevel > 2:
        logger.info("\nAll Pcells of '{}' as {{name:parameters}}".format(blackBasename))
        pprint(allPcellParams)
    return allPcellParams


def _createWhitePcellLib(
    gdsin,
    whitelibrary=None,
    black2whiteMap=None,
    prefixb='black_',
    prefixw='wb_',
    infolevel=0,
    layermap=None,
    layermapmode=None,
):
    """Create a static whitebox gds library from black Pcells found in <gdsin>.

    Args:
        gdsin (str): file name of input GDS.
        gdsout (str): optional file name of output GDS.
        whitelibrary (str): optional file name of generated white cell GDS library.
        black2whitemap (dict): mapping black cell names to white cell functions {name: function}.
        prefixb: blackbox prefix (default = 'black_')
        prefixw (str): whitebox prefix  (default = 'wb_')
        infolevel (int): amount of debug info printed (default = 0)

    Returns:
        str, dict: gds filename of white Pcell library, {black-cellname: white-cellname}
    """

    if whitelibrary is None:
        timestr = time.strftime("%Y-%m-%d")
        whitelibrary = "{}_white_lib_{}.gds".format(gdsin[:-4], timestr)

    # Using the black2whiteMap dictionary:
    # - create a list of all black cellnames in the map,
    # - extract their Pcell parameters,
    # - create a matching list of newly generated white cells using the params.
    # Note the white cell is (should be) generated with a black cell inside.
    # Note that the replaced white cells have the same name as the original
    # black cells for matching the gds instantiation reference.
    gdsinstream = nd.GDSII_stream(gdsin)
    allInputCellNames = gdsinstream.cells.keys()
    allBlackCellNames = []
    allWhiteCells = []
    allWhiteCellNames = []
    excludes = _exclude_substrings(list(black2whiteMap.keys()))
    for blackBasename, whiteFunction in black2whiteMap.items():
        #print('whiteFunction', blackBasename, whiteFunction, )
        if whiteFunction is None:
            continue
        cells_x_Params = _readPcellParams(
            blackBasename=blackBasename,
            gdsin=gdsinstream,
            allInputCellNames=allInputCellNames,
            excludes=excludes,
            infolevel=infolevel
        )
        for cellname, parameters in cells_x_Params.items():
            try:
                whiteCell = whiteFunction(**parameters)
                allWhiteCells.append(whiteCell)
                allBlackCellNames.append(cellname) #inside for loop to sync order of white and black
            except Exception as Error:
                logger.exception(
                    "White box function or whitbox function call.\n"\
                    " - can not generate cell '{}'\n"\
                    " - whitebox function: {}\n"\
                    " - parameters: {}\n".\
                    format(cellname, whiteFunction, parameters)
                )
                #traceback.print_exc(file=sys.stdout)
                raise
                #TODO: add empty/error cell to keep black and white list in sync
        #if infolevel > 3:
        #    print('\nblackBasename = ', blackBasename)
        #    print('whiteCellparams:\n', allWhiteCells)

    if infolevel > 3:
        logger.info("\n{:30}{}".format("allBlackCells", "allWhiteCells"))
        for base, cell in zip(allBlackCellNames, allWhiteCells):
            try :
                logger.info("{:40}{}".format(base, cell.cell_name))
            except:
                nd.cfg.print_except("bare-except:")
                logger.info("{:40}{}".format(base, cell))

    #generate cell names for all black and white cells
    allBlackCellNamesOut = [prefixb+name for name in allBlackCellNames]
    allWhiteCellNames = [cell.cell_name for cell in allWhiteCells]

    #map input black cellname to output black cellname (for gds debugging).
    blackmap = {black:newblack for black, newblack in zip(allBlackCellNames, allBlackCellNamesOut)}
    #map white cellname to black cellname
    whitemap = {white:black for white, black in zip(allWhiteCellNames, allBlackCellNames)}

    if infolevel > 1:
        logger.info('\nMapping for renaming black cells:')
        pprint(blackmap)
        logger.info('\nMapping to copy white cellnames to original (black) gdsin cellnames:')
        pprint(whitemap)
        #print('\n')

    #Export all white cells into the <whitelibrary> gds file
    nd.export_gds(allWhiteCells, filename=whitelibrary)

    # rename black cells:
    g = nd.GDSII_stream(whitelibrary, cellmap=blackmap)
    g.GDSII_write(whitelibrary)

    # rename white cells into original black cell names:
    #g = nd.GDSII_stream(whitelibrary, cellmap=whitemap)
    #g.GDSII_write(whitelibrary)

    #white to black cell name map:
    b2w = dict(zip(allBlackCellNames, allWhiteCellNames))
    return whitelibrary, b2w


def replaceCells(
        gdsin,
        gdsout=None,
        PcellFunctionMap=None,
        ScellMapping=None,
        md5=False,
        md5suffix='.md5',
        suffix='_white',
        addtime=True,
        infolevel=0,
        layermap=None,
        layermapmode=None,
):
    """Replace black with white cells in file <gdsin> and export result to <gdsout>.

    The replacement is Pcell (parametric cell) or a Scell (static cell) based.
    Note: Only apply one mapping per function call, Pcell *or* Scell,
    to stay in control of the mapping order. Sequential calls to
    this function are perfectly fine, but the order may matters for the outcome.

    * Pcell replacement: Needs a mapping of the Pcell basename to a cell function:
        {<black_cell_basename>, <white_cell_function_pointer>}
    * Scell replacement: Needs a gds library <ScellFile> and a map <ScellMap>:
        * ScellFile (str): filename of gds library
        * ScellMap (dict): {black_cellname: white_cellname}.

    Args:
        gdsin (str): input filename
        gdsout (str): optional output filename
        PcellMap (dict): black2white dictionary {black_cellbasename: white_Pcell_function}
        ScellMapping (dict): {<ScellFile>: {<ScellMap>}} for black2white mapping
            of one or more cell libraries.
        md5 (bool): save md5sum. Default=False
        md5suffix (str): suffix for md5sum file. Default = "_md5"
        suffix (str): suffix of the output gds. Default = "_white"
        addtime (bool): add time (date) as last suffix. Default = True
        infolevel (int): amount of debug info printed

    Returns:
        str: gds output filename of file with replaced cells
    """
    nd.cfg.gdsload = gdsin
    # avoid empty sufices.
    if md5suffix == "":
        md5suffix = ".md5"
    if suffix == "":
        suffix = "_white"

    if gdsout is None:
        if addtime:
            timestr = time.strftime("%Y-%m-%d")
            gdsout = f"{gdsin[:-4]}{suffix}_{timestr}.gds"
        else:
            gdsout = f"{gdsin[:-4]}{suffix}.gds"
    if PcellFunctionMap is not None:
        whitePcellLibname, PcellMap = (
            _createWhitePcellLib(
                gdsin,
                black2whiteMap=PcellFunctionMap,
                infolevel=infolevel,
                layermap=layermap,
                layermapmode=layermapmode,
            )
        )
        whitePcellstream = nd.GDSII_stream(whitePcellLibname)
        whitePcells = Counter(whitePcellstream.cells.keys())
        if infolevel > 0:
            logger.info("\nCreated white Pcell gds library '{}'".
                format(whitePcellLibname))
    else:
        whitePcells = Counter([])

    whiteScellstreams = {}
    ScellMap = {}     # {black_name: white_name}
    ScellstrmMap = {} # {white_name: strm_no} map the black box name to a lib stream
    if ScellMapping is not None:
        for i, (libfile, mapping) in enumerate(ScellMapping.items()):
            if infolevel > 0:
                logger.info("Processing cell library '{}'".format(libfile))
            whiteScellstreams[i] = nd.GDSII_stream(libfile)

            # validate mapping:
            whiteScellNames = whiteScellstreams[i].cells.keys()
            for whitemapname in mapping.values():
                if whitemapname not in whiteScellNames:
                    mes = "Provided white cell mapping is not consistent with"\
                        " cells in library file '{}':"\
                        " cell '{}' not found.".format(libfile, whitemapname)
                    logger.exception(mes)
                    raise Exception(mes)

            # append dictionaries of libs:
            ScellMap = {**ScellMap, **mapping}
            ScellstrmMap = {**ScellstrmMap, **{white:i for white in mapping.values()}}
            # White boxes are needed in the mappig because in the final replacement
            # the black-to-white cellmapping has been applied and the lookup
            # uses white names:
            # TODO: check there is no whitebox overlap in libs
            # TODO: check if there is remapping across libs A->B->C
        allWhiteScells = ScellstrmMap.keys()
        if infolevel > 0:
            logger.info("Detected {} white cell(s) in libraries.".
                format(len(allWhiteScells)))

    gdsinstream = nd.GDSII_stream(gdsin)
    gdsincells = Counter(gdsinstream.cells.keys())

    if len(ScellMap) > 0 and len(whitePcells) > 0:
        logger.error("Use only Pcells or Scells in a call to replaceCells().")
        return None

    # Rename pcells and scells variables for common post-processing code.
    if len(ScellstrmMap) > 0:
        whiteCellstreams = whiteScellstreams
        cellMap = ScellMap
        cellstrmMap = ScellstrmMap
    elif len(whitePcells) > 0:
        whiteCellstreams = {0: whitePcellstream}
        cellMap = PcellMap
        cellstrmMap = {white:0 for black, white in cellMap.items()}

    nd.cfg.load_gds = False

#==============================================================================
# Find which cells to replace and where needed rename
#==============================================================================
    replace = {}
    noreplace = []
    for Bname, Wname in cellMap.items():
        if Bname in gdsincells.keys():
            replace[Bname] = Wname
        else:
            noreplace.append(Bname)

    # list unmapped cellsin conflict with new cell names after map.
    stayAll = list(gdsincells - Counter(replace.keys()))
    stayNOTOK = list( (Counter(stayAll) & Counter(replace.values())).elements())
    stayOK = list(Counter(stayAll) - Counter(stayNOTOK))

    if infolevel > 0:
        logger.info(
            '{} validated replacement(s) black_gdsin_cell: white_library_cell:'
            .format(len(replace)))
        for Bname, Wname in replace.items():
            logger.info("{}: {}".format(Bname, Wname))

    # create a cellmap for gdsinstream and import the gds using the map.
    # TODO: in memomry
    cellmap = dict(replace)
    for cell in stayNOTOK:
        newname = cell+'$' #TODO: check if this one is unique too
        cellmap[cell] = newname
        stayOK.append(newname)

    #reload gdsinstream applying the cellmapping to obtain SREF with white names
    del(gdsinstream)
    gdsinstream = nd.GDSII_stream(gdsin, cellmap=cellmap)


#==============================================================================
# Replace cells
#==============================================================================
    #print('\nGenerate gdsout:')
    with open(gdsout, 'wb') as f:
        f.write(b''.join(rec.stream for rec in gdsinstream.header))

        #Walk tree for each top cell rather then linear loop.
        visited = set()
        def gdsin_iter(name, level=0):
            """Iterator over cell with <name>."""
            nonlocal stayOK
            visited.add(name)
            if name in stayOK:
                yield name, -1
            else:
                yield name, cellstrmMap[name]
                return None
            for sub in gdsinstream.cells[name].snames:
                if sub not in visited:
                    yield from gdsin_iter(sub, level+1)

        for top in gdsinstream.topcell():
            for cellname, strm_no in gdsin_iter(top):
                if infolevel > 1:
                    logger.info("strm_no, BB: {}, {}".format(strm_no, cellname))
                if strm_no == -1:
                    f.write(gdsinstream.cells[cellname].stream)
                else:
                    f.write(whiteCellstreams[strm_no].GDSII_stream_cell(cellname))

        #for cellname in whitePcells:
        #    f.write(whitePcellstream.GDSII_stream_cell(cellname))

        # for Bname, Wname in replace.items():
        #     f.write(whiteScellstream.GDSII_stream_cell(Wname))

        # for cellname in stayOK:# gdsinstream.cells.keys():
        #     f.write(gdsinstream.cells[cellname].stream)

        f.write(gdsinstream.footer.stream)

    print("Wrote white gds file '{}'".format(gdsout))

    if md5:
        md5sum = nd.util.file2md5(os.path.join(gdsout), save=False)
        logger.info("md5sum of gds output: {}".format(md5sum))
        with open(f"{gdsout[:-4]}{md5suffix}", "w") as F:
            F.write(
                "{}  {}".format(
                    md5sum, os.path.join(os.path.basename(gdsout))
                )
            )  # do not use a \n

    return gdsout
