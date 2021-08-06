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
# (c) 2016 Xaveer Leijtens
# (c) 2017 Ronald Broeke

"""Module to read a GDSII file record by record.

Substitute cell names and/or layer names.
Keep cell constructions.

The hierarchy is:

-<structure> (gds cell)
    -<element>:
        -element-type: polygon (boundary), polyline (path), annotation, sref (instantiation), etc.
        -element-property: layer, XY, etc.
        -data

Note that the cell identifier, element-type and element-property are all records.

GDSII files in nazca are a subset of the GDSII spec:

gdsii:
    HEADER          \
    BGNLIB          | this is "header" in nazca
    LIBNAME         |
    UNITS           /
    {<structure>}*  > this is "cell" in nazca
    ENDLIB          > this is "footer" in nazca

<structure>:
    BGNSTR
    STRNAME
    {<element>}*
    ENDSTR

<element>:
    {<boundary>|<path>|<sref>|<aref>|<text>|<otherstuff>}
    {data}
    ENDEL

 GDSII cells are a set of records:
 - start record
 - name record
 - various types of records, which can also be references to other cells
 - end record

 All GDSII records are structured as:
 - 16 bit int RECORD_LENGTH (number of bytes in record)
 -  8 bit int RECORD_TYPE
 -  8 bit int DATA_TYPE
 -  the data: (RECORD_LENGTH - 4) / (DATA_TYPE length) elements
"""

import struct
import io
import sys
from collections import OrderedDict, defaultdict
from . import gds_base as gb
from nazca import cfg as cfg
from nazca.mask_layers import get_layer

elm_open = {
    gb.GDS_record.BOUNDARY,
    gb.GDS_record.BOX,
    gb.GDS_record.PATH,
    gb.GDS_record.SREF,
    gb.GDS_record.AREF,
    gb.GDS_record.TEXT
}

elm_data = {
    gb.GDS_record.LAYER,
    gb.GDS_record.DATATYPE,
    gb.GDS_record.TEXTTYPE,
    gb.GDS_record.PATHTYPE,
    gb.GDS_record.WIDTH,
    gb.GDS_record.BOXTYPE,
    gb.GDS_record.XY,
    gb.GDS_record.SNAME,
    gb.GDS_record.MAG,
    gb.GDS_record.ANGLE,
    gb.GDS_record.STRING,
    gb.GDS_record.STRANS,
    gb.GDS_record.COLROW
}


class GDSII_Error(Exception):
    pass


def DFS(graph):
    """Check if a cell is a DAG.

    GDSII cells form a sorted graph that should not contain loops. This
    routine does a Depth First Search of all cells and returns an
    topologically sorted list.
    """
    def DFS_visit(cell):
        if graph.get(cell, None) is None:
            raise ValueError(f"Reference found to unavailable cell '{cell}' in gds file.")
        for c in graph[cell]:
            if c in start:
                eprint("Error in GDS\nGraph:", graph)
                raise GDSII_Error("Cyclic reference at cell '{}'.".format(c))
            if c not in parent:
                parent[c] = cell
                start.add(c)
                DFS_visit(c)
                start.remove(c)
                order.append(c)
    parent = {}
    start = set()
    order = []
    cells = graph.keys()
    for c in cells:
        if c not in parent:
            parent[c] = None
            start.add(c)
            DFS_visit(c)
            start.remove(c)
            order.append(c)
#    order.reverse() # Start with top first.
    return order

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def unpack_uint8(byte1):
    return int(struct.unpack('>B', byte1)[0])

def unpack_uint16(byte2): # unsigned
    return int(struct.unpack('>H', byte2)[0])

def unpack_int16(byte2):
    return int(struct.unpack('>h', byte2)[0])

def unpack_int32(byte4):
    return int(struct.unpack('>l', byte4)[0])

# 4-byte real is not used.


def unpack_real8(byte8):
    if byte8 == b'00000000':
        return 0.0
    # value = M / 2**56 * 16**(E - 64)
    #       = M * 16**(E - 64 - 56/4)
    #       = M * 16**(E - 78)
    E, M6, M54, M3210 = struct.unpack(">BBHL",byte8) # 1, 1, 2, 4 bytes
    if E & 0x80:
        sign = -1
    else:
        sign = 1
    E = (E & 0x7f) - 78
    M = M6 * 0x1000000000000 + M54 * 0x100000000 + M3210
    return float(sign * M * (16**E))

def nx_uint8(strm):
    f = io.BytesIO(strm)
    while True:  # Can be a stream with a number of 1-byte values
        byte1 = f.read(1)
        if byte1:
            yield unpack_uint8(byte1)
        else:
            break

def nx_int16(strm):
    f = io.BytesIO(strm)
    while True:  # Can be a stream with a number of 2-byte values
        byte2 = f.read(2)
        if byte2:
            yield unpack_int16(byte2)
        else:
            break

def nx_int32(strm):
    f = io.BytesIO(strm)
    while True:  # Can be a stream with a number of 4-byte values
        byte4 = f.read(4)
        if byte4:
            yield unpack_int32(byte4)
        else:
            break

def nx_real8(strm):
    f = io.BytesIO(strm)
    while True:  # Can be a stream with a number of 4-byte values
        byte8 = f.read(8)
        if byte8:
            yield unpack_real8(byte8)
        else:
            break


class GDSII_cell:
    """Class for storing GDS cell content.

    A cell contains three attributes that constitute all cell content
        -header stream
        -list of element objects (sref, polygon, polyline, etc)
        -set of cell references <snames>
        -footer stream
    """

    def __init__(self, pos=None):
        """Initialize a GDSII_cell.

        Returns:
            None
        """
        self.snames = set()
        self.count_srefs = defaultdict(int)  # keep track of number of references to cells
        self.header = []
        self.elements = []
        self.footer = GDSII_record(gb.gds_endstr())
        self.pos = pos

    def __str__(self):
        """String representation of the cell."""

        return str('\n'+'\n'.join(str(h) for h in self.header) + '\n' +
            '\n'.join(str(elem) for elem in self.elements) + '\n' +
            str(self.footer))

    def references(self, sname):
        """Add a cell reference <sname> to the cell.

        Args:
            sname (str): name of referenced cell

        Returns:
            None
        """
        self.snames.add(sname)
        self.count_srefs[sname] += 1

    def addelem(self, elem):
        """Add an GDSII-element to the cell.

        Args:
            elem (GDSII_element):

        Returns:
            None
        """
        self.elements.append(elem)

    @property
    def name(self):
        """Get cell name.

        Returns:
            str: name of the cell
        """
        return self.header[1].data

    @property
    def stream(self):
        """Create the full stream of the cell.

        Concatenate the list of streams of headers, elements and the footer.

        Returns:
            bytearray: stream of the cell
        """
        return b''.join(h.stream for h in self.header) + \
               b''.join(e.stream for e in self.elements) + \
               self.footer.stream


class GDSII_element:
    """Class for GDS elements.

    A GDSii_element object stores a list of GDSII_records and a number of
    methods to read out the properties (layer and xy position) of each record.

    The byte-content of the element is a number of records.
    The content start and end are identified by a record, i.e.
    they are the first (element type) and last record (endel) of the element.

    Element types are polylines, polygons, sref, aref, etc.
    Element properties are identifiers like layer, XY, etc., followed by data
    of the property.
    """

    def __init__(self, records=None):
        """Initialize a GDS element.

        Args:
            records (list of GDSII_record): records to store in the element
            (default = None)

        Returns:
            None
        """
        if records is None:
            records = []
        self.records = list(records)  # force copy

    def __str__(self):
        return '\n'.join(str(rec) for rec in self.records)

    def addrecord(self, record):
        """Add a GDSII-record to the GDSII_element.

        Args:
            record (GDSII_record): record to add

        Returns:
            None
        """
        self.records.append(record)

    @property
    def stream(self):
        """Create a stream of all records in the GDSII_element.

        Returns:
            bytearray: gds stream
        """
        return b''.join(rec.stream for rec in self.records)

    @property
    def etype(self):
        """Return element type.

        The element type is the first record in the 'records' list.
        Element types are in set 'elm_open'.
        """
        return self.records[0].rtype

    @property
    def annotation(self):
        """Get the properties of an annotation element.

        Returns None if the element is not an annotation.

        Returns:
            [int, tuple, str]: layer, position (x, y), text
        """
        # GDS: TEXT [ELFLAGS] [PLEX] LAYER <textbody>
        if self.etype != gb.GDS_record.TEXT:
            return None
        layer, texttype, XY, text = None, None, None, None
        for r in self.records[1:]:
            if r.rtype == gb.GDS_record.LAYER:
                layer = r.data[0]
            elif r.rtype == gb.GDS_record.TEXTTYPE:
                texttype = r.data[0]
            elif r.rtype == gb.GDS_record.XY:
                XY = r.data
            elif r.rtype == gb.GDS_record.STRING:
                text = r.data
        # annotations have a TEXTTYPE not a DATATYPE
        return [(layer, texttype), XY, text]

    @property
    def polyline(self):
        """Get the properties of a polyline.

        Returns None if the element is not a polyline.

        Returns:
            [int, tuple]: layer, position (x, y)
        """
        # GDS: PATH [ELFLAGS] [PLEX] LAYER DATATYPE [PATHTYPE][WIDTH] XY
        if self.etype != gb.GDS_record.PATH:
            return None
        layer, XY, pathtype, width = None, None, None, None
        for r in self.records[1:]:
            if r.rtype == gb.GDS_record.LAYER:
                layer = r.data[0]
            elif r.rtype == gb.GDS_record.DATATYPE:
                datatype = r.data[0]
            elif r.rtype == gb.GDS_record.PATHTYPE:
                pathtype = r.data[0]
            elif r.rtype == gb.GDS_record.WIDTH:
                width = r.data[0]
            elif r.rtype == gb.GDS_record.XY:
                XY = r.data
        return [(layer, datatype), XY, pathtype, width]

    @property
    def polygon(self):
        """Get the properties of a polygon.

        Returns None if the element is not a polygon.

        Returns:
            [int, tuple]: layer, position (x, y)
        """
        # GDS: BOUNDARY [ELFLAGS] [PLEX] LAYER DATATYPE XY
        if self.etype != gb.GDS_record.BOUNDARY:
            return None
        layer, XY = None, None
        for r in self.records[1:]:
            if r.rtype == gb.GDS_record.LAYER:
                layer = r.data[0]
            elif r.rtype == gb.GDS_record.DATATYPE:
                datatype = r.data[0]
            elif r.rtype == gb.GDS_record.XY:
                XY = r.data
        return [(layer, datatype), XY]

    @property
    def box(self):
        """Get the properties of a box element.

        Returns None if the element is not an annotation.

        Returns:
            [int, tuple, str]: layer, position (x, y)
        """
        # GDS: BOX [ELFLAGS] [PLEX] LAYER BOXTYPE XY
        if self.etype != gb.GDS_record.TEXT:
            return None
        layer, boxtype, XY = None, None, None
        for r in self.records[1:]:
            if r.rtype == gb.GDS_record.LAYER:
                layer = r.data[0]
            elif r.rtype == gb.GDS_record.BOXTYPE:
                boxtype = r.data[0]
            elif r.rtype == gb.GDS_record.XY:
                XY = r.data
        # box has a BOXTYPE, not a DATATYPE
        return [(layer, 0), XY]

    @property
    def instance(self):
        """Get the properties of a SREF element.

        Returns:
            [int, tuple, str]: layer, position (x, y), text
        """
        # GDS: SREF [ELFLAGS] [PLEX] SNAME [<strans>] XY
        #     <strans: STRANS [MAG] [ANGLE]
        if self.etype != gb.GDS_record.SREF:
           return None
        sname, strans, mag, angle, XY = None, 0, 1.0, 0.0, None
        for r in self.records[1:]:
            if r.rtype == gb.GDS_record.SNAME:
                sname= r.data
            elif r.rtype == gb.GDS_record.STRANS:
                if r.data[0] & 0x8000:
                    strans = 1
            elif r.rtype == gb.GDS_record.MAG:
                mag = r.data[0]
            elif r.rtype == gb.GDS_record.ANGLE:
                angle = r.data[0]
            elif r.rtype == gb.GDS_record.XY:
                XY = r.data
        return [sname, strans, mag, angle, XY]

    @property
    def array(self):
        """Get the properties of an AREF element.

        Returns:
            [int, tuple, str]: layer, position (x, y), text
        """
        # GDS: SREF [ELFLAGS] [PLEX] SNAME [<strans>] XY
        #     <strans: STRANS [MAG] [ANGLE]
        if self.etype != gb.GDS_record.AREF:
            return None
        sname, strans, mag, angle, col, row, XY = None, None, 1.0, None, None, None, None
        for r in self.records[1:]:
            if r.rtype == gb.GDS_record.SNAME:
                sname= r.data
            elif r.rtype == gb.GDS_record.STRANS:
                strans = r.data
            elif r.rtype == gb.GDS_record.MAG:
                mag = r.data[0]
            elif r.rtype == gb.GDS_record.ANGLE:
                angle = r.data[0]
            elif r.rtype == gb.GDS_record.COLROW:
                col, row = r.data
            elif r.rtype == gb.GDS_record.XY:
                XY = r.data
        return [sname, strans, mag, angle, col, row, XY]


class GDSII_record:
    """Class for storing a GDSii record in byte stream format.

    Note that cells and elements are (build from) records.
    Records have the following structure:

    byte 0, 1: record byte length

    byte 2, 3: record-type (GDS-records)

    byte 4, 5: data-type (GDS_datatype). NODATA when there is no data.

    byte 6, rlen: data, or nothing when data-type is NODATA
    """

    def __init__(self, strm, pos=0):
        """Construct a GDSII_record.

        Args:
            strm (bytearray): stream
            pos (int): start of record in the stream

        Returns:
            None
        """
        self.strm = strm  # byte array
        self.pos = pos  # Start of record

    def __str__(self):
        """String representation of the record."""
        if self.rtype in elm_open:
            e = '\n┌'
            d = '│'
        elif self.rtype in elm_data:
            e = '├'
            d = '│'
        elif self.rtype == gb.GDS_record.ENDEL:
            e = '└'
        else:
            e = ''
            d = ''

        s = '{}{}, {}, 4+{} bytes'.format(
            e,
            gb.GDS_record.name[self.rtype],
            gb.GDS_datatype.name[self.dtype],
            self.rlen - 4)
        if self.dtype == gb.GDS_datatype.ASCII:
            s += '\n{} "{}"'.format(d, self.data)
        elif self.dtype != gb.GDS_datatype.NODATA:
            s += '\n{} {}'.format(d, self.data)

        return s

    @property
    def rlen(self):
        """Return record length.

        GDSII spec defines length as signed 2-byte int, but the sign does not
        make sense and there are GDSII files with large (> 0x8000) record
        lengths. By using unsigned here (uint16), we can also read those files.
        """
        return unpack_uint16(self.strm[self.pos:self.pos+2])

    @property
    def rtype(self):
        """get GDS-record-type

        Returns:
            uint8: record-type value
        """
        return unpack_uint8(self.strm[self.pos+2:self.pos+3])

    @property
    def dtype(self):
        """Get GDS data-type.

        Returns:
            uint8: data-type value
        """
        return unpack_uint8(self.strm[self.pos+3:self.pos+4])

    @property
    def data(self):
        """Convert the gds stream into proper data.

        Returns:
            None
        """
        if self.dtype == gb.GDS_datatype.NODATA:
            return None
        elif self.dtype == gb.GDS_datatype.BITARRAY:
            return list(nx_int16(self.strm[self.pos+4:self.pos+6]))
        elif self.dtype == gb.GDS_datatype.INT16:
            return list(nx_int16(self.strm[self.pos+4:self.pos+self.rlen]))
        elif self.dtype == gb.GDS_datatype.INT32:
            return list(nx_int32(self.strm[self.pos+4:self.pos+self.rlen]))
        elif self.dtype == gb.GDS_datatype.REAL8:
            return list(nx_real8(self.strm[self.pos+4:self.pos+self.rlen]))
        elif self.dtype == gb.GDS_datatype.ASCII:
            if self.strm[self.pos+self.rlen-1]:  # Remove last 0 byte padding
                e = self.pos + self.rlen
            else:
                e = self.pos + self.rlen - 1
            return self.strm[self.pos+4:e].decode('utf-8')
        else:
            raise GDSII_Error('Unknown data type {}'.format(self.dtype))

    @property
    def stream(self):
        """Return the proper bytestream for this record.

        Return:
            bytearray: record
        """
        return self.strm[self.pos:self.pos + self.rlen]


class GDSII_stream:
    """Class to read, modify and write a GDSII stream.

    A stream consists of a header, cells with elements and a footer.
    """

    def __init__(
        self,
        filename,
        cellmap=None,
        layermap=None,
        layermapmode=None,
        parse=True,
    ):
        """Initiliaze a GDSII_stream.

        Args:
            filename (str|bytes): file name to read or buffer with contents
                of gds file, in bytes or bytesarray.
            cellmap (dict): {celname_old:cellname_new} optional dictionary to
                map cell names in <filename> onto new names in the stream
            layermap (dict): {layer_number_old:layer_number_new} optional
                dictionary to map layer numbers in <filename> onto new layer
                numbers in the stream

        Returns:
            None
        """
        if cellmap is None:
            cellmap = {}
        if layermap is None:
            layermap = {}
        if cfg.layermap:  # add the global layer as a start
            layermap = {**cfg.layermap, **layermap}

        self.header = []
        self.cells = OrderedDict()  # All cells: {cellname: GDSII_cell}
        self.footer = GDSII_record(gb.gds_endlib())

        # Hold the GDSII file structure
        self.snames = set()  # Names of all referenced cells.
        self.graph = dict()  # Holds the DAG structure of the gds.
        self.ref_del = set()  # Names of cells referenced from deleted cells.
        self.layer_remove = set()  # Names of layers to be removed.
        self.cell_remove = set()  # Names of cells to be removed.
        self.layermap = layermap
        self.layermapmode = layermapmode
        self.cellmap = {}
        self.gds_db_user = None
        self.gds_db_unit = None

        if isinstance(filename, str):
            self.filename = filename
            with open(filename, 'rb') as f:
                self.gds_stream = bytearray(f.read())
            # eprint("Read GDS from '{}'".format(filename))
        else:
            self.filename = '<buffer>'
            self.gds_stream = filename
            # eprint("Read GDS from buffer")
        self.stream_length = len(self.gds_stream)  # stream length in bytes.
        self.create_layermap(layermap, layermapmode)
        self.create_cellmap(cellmap)
        if parse:
            self.parse()
            for c in self.cells:
                # gather all referenced cells.
                self.snames |= self.cells[c].snames
                # Build the graph
                self.graph[c] = self.cells[c].snames
            # This will check for cyclic references as well
            self.order = DFS(self.graph)


    def create_layermap(self, layermap, layermapmode):
        """Create a valid internal layermap from the provided <layermap>.

        The internal layermap is stored in dict self.layermap

        Args:
            layermap (dict:)

        Returns:
            None
        """
        layermap_copy = layermap.copy()
        for _layi, _layo in layermap.items():
            layi, layo = _layi, _layo
            if isinstance(layi, int):
                layi = (layi, 0)
            if isinstance(layo, int):
                layo = (layo, 0)
            if _layi is None or _layo is None:
                self.layer_remove.add(layi)
            else:
                assert isinstance(layi, tuple)
                assert isinstance(layo, (tuple, str))
                layermap_copy[layi] = layo
        self.layermap = layermap_copy
        return None


    def create_cellmap(self, cellmap):
        """Create a valid internal cellmap from the provided <cellmap>.

        The internal cellmap is stored in dict self.cellmap.

        Args:
            cellmap (dict:)

        Returns:
            None
        """
        for cell in cellmap:
            if cellmap[cell] is None:
                self.cell_remove.add(cell)
            else:
                self.cellmap[cell] = cellmap[cell]
        return None


    def count_srefs(self):
        """Print how many times cells have been instantiated to stdout."""
        for cellname, cell in sorted(self.cells.items()):
            print("{}:".format(cellname))
            for sref, count in sorted(cell.count_srefs.items()):
                print("   {:4d}x {}".format(count, sref))


    def GDSII_write(self, filename):
        """Write a GDSII stream to file.

        Args:
            filename (str): output filename for binary gds

        Return:
            None
        """
        with open(filename, 'wb') as f:
            f.write(b''.join(rec.stream for rec in self.header))
            if self.cell_remove == set(): # No cells were removed
                for cellname in self.cells:
                    f.write(self.cells[cellname].stream)
            else:
                for cellname in self.topcell():
                    # Write the topcell trees
                    self.GDSII_write_cell_and_under(f, cellname)
            f.write(self.footer.stream)


    def ASCII_write(self, filename=False):
        """Write the GDS in a human readable format.

        If no filename is given (False) output is directed to stdout.

        Args:
            filename (str): output filename for ascii representation of gds
                (default = False)

        Returns:
            str: gds in ascii representation
        """
        buf = io.StringIO()
        for rec in self.header:
            buf.write(str(rec)+'\n')
        for cellname in self.cells:
            if cellname not in self.cell_remove:
                buf.write(str(self.cells[cellname])+'\n')
        buf.write(str(self.footer)+'\n')
        if filename is False: # Write to stdout
            sys.stdout.write(buf.getvalue())
        elif filename is None: # Don't write, only return the string.
            pass
        else:
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(buf.getvalue())
        return buf.getvalue()


    def _GDSII_stream_cell_and_under(self, strm, cellname, init=True):
        """Obtain a GDS stream of cell named <cellname>, and its children.

        Iterative part of method CDSII_stream_cell.

        Args:
            strm (bytearrray): gds stream being build, starting as bytearray()
            cellname (str): cell name of top cell
            init (bool): flag to check is a top cell level (default:True)

        Returns:
            bytearray: gds stream
        """
        global done
        if init:
            done = []
        if cellname not in done:
            try:
                strm.extend(self.cells[cellname].stream)
            except:
                raise Exception("Error: Looking up a non-existing cellname '{}' "\
                    "in file '{}'. Valid cellnames are {}."\
                    .format(cellname, self.filename, list(self.cells.keys())))
            done.append(cellname)

        for subcellname in self.cells[cellname].snames:
            self._GDSII_stream_cell_and_under(strm, subcellname, init=False)
        return strm


    def GDSII_stream_cell(self, cellname):
        """Return the GDSII stream of the cell with name <cellname> and below.

        Args:
             cellname (str): name of the cell to get the gds stream off

        Returns:
            bytearray: gds stream
        """
        return self._GDSII_stream_cell_and_under(bytearray(), cellname,
            init=True)


    def _GDSII_write_cell_and_under(self, f, cellname, init=True):
        """Write gds stream of cell named <cellname> to file <f>.

        Iterative part of method GDSII_write_cell.

        Args:
            f (file): file handle of file to write to
            cellname (str): name of the cell to write to file
            init (bool): flag to check if in top cell

        Returns:
            None
        """
        global done
        if init:
            done = []
        if cellname in self.cell_remove:
            return
        if cellname not in done:
            f.write(self.cells[cellname].stream)
            done.append(cellname)
        for subcellname in self.cells[cellname].snames:
            self._GDSII_write_cell_and_under(f, subcellname, False)


    def GDSII_write_cell(self, cellname, filename):
        """Write a cell (and the cells referenced by it).

        Args:
            cellname (str): name of the cell to write
            filename (str): name of the file the gds stream is saved in

        Returns:
            None
        """
        with open(filename, 'wb') as f:
            f.write(b''.join(rec.stream for rec in self.header))
            self.GDSII_write_cell_and_under(f, cellname)
            f.write(self.footer.stream)


    def topcell(self):
        """Get all topcells in the stream.

        Returns:
            list of str: list of all cell names that are top cells (are not referenced)
        """
        # This does not change if cells are deleted, which is used to write
        # out the GDSII with deep-deleted cells.
        cells = set(self.cells.keys())
        return cells - self.snames


    def cell_branch(self, cellname, cellnames=None, level=0):
        """Create a set of cellnames of all cells in branch <cellname>.

        Returns:
            set: cellnames
        """
        #TODO: check if is this the same function(ality) as GDSII_stream_cell?
        if level == 0:
            cellnames = set()
        cellnames.add(cellname)

        if cellname not in list(self.cells.keys()):
            raise Exception("Error: Looking up a non-existing cellname '{}' "\
                "in file '{}'. Valid cellnames are {}.".\
                    format(cellname, self.filename, list(self.cells.keys())))
        for subname in self.cells[cellname].snames:
            self.cell_branch(subname, cellnames, level+1)
        level -= 1
        if level == -1:
            return cellnames


    def _print_structure(self, name, level=0, sort=False):
        """Print the cell tree in ascii format.

        Args:
            name (str): cellname
            level (int): function internal recursive counter (default = 0)
            sort (bool): sort cellnames alphabetically

        Returns:
            None
        """
        if level == 0:
            print("□", name)
        else:
            print("  {}└{}".format("│ " * (level-1), name))
            #TODO remove one '|' if last cell in a branch starts a deeper level.
        if sort:
            snames = sorted(self.cells[name].snames)
        else:
            snames = self.cells[name].snames
        for sub in snames:
            self._print_structure(sub, level+1)


    def print_structure(self, name=None, sort=False):
        """Print the cell tree in ascii format.

        Args:
            name (str | list of str): cellname(s)
            sort (bool): sort cellnames alphabetically

        Returns:
            None
        """
        if name is None:
            names = self.topcell()
        elif isinstance(name, str):
            names = [name]
        for i, name in enumerate(names):
            if len(names) > 1:
                print('topcell-{}:'.format(i))
            self._print_structure(name, sort)


    def _gds_cell_iter(self, rec_iter):
        """Generator to iterate over all cells with cell and layer mapping applied.

        The rec_iter is the thread along which all cells and elements are
        found in a GDS byte-stream. The stream is split and stored in
        internal cell and element dictionaries.

        This method takes care of cell renaming.
        For each cell it calls method 'gds_elem_iter' which splits
        and stores all elements per the cell while applying cell-reference
        mapping and layer, datatype mapping.

        Removing a cell (mapping to None) is done by removing all
        references to that cell. The cell itself remains behind, making it
        a topcell (= unreferenced cell). When writing the gds file this
        topcell is not written. Instances down-stream the cell tree are
        also discarded, unless referenced in other topcell trees.

        Args:
            rec_iter (GDSII_record iterator): scan through record for cells

        Yields:
            GDSII_cell: next cell object in rec_iter after mapping
        """
        for rec, pos in rec_iter:
            if rec.rtype == gb.GDS_record.ENDLIB: # End of file.
                break
            cell = GDSII_cell(pos)

            # Read BGNSTR
            if rec.rtype != gb.GDS_record.BGNSTR:
                raise GDSII_Error("'BGNSTR' record expected. Got '{}'.".format(
                    gb.GDS_record.name[rec.rtype]))
            cell.header.append(rec)

            # Read STRNAME
            rec, pos = next(rec_iter)
            if rec.rtype != gb.GDS_record.STRNAME:
                raise GDSII_Error("STRNAME record expected")
            cellname = rec.data

            if cellname in self.cellmap:
                if self.cellmap[cellname] is not None:
                    cellname = self.cellmap[cellname]
                    # Substitute new cellname: make new STRNAME record
                    rec = GDSII_record(gb.gds_strname(cellname))
            cell.header.append(rec)
            # The optional STRCLASS record is assumed to be absent.

            # Loop over this cell's elements
            for elem in self._gds_elem_iter(rec_iter, cell):
                cell.addelem(elem)
            yield cell


    def _gds_elem_iter(self, rec_iter, cell):
        """Generator over all GDS elements in <rec_iter> within <cell>

        If a layermap or cellmap has been provided in the load,
        they will be applied here.

        Record iterator rec_iter is the thread along which GDS elements are
        rebuild and stored internally in an element dictionary.
        This Generator provides an element iterator and maps layer, datatype
        and cell reference names as specified in class attributes
        'layermap' and 'cellmap', before yielding.


        Args:
            rec_iter (GDSII_record iterator): records to scan through for elements.
                rec_iter starts off just beyond the cell header and stops at
                the cell ending.
            cell (GDSII_cell): active cell being copied.
                Needed for adding instances to the Nazca cell reference list.

        Yields:
            GDSII_element: next element object in rec_iter after mapping
        """
        elem = GDSII_element()
        for rec, pos in rec_iter:
            if rec.rtype == gb.GDS_record.LAYER:
                lay = rec.data[0]
                continue

            # DATATYPE after LAYER: boundary, path
            # TEXTTYPE after LAYER: annotation
            # BOXTYPE after LAYER: box
            elif (rec.rtype == gb.GDS_record.DATATYPE) or \
                 (rec.rtype == gb.GDS_record.TEXTTYPE) or \
                 (rec.rtype == gb.GDS_record.BOXTYPE):
                dtype = rec.data[0]
                layID = (lay, dtype)
                if (
                    layID in self.layer_remove
                    or (layID not in self.layermap and self.layermapmode == 'none')
                ):
                    # Remove this element: read until ENDEL and delete/recreate elem
                    for rec2, pos in rec_iter:
                        if rec2.rtype == gb.GDS_record.ENDEL:
                            elem = GDSII_element()
                            break
                    continue
                if layID in self.layermap: # Replace with new layer number.
                    LDT = self.layermap[layID]
                    if cfg.gds_over_getlayer:
                        LDT = cfg.layername2LDT[get_layer(LDT)][0:2]
                    rec1 = GDSII_record(gb.gds_layer(LDT[0]))
                    dat1 = LDT[1]
                else:
                    rec1 = GDSII_record(gb.gds_layer(lay))
                    dat1 = dtype
                if rec.rtype == gb.GDS_record.DATATYPE:
                    rec2 = GDSII_record(gb.gds_datatype(dat1))
                elif rec.rtype == gb.GDS_record.TEXTTYPE:
                    rec2 = GDSII_record(gb.gds_texttype(dat1))
                elif rec.rtype == gb.GDS_record.BOXTYPE:
                    rec2 = GDSII_record(gb.gds_boxtype(dat1))
                elem.addrecord(rec1)
                elem.addrecord(rec2)
                lay = -1 #reset layer for next loop
                continue

            #instance name mapping:
            elif rec.rtype == gb.GDS_record.SNAME: # Reference to other cell.
                if rec.data in self.cell_remove: # Remove reference
                    cell.references(rec.data) # But still record the reference.
                    # Remove this element: read until ENDEL and delete.
                    for rec2, pos in rec_iter:
                        if rec2.rtype == gb.GDS_record.ENDEL:
                            elem = GDSII_element()
                            break
                    continue
                if rec.data in self.cellmap:
                    rec = GDSII_record(gb.gds_sname(self.cellmap[rec.data]))
                cell.references(rec.data)

            elif rec.rtype == gb.GDS_record.ENDEL: # End of element.
                elem.addrecord(rec)
                yield elem
                elem = GDSII_element() # New element
                continue

            elif rec.rtype == gb.GDS_record.ENDSTR: # Last element (cell footer)
                break # End of elements for this cell.

            elem.addrecord(rec) #transparent record copy


    def gds_record_iter(self, strm, strmlen, pos=0):
        """Generator over the cell records.

        Args:
            strm (bytearray): GDS stream
            strmlen (int): length of the stream to iterate over
            pos (int): start position in the stream (default = 0)

        Yields:
            GDSII-record: next record object, be it cell, element or element-property
        """
        while pos < strmlen:
            rec = GDSII_record(strm, pos)
            yield rec, pos
            pos += rec.rlen


    def parse(self):
        """Find all the cells in a GDS and add them to an internal dictionary.

        This method uses a record iterator to scan its GDS stream and
        construct an internal dictionary of GDS cells and elements with
        cell mapping and layer, datatype mapping applied.

        Returns:
            None
        """
        rec_iter = self.gds_record_iter(self.gds_stream, self.stream_length)

        # Read the header: all records until (including) UNITS
        for rec, pos in rec_iter:
            self.header.append(rec)
            if rec.rtype == gb.GDS_record.UNITS:
                self.gds_user = rec.data[0]
                self.gds_db_unit = rec.data[1]
                break # last record of header

        # Read all the cells, ends on ENDLIB record (end of file).
        for cell in self._gds_cell_iter(rec_iter):
            # Add cell reference
            self.cells[cell.name] = cell


    @property
    def libname(self):
        for rec in self.header:
            if rec.rtype == gb.GDS_record.LIBNAME:
                return str(rec.data)
        else:
            return None # No LIBNAME record in GDS file header


    @property
    def gdsversion(self):
        # First header record should be HEADER with version number
        if self.header[0].rtype == gb.GDS_record.HEADER:
            return self.header[0].data[0]
        else:
            return None
