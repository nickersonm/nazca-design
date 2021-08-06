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
# Low-level GDSII routines
#
# (c) 2011 Xaveer Leijtens
#
"""Low level GDSII routines

Gds_base is a module that holds routines which implement some of the basic
GDSII functions. Currently it is limited to the following GDSII data types:
- polyline (path)
- polygon (boundary)
- string (annotation)
"""
import struct
import array
from numpy import int64
import datetime
from itertools import groupby
from math import floor
from nazca.logging import logger


# Size of databse unit (DBU) in meters (nanometer)
# Multiply the gds integer with the db_unit to get the real size.
gds_db_unit = 0.000000001

# Size of DB unit in user units (gds_db_unit / 1 um: 0.001 micron)
gds_db_user = 0.001  
# In Nazca the ratio db_unit / db_user has to remain 1e-6 because it has an internal 
# of unit 1 um.
# The gds_db_user is not really important as DBU already sets an absolute scale.
#
# Example:
#   if you go to a 5 nm grid use:
#   gds_db_unit = 0.000000005
#   gds_db_user = 0.005


# Maximum and minimum coordinates that may still display in Klayout.
# For default settings this corresponds to +/- 1 meter.
# This may give visual feedback, so users can figure out what goes wrong.
maxint32 = struct.pack(">l", 1000000000)
minint32 = struct.pack(">l", -1000000000)


# GDSII record types
# fmt: off
class GDS_record:
    """Define gds record names."""
    (
        HEADER, BGNLIB, LIBNAME, UNITS, ENDLIB, BGNSTR, STRNAME, ENDSTR, BOUNDARY, PATH,
        SREF, AREF, TEXT, LAYER, DATATYPE, WIDTH, XY, ENDEL, SNAME, COLROW, TEXTNODE,
        NODE, TEXTTYPE, PRESENTATION, SPACING, STRING, STRANS, MAG, ANGLE, UINTEGER,
        USTRING, REFLIBS, FONTS, PATHTYPE, GENERATIONS, ATTRTABLE, STYPTABLE, STRTYPE,
        ELFLAGS, ELKEY, LINKTYPE, LINKKEYS, NODETYPE, PROPATTR, PROPVALUE, BOX, BOXTYPE,
        PLEX, BGNEXTN, ENDEXTN, TAPENUM, TAPECODE, STRCLASS, RESERVED, FORMAT, MASK,
        ENDMASKS, LIBDIRSIZE, SRFNAME, LIBSECUR, BORDER, SOFTFENCE, HARDFENCE, SOFTWIRE,
        HARDWIRE, PATHPORT, NODEPORT, USERCONSTRAINT, SPACER_ERROR, CONTACT,
    ) = list(range(70))
    name = (
        "header", "bgnlib", "libname", "units", "endlib", "bgnstr", "strname", "endstr",
        "boundary", "path", "sref", "aref", "text", "layer", "datatype", "width", "xy",
        "endel", "sname", "colrow", "textnode", "node", "texttype", "presentation",
        "spacing", "string", "strans", "mag", "angle", "uinteger", "ustring", "reflibs",
        "fonts", "pathtype", "generations", "attrtable", "styptable", "strtype",
        "elflags", "elkey", "linktype", "linkkeys", "nodetype", "propattr", "propvalue",
        "box", "boxtype", "plex", "bgnextn", "endextn", "tapenum", "tapecode",
        "strclass", "reserved", "format", "mask", "endmasks", "libdirsize", "srfname",
        "libsecur", "border", "softfence", "hardfence", "softwire", "hardwire",
        "pathport", "nodeport", "userconstraint", "spacer_error", "contact",
    )
# fmt: on


# GDSII data values
class GDS_datatype:
    """Define gds record values."""
    NODATA, BITARRAY, INT16, INT32, REAL4, REAL8, ASCII = list(range(7))
    name = ("nodata", "bitarray", "int16", "int32", "real4", "real8", "ascii")


def round_to_db_unit(x):
    return floor(x / gds_db_user + 0.5)


# .pack will check for out of range values.
def pack_uint8(B):  # unsigned
    return struct.pack(">B", B)


def pack_int16(h):  # signed
    return struct.pack(">h", h)


def pack_uint16(H):  # unsigned
    return struct.pack(">H", H)


def pack_int32(l):  # signed
    try:  # Use try clause to not cause any overhead.
        return struct.pack(">l", l)
    except struct.error:
        logger.error(
            f"GDS int32 overflow: {l}. Coordinate out of range.\n"
            "  Number should be -2147483648 <= number <= 2147483647\n"
            f"  You try to write structures larger than {l*gds_db_unit:g} meter!\n"
            "  GDS output is INVALID!"
        )
        return maxint32 if l > 0 else minint32


def pack_real8(D):
    assert -7e75 < D < 7e75, "Floating point out of range"
    rec = [0] * 8
    if D < 0:
        rec[0] = 0x80
        D = -D
    E = 64
    if D < 8e-78:  # 16^-64
        rec[0] |= 64
        return array.array("B", rec).tobytes()
    # Mantissa in [1/16,1>
    if D >= 1:
        while D >= 1:
            D /= 16
            E += 1
    elif D < 0.0625:
        while D < 0.0625:
            D *= 16
            E -= 1
    rec[0] |= E & 0x7F
    M = int64(D * 0x100000000000000 + 0.5)
    for i in range(7, 0, -1):
        rec[i] = M & 0xFF
        M >>= 8
    return array.array("B", rec).tobytes()


def pack_datime(lcltm=None):  # Local date and time
    if lcltm is None:
        lcltm = datetime.datetime.today()
    return (
        pack_int16(lcltm.year)
        + pack_int16(lcltm.month)
        + pack_int16(lcltm.day)
        + pack_int16(lcltm.hour)
        + pack_int16(lcltm.minute)
        + pack_int16(lcltm.second)
    )


# If needed, padd NULL to an even byte count
def pack_padstring(bytestring):
    if len(bytestring) & 1:
        return bytestring + b"\0"
    else:
        return bytestring


# GDS records
def gds_header():  # Record length = 4 + 2
    return (
        pack_int16(6)
        + pack_uint8(GDS_record.HEADER)
        + pack_uint8(GDS_datatype.INT16)
        + pack_int16(0x0258)
    )  # software version


def gds_bgnlib(lcltm=None):  # Record length = 4 + 2*12
    return (
        pack_int16(28)
        + pack_uint8(GDS_record.BGNLIB)
        + pack_uint8(GDS_datatype.INT16)
        + pack_datime(lcltm)
        + pack_datime(lcltm)
    )  # Creation and modification dates


def gds_libname(string):  # Record length = 4 + strlen
    name = bytes(string, encoding="UTF-8")
    l = len(name)
    if l & 1:
        l += 1
    return (
        pack_int16(4 + l)
        + pack_uint8(GDS_record.LIBNAME)
        + pack_uint8(GDS_datatype.ASCII)
        + pack_padstring(name)
    )


def gds_units():  # Record length = 4 + 2*8
    return (
        pack_int16(20)
        + pack_uint8(GDS_record.UNITS)
        + pack_uint8(GDS_datatype.REAL8)
        + pack_real8(gds_db_user)
        + pack_real8(gds_db_unit)
    )


def gds_nodata(rec):  # Record length = 4
    return pack_int16(4) + pack_uint8(rec) + pack_uint8(GDS_datatype.NODATA)


def gds_bgnstr(lcltm=None):  # Record length = 4 + 2*12
    return (
        pack_int16(28)
        + pack_uint8(GDS_record.BGNSTR)
        + pack_uint8(GDS_datatype.INT16)
        + pack_datime(lcltm)
        + pack_datime(lcltm)
    )  # Creation and modification dates


def gds_sref():
    return gds_nodata(GDS_record.SREF)


def gds_aref():
    return gds_nodata(GDS_record.AREF)


def gds_endstr():
    return gds_nodata(GDS_record.ENDSTR)


def gds_path():
    return gds_nodata(GDS_record.PATH)


def gds_boundary():
    return gds_nodata(GDS_record.BOUNDARY)


def gds_endel():
    return gds_nodata(GDS_record.ENDEL)


def gds_endlib():
    return gds_nodata(GDS_record.ENDLIB)


def gds_text():
    return gds_nodata(GDS_record.TEXT)


def gds_strname(string):  # Record length = 4 + strlen
    name = bytearray(string, encoding="UTF-8")
    l = len(name)
    if l & 1:
        l += 1
    return (
        pack_int16(4 + l)
        + pack_uint8(GDS_record.STRNAME)
        + pack_uint8(GDS_datatype.ASCII)
        + pack_padstring(name)
    )


def gds_string(string):  # Record length = 4 + strlen
    name = bytes(string, encoding="UTF-8")
    l = len(name)
    if l & 1:
        l += 1
    return (
        pack_int16(4 + l)
        + pack_uint8(GDS_record.STRING)
        + pack_uint8(GDS_datatype.ASCII)
        + pack_padstring(name)
    )


def gds_layer(lay):  # Record length = 4 + 2
    return (
        pack_int16(6)
        + pack_uint8(GDS_record.LAYER)
        + pack_uint8(GDS_datatype.INT16)
        + pack_uint16(lay)
    )


def gds_datatype(datatype):  # Record length = 4 + 2
    return (
        pack_int16(6)
        + pack_uint8(GDS_record.DATATYPE)
        + pack_uint8(GDS_datatype.INT16)
        + pack_int16(datatype)
    )


def gds_width(w):  # Record length = 4 + 4
    return (
        pack_int16(8)
        + pack_uint8(GDS_record.WIDTH)
        + pack_uint8(GDS_datatype.INT32)
        + pack_int32(round_to_db_unit(w))
    )


def gds_colrow(col, row):  # Record length = 4 + 4
    return (
        pack_int16(8)
        + pack_uint8(GDS_record.COLROW)
        + pack_uint8(GDS_datatype.INT16)
        + pack_uint16(col)
        + pack_uint16(row)
    )


def gds_xy(XY, close, min_length=1, unique=True):
    """Create xy record.

    Args:
        XY (list (float, float)): list of points (x, y)
        close:
        min_length (int): default=1
        unique (bool): remove consecutive duplicates (default = True).
             Keep False for e.g. AREF definitions.
    """
    # Convert to database units.
    XY = [(round_to_db_unit(x), round_to_db_unit(y)) for x, y in XY]
    # Add start to the end if close required.
    if close and XY[0] != XY[-1]:
        XY.append(XY[0])
    # Remove consecutive duplicates.
    if unique:
        xy = [coor for coor, group in groupby(XY)]
    else:
        xy = XY
    n = len(xy)
    if n < min_length:
        return None
    rec = bytearray(pack_uint16(4 + 8 * n))  # Record length = 4 + 4*2*n
    rec.extend(pack_uint8(GDS_record.XY) + pack_uint8(GDS_datatype.INT32))
    for (x, y) in xy:
        rec.extend(pack_int32(x) + pack_int32(y))
    return rec


# Name of the referenced structure
def gds_sname(string):  # Record length = 4 + strlen
    name = bytes(string, encoding="UTF-8")
    l = len(name)
    if l & 1:
        l += 1
    return (
        pack_int16(4 + l)
        + pack_uint8(GDS_record.SNAME)
        + pack_uint8(GDS_datatype.ASCII)
        + pack_padstring(name)
    )


def gds_strans(bits):  # Record length = 4 + 2
    return (
        pack_int16(6)
        + pack_uint8(GDS_record.STRANS)
        + pack_uint8(GDS_datatype.BITARRAY)
        + pack_uint16(bits)
    )


def gds_mag(mag):  # Record length = 4 + 8
    return (
        pack_int16(12)
        + pack_uint8(GDS_record.MAG)
        + pack_uint8(GDS_datatype.REAL8)
        + pack_real8(mag)
    )


def gds_angle(ang):  # Record length = 4 + 8
    return (
        pack_int16(12)
        + pack_uint8(GDS_record.ANGLE)
        + pack_uint8(GDS_datatype.REAL8)
        + pack_real8(ang % 360)
    )


def gds_texttype(texttype):  # Record length = 4 + 2
    return (
        pack_int16(6)
        + pack_uint8(GDS_record.TEXTTYPE)
        + pack_uint8(GDS_datatype.INT16)
        + pack_int16(texttype)
    )


def gds_pathtype(pathtype):  # Record length = 4 + 2
    return (
        pack_int16(6)
        + pack_uint8(GDS_record.PATHTYPE)
        + pack_uint8(GDS_datatype.INT16)
        + pack_int16(pathtype)
    )


def gds_annotation(lay, xy, string, texttype=0):
    # here xy is a single [x,y] pair, not a list of pairs.
    # Note that datatype becomes texttype for annotations
    # Klayout throws an error on datatype in text.
    return (
        gds_text()
        + gds_layer(lay)
        + gds_texttype(texttype)
        + gds_xy([xy], False)
        + gds_string(string)
        + gds_endel()
    )


def gds_polyline(xy, w, lay, close=False, datatype=0, pathtype=0):
    xy_rec = gds_xy(xy, close, 2)
    if xy_rec is None:  # Less than 2 points in closed polyline.
        return bytearray(b"")
    return (
        gds_path()
        + gds_layer(lay)
        + gds_datatype(datatype)
        + gds_pathtype(pathtype)
        + gds_width(w)
        + xy_rec
        + gds_endel()
    )


def gds_polygon(xy, lay, datatype=0):
    xy_rec = gds_xy(xy, True, 4)
    if xy_rec is None:  # Less than 3 points in closed boundary.
        return bytearray(b"")
    return (
        gds_boundary() + gds_layer(lay) + gds_datatype(datatype) + xy_rec + gds_endel()
    )
