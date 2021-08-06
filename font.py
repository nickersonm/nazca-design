#!/usr/bin/env python3
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
# Polygon font output
#
# (c) 2017-2020 Xaveer Leijtens
# (c) 2017-2019 Ronald Broeke

"""
Module to add text to a layout.

Example:
    A simple way to put text in the layout using default settings::

        import nazca as nd

        nd.text(text='Hello world').put(0)
        nd.export_plt()
"""

import pickle
import warnings
import os

from . import netlist as net

from .util import md5
from .cfg import default_font as default_fnt  # prevent name clash with function


lss = 1.10  # line spacing scale
els = 1.10  # exclusion layer scale

empty_char = (  # Stub for unknown character
    35,  # Character width
    (
        (
            (0, 22.5),  # Thin rectangle with left-side opening
            (0, 42.5),
            (5, 42.5),
            (5, 27.5),
            (30, 27.5),
            (30, 72.5),
            (5, 72.5),
            (5, 57.5),
            (0, 57.5),
            (0, 77.5),
            (35, 77.5),
            (35, 22.5),
        ),
    ),
)

# Path to this file
modulepath = os.path.dirname(__file__)


class Font:
    """Class to define text object in various fonts.

    Nazca font files have the following structure:
    version: font file version number
    height: total height of the font
    space: interword space (=0 for fixed fonts)

    font: dictionary of characters:
        key: character
        w: width of that character
        p: list of polygons, of which each is a list of (x,y) coordinate points
    """

    def __init__(self, fontfile=default_fnt, box_layer=None, box_buf=None):
        """Initialize a font object.

        Args:
            fontfile (str): name of the fontfile (.nft)
            box_layer (str | int | tuple): layer reference to generate a square
                box in behind the text, e.g. for NOFILL (default = None)
            box_buf (float): extra buffer for the box_layer in um

        Returns:
            None
        """
        self.fontfile = fontfile
        self.box_layer = box_layer
        self.box_buf = box_buf

        if len(fontfile) > 4 and fontfile[-4:] == ".nft":
            ff = fontfile
        else:
            ff = fontfile + ".nft"
        if not os.path.exists(ff):
            ff = os.path.join(modulepath, "fonts", ff)
            if not os.path.exists(ff):
                print("Can't open font '{}'. Using '{}'.".format(fontfile, default_fnt))
                ff = os.path.join(modulepath, "fonts", default_fnt + ".nft")
        # short name w/o extension
        self.ff = os.path.splitext(os.path.basename(ff))[0]
        with open(ff, "rb") as f:
            self.version, self.lh, self.spc, self.tfont = pickle.load(f)

    def textheight(self, text, height=50):
        """Returns:
        float: total textheight of line(s) in text
        """
        return (text.count("\n") * lss + 1) * height

    def linelength(self, text, height=50):
        """Returns:
        float: total linelength of longest line in (multiline) text
        """
        scale = height / self.lh
        maxl = 0
        for line in text.split("\n"):
            # return width of line
            l = sum(self.tfont.get(char, empty_char)[0] for char in line) * scale
            l += (len(line) - 1) * self.spc * scale
            maxl = max(maxl, l)

        return l

    def text(
        self,
        text="Text",
        height=None,
        strokewidth=None,
        layer=1,
        align="lb",
        box_layer=None,
        box_buf=None,
        instantiate=False,
        cellname=None,
    ):
        """Output a cell with text of specified dimensions and alignment.

        In case of an exclusion (NOFILL) box_layer, e.g. for tiling protection,
        the box will be scaled with the text height unless a box_buf value
        is provided. The box_buf represents an absolute
        buffer between the text field and the box edge in um.

        Args:
            text (str): text to display
            height (float): height of the text in um (default 50)
            strokewidth (float): approximate font stroke width in um (height
                takes precedence if both are specified)
            layer (int | string): layer number or name to place text in
            align (str): relative placement: regex: [lcr][bct] (default = 'lb')
                examples: 'lc', 'cc', 'tr'
            box_layer (str | int | tuple): layer reference to generate a square
                box in behind the text, e.g. for NOFILL (default = None)
            box_buf (float): extra buffer for the box_layer in um
            instantiate (bool): instantiate cell (default=False)
            cellname (str): overrule the standard hash text cell naming (default=None)
                Set instantiate=True to get the text as a cell in the gds export.
            xs (str): not yet implemented. optional xsection to place text in. Useful to copy text to
               multiple layers, e.g. a logical layer.

        Returns:
            Cell: cell with text as provided in <text>

        Example::

            import nazca as nd

            message = nd.text(text='Hello world', align='cc')
            message.put(0)
            nd.export_plt()
        """
        if height is None:
            if strokewidth is None:
                height = 50
            else:
                # Scale font according to the width of the "|" character.
                polygon = self.tfont["|"][1][0]  # Just the first polygon
                _x = [x for x, y in polygon]
                height = self.lh * strokewidth / (max(_x) - min(_x))
                print(strokewidth, self.lh, height)
        if box_layer is None:
            box_layer = self.box_layer
        if box_buf is not None:
            box_buf = self.box_buf
        if box_buf is None:
            box_buf = height * (els - 1)

        if not isinstance(layer, (list, tuple)):
            layers = [layer]
        else:
            layers = layer

        halign, valign = align.lower()
        if halign not in "lcr":
            warnings.warn("text halign parameter not one of lcr, default to 'l'.")
            halign = "l"
        if valign not in "bct":
            warnings.warn("text valign parameter not one of bct, default to 'b'.")
            valign = "b"

        # Character parameters
        scale = height / self.lh
        ls = lss * height  # line spacing

        id5 = md5(text + align.lower() + "{}{}{}".format(self.ff, height, layer))
        if cellname is None:
            cellname = "text_{}_{}".format(self.ff, id5)
        with net.Cell(cellname, instantiate=instantiate) as txt:
            txt.auxiliary = True
            bby = self.textheight(text, height)
            if valign == "c":
                ypos = bby / 2 - height
            elif valign == "t":
                ypos = -height
            else:
                ypos = bby - height

            for line in text.split("\n"):
                bbx = self.linelength(line, height)
                if halign == "c":
                    xpos = -bbx / 2
                elif halign == "r":
                    xpos = -bbx
                else:
                    xpos = 0
                xpos0 = xpos

                for char in line:
                    # w, polys = self.tfont[char]
                    w, polys = self.tfont.get(char, empty_char)
                    for poly in polys:
                        XY = [(x * scale, y * scale) for x, y in poly]
                        if len(XY) > 3:
                            for lay in layers:
                                net.Polygon(points=XY, layer=lay).put(xpos, ypos, 0),
                    xpos = xpos + (w + self.spc) * scale
                if box_layer is not None:
                    net.Polygon(
                        layer=box_layer,
                        points=[
                            (0, 0),
                            (0, height + 2 * box_buf),
                            (bbx + 2 * box_buf, height + 2 * box_buf),
                            (bbx + 2 * box_buf, 0),
                        ],
                    ).put(xpos0 - box_buf, ypos - box_buf, 0)
                ypos -= ls

            net.Pin(name="b0", xs=None).put(xpos, ypos + ls, 0)
        return txt

    def sample(self, height=50, layer=1, instantiate=False):
        """Return a cell with all characters present in a font.

        Args:
            height (float): height off the text in um
            layer (int | string): layer number or name to place text in
            instantiate (bool): instantiate cell (default False)

        Returns:
            Cell: cell with a sample of all characters in the font.

        Example::

            cousine = nd.Font("cousine")
            sample = cousine.sample()
            sample.put()
            nd.export_plt()
        """
        id5 = md5("{}_{}_{}".format(self.ff, height, layer))
        with net.Cell(
            "sample_{}_{}".format(self.ff, id5), instantiate=instantiate
        ) as S:
            s = self.characters()
            n = 15
            lines = "\n".join([s[i : i + n] for i in range(0, len(s), n)])
            self.text(
                lines, height=height, align="cc", layer=layer, instantiate=False
            ).put()
        return S

    def characters(self):
        """Get all font characters.

        Returns:
            str: string with all characters present in a font
        """
        return "".join(sorted(list(self.tfont.keys())))


# Implement default font callable via nazca.text()
fdef = Font()


def default_font(f=default_fnt):
    global fdef
    fdef = Font(f)


def text(
    text="Text",
    height=None,
    strokewidth=None,
    layer=1,
    stretch=1,
    align="lb",
    linewidth=1,
    instantiate=False,
    box_layer=None,
    box_buf=None,
    cellname=None,
):
    return fdef.text(
        text=text,
        height=height,
        strokewidth=strokewidth,
        layer=layer,
        align=align,
        box_layer=box_layer,
        instantiate=instantiate,
        cellname=cellname,
    )


def linelength(text, height=50):
    return fdef.linelength(text, height)
