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
# (c) 2016-2019 Ronald Broeke
#

"""
Generate layout to gds, plt, svg or python script.

Example::

    import nazca as nd

    nd.export_gds()
    nd.export_plt()
    nd.export_svg()
"""

from collections import defaultdict
from math import sin, cos, radians
from shutil import copy2
import os
from IPython.core.display import display, HTML
from collections import defaultdict
from pprint import pprint

# for uPDK
import yaml
import json

import matplotlib.pyplot as mplt
from matplotlib.patches import Polygon as matPolygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import rgb2hex

from . import gds as gdsmod
from . import gds_base as gbase
from . import gds_import as gstream
from .netlist import Netlist, Cell, Pin, Polygon, Polyline, Annotation
from .mask_layers import add_layer, get_layer, get_xsection, set_layercolor

from . import cfg
from . import util
from .logging import logger
from .polysplitter import limit_polygon_points

try:
    import svgwrite
    cfg.SVGWRITE = True
except:
    cfg.SVGWRITE = False

try:
    from plotly.offline import plot  # creates new browser tab
    from plotly.offline import init_notebook_mode, iplot  # inline plots
    init_notebook_mode(connected=True)
    import plotly.graph_objects as go
    cfg.PLOTLY = True
except:
    cfg.PLOTLY = False


# =============================================================================
# Nazca cell
# =============================================================================
class ClsNazca:
    """Helper class to handle Nazca cell tree transformations.

    Note that this class does not reconstruct the netlist in the
    cell it returns, but manipulates it on element and layer level for the
    purpose of GDS export.

    Use cases:
    flatten a celltree, change and/or delete or filter layers, cell names,
    and masl elements.
    """

    def __init__(
        self,
        instantiate=True,
        cellmap=None,
        layermap=None,
        flat=None,
        hierarchy=None,
        infolevel=0,
        prefix = '',
        suffix = '_N',
    ):
        """Initialize a Clsazca object.

        Args:
            instantiate (bool):
            cellmap (dict): mapping of old_cellname: new_cellname
            layermap(dict): mapping of old_layer: new_layer
            flat (bool): flatten mask output
            infolevel (int):

        Returns:
            None
        """
        self.cellsopen = []
        self.flat = flat
        self.instantiate = instantiate
        self.layermap = layermap
        self.cellmap = cellmap
        self.oldcellrefs = {}  # link old cell to new cell
        self.newcellrefs = {}  # link new cell to old cell
        self.topcell = None
        self.infolevel = infolevel
        self.hierarchy = hierarchy
        self.prefix = prefix
        self.suffix = suffix

    def open(self, params, new_cell_name=None):
        """Open new Cell object based on namedtuple <cellinfo>.

        Args:
            cell (Cell): original cell to be copied/filtered
            level (int): cell level in hierarchy (0 is topcell)
            new_cell_name (str): overrule params.new_cell_name

        Returns:
            None
        """
        self.level = params.level

        if new_cell_name is None:
            cell_name = self.prefix + params.cell.cell_name + self.suffix
        else:
            cell_name = new_cell_name
        if self.infolevel > 0:
            print(
                "{:02d} {}ClsNazca: open  '{}' source:'  a1_black {}'".format(
                    self.level, "  " * self.level, cell_name, params.cell.cell_name
                )
            )
        if params.hierarchy == 'full':
            instantiate = True
        else:
            instantiate = params.cell.instantiate
        newcell = Cell(name=cell_name, instantiate=instantiate)
        self.cellsopen.append(newcell)
        self.oldcellrefs[newcell] = params.cell
        self.newcellrefs[params.cell] = newcell
        if params.level == 0:
            self.topcell = newcell
        # copy all pins from original into the new cell:
#        for name, node in params.cell.pin.items():
#            Pin(name, width=node.width, xs=node.xs, type=node.type,
#                chain=node.chain, show=node.show, remark=node.remark
#            ).put(node.xya())
        return newcell

    def close(self, params):
        """Close Cell object.

        Returns:
            ClsNazca:
        """
        newcell = self.cellsopen.pop()
        newcell.close()
        if self.infolevel > 0:
            print(
                "{:02d} {}ClsNazca: close '{}'".format(
                    self.level, "  " * self.level, newcell.cell_name
                )
            )
        return newcell

    def add_polygon(self, layer, xy):
        """Add a polygon.

        Args:
            layer (str): layername
            xy (list of (float, float)): list of (x, y)

        Returns:
            None
        """
        Polygon(points=xy, layer=layer).put(0)
        # print("ClsNazca: add_polygon to", cfg.cells[-1].cell_name, xy[:3], layer)

    def add_polyline(self, layer, xy, width=1.0):
        """Add a polyline.

        Args:
            layer (str): layername
            xy (list of (float, float)): list of (x, y)

        Returns:
            None
        """
        Polyline(points=xy, layer=layer, width=width).put(0)
        # print("ClsNazca: add_polylines")

    def add_annotation(self, text, layer, pos):
        """Add an annitation.

        Args:
            text (str): text to use as annotation
            layer (str): layer name
            pos ((float, float)): annotation posistion (x, y)

        Returns:
            None
        """
        # print("ClsNazca: add_annotation")
        Annotation(text=text, layer=layer).put(*pos)

    def add_instance(self, inode, xya, flip):
        """Add instance.

        Args:
            inode (instance Node): Node as Instance reference
            xya ((float, float, float)): coordinate to put instance
            flip (bool): flip or mirror state of the instance placement

        Returns:
            None
        """
        if self.infolevel > 0:
            print(
                "{:02d} {}ClsNazca: add_instance '{}' xya:{}, flip:{}".format(
                    self.level, "  " * self.level, inode.cell.cell_name, xya, flip
                )
            )
        cell = inode.cell
        #if cell.instantiate:  # and not flat:
            # TODO, flatten array's here if not instantiated?

        newcell = self.newcellrefs.get(cell, None)
        if newcell is None:
            raise Exception(f"Can not put cell '{cell.cell_name}' as it has not been opened yet:"
                " Do not use hierarchy='apply' in output to a Nazca tree."
            # TODO: Can be solved by adding instances at cell closure rather than opening?
            )
        else:
            newcell.put(
                'org',
                *xya,
                'org',
                flip=inode.flip,
                array=inode.array,
                scale=inode.scale
            )

    def add_gds(self):
        """Not implemented.

        Returns:
            None
        """
        if self.infolevel > 0:
            print("ClsNazca: add_gds")



# =============================================================================
# GDS
# =============================================================================
cfg.force_nazca_version = None  # override auto nazca versioning if not None.


class ClsGDS:
    """Helper class to handle gdsii compatible export of masks."""

    def __init__(self):
        """Initialize a ClsGDS object.

        Returns:
            None
        """
        pass

    def open(self, filebasename):
        """Open ClsGDS file output.

        Args:
            filebasename (str): base filename

        Returns:
            None
        """
        self.filename = filebasename + ".gds"
        self.tempname = os.path.join(
            os.path.dirname(filebasename), ".tmp_" + os.path.basename(filebasename) + ".gds"
        )
        self.outfile = open(self.tempname, "wb")
        self.outfile.write(gdsmod.layout_open(name=cfg.force_nazca_version))
        self.content = defaultdict(list)  # collect cell content per (instaniated) level.

    def write(self, level):
        """Write to ClsGDS output file.

        Args:
            level (int): level in hierarchy to write to

        Returns:
            None
        """
        self.outfile.write(b"".join(self.content[level]))

GDS = ClsGDS()

# ==============================================================================
# svg
# ==============================================================================
class ClsSVG:
    """Helper class to handle svg compatible export of masks."""

    def __init__(self):
        """Initilaize a ClsSVG object.

        Returns:
            None
        """
        # if filename is not None:
        #    filename = path.splitext(filename)[0] + '.svg'
        self.minx = 1e6
        self.miny = 1e6
        self.maxx = -1e6
        self.maxy = -1e6

    def open(self, filebasename="nazca_export"):
        """Open SVG.

        Open after layer colors have been loaded.

        Args:
            basename (str): base filename.

        Returns:
            None
        """
        self.filename = filebasename + ".svg"
        # create a map of defined layer colors to speed up lookup
        self.drawing = svgwrite.Drawing(
            filename=self.filename, size=("100%", "100%"), profile="tiny", debug=False
        )
        self.layer2color = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ["color_name"]].iterrows():
                self.layer2color[row[0]].append(dict(cfg.colors.iloc[i]))
        self.cmap = cfg.plt_cmap
        self.cmaplen = len(self.cmap)
        self.alpha = cfg.plt_alpha

    def add_polygon(self, layer, points, bbox):
        """Add a polygon to the SVG drawing.

        Args:
            layer (str): layername
            points (list of (float, float)): polygon points [(x, y), ...]
            bbox ((list of (float, float)): bounding box the points [(x, y), ...]

        Returns:
            None
        """
        L, D, T = cfg.layername2LDT[get_layer(layer)]
        color_name = "{}/{}".format(L, D)
        for colors in self.layer2color[color_name]:
            if not colors["visible"]:
                continue
            edgecolor = colors["frame_color"]
            facecolor = colors["fill_color"]
            alpha = colors["alpha"]
            lw = colors["width"]
            if colors["dither_pattern"] == "I1":
                fill_opacity = 0
            else:
                fill_opacity = alpha
            self.g = self.drawing.g(
                stroke=edgecolor,
                stroke_opacity=alpha,
                fill=facecolor,
                fill_opacity=fill_opacity,
                stroke_width=lw,
            )
            self.mask = self.drawing.add(self.g)

            p = self.drawing.polygon(points=points, transform="scale(1,-1)")
            self.mask.add(p)

            if self.minx > bbox[0]:
                self.minx = bbox[0]
            if self.maxx < bbox[2]:
                self.maxx = bbox[2]
            if self.miny > bbox[1]:
                self.miny = bbox[1]
            if self.maxy < bbox[3]:
                self.maxy = bbox[3]

    def close(self):
        """Close the svg mask layout output file.

        Returns:
            None
        """
        x, y, w, h = self.minx, self.miny, self.maxx - self.minx, self.maxy - self.miny
        self.drawing.viewbox(minx=x, miny=-y - h, width=w, height=h)
        self.drawing.save()

SVG = ClsSVG()


# ==============================================================================
# matplotlib export
# ==============================================================================
class ClsMatplotlib:
    """Helper class to handle matplotlib export of masks."""

    def __init__(self):
        """Construct an object handling Matplotlib exports.

        Returns:
            None
        """

    #  try:
    #      font = cfg.matplotlib_font
    #
    #  except:
    #  font = {
    #      # 'family'  : 'normal',
    #      'style'  : 'normal',
    #      'weight' : 'light', # 'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
    #      'size'   : cfg.plt_fontsize}
    #  mplt.rc('font', **font)

    def open(self, cell=None, title=None, show_pins=True):
        """Inititialize Matplotlib mask output.

        Args:
            cell (Cell): Cell to draw
            title (str): title in Matplotlib output

        Returns:
            None
        """
        font = {
            #'family'  : 'normal',
            "style": "normal",
            "weight": "light",  #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
            "size": cfg.plt_fontsize,
        }
        mplt.rc("font", **font)

        self.cell = cell

        if title is not None:
            # self.title = '{} - {}'.format(title, cell.cell_name)
            self.title = "{}".format(title)
        else:
            self.title = cell.cell_name
            self.title = "cell: {}".format(self.title)

        self.cellname = cell.cell_name
        self.patches = []
        self.colors = []
        self.minx = 1e6
        self.miny = 1e6
        self.maxx = -1e6
        self.maxy = -1e6
        self.figsize = cfg.plt_figsize
        self.show_pins = show_pins

        # Create a layer-to-color map of type list of
        #   defined layer colors to speed up look-ups.
        self.layer2color = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ["color_name"]].iterrows():
                self.layer2color[row[0]].append(dict(cfg.colors.iloc[i]))

        # backup colormap for layers without color setting
        self.cmap = cfg.plt_cmap
        self.cmaplen = len(self.cmap)
        self.alpha = cfg.plt_alpha

    def add_polygon(self, layer, points, bbox):
        """Add a polygon to the plot.

        Note that using linewidth > 0 significantly slows down drawing in Matplotlib.

        Args:
            layer (int, int): mask layer
            points (list of (float, float)): polygon points
            bbox (tuple of floats): bounding box of polygon (x1, y1, x2, y2)

        Returns:
            None
        """
        if layer == cfg.default_layers["docu_pin"]:
            return None
        polygon = None

        layerID = get_layer(layer)
        L, D, T = cfg.layername2LDT[layerID]
        color_name = "{}/{}".format(L, D)
        if color_name not in self.layer2color.keys():  # add back-up color
            col = rgb2hex(self.cmap[L % self.cmaplen])
            addcolor = set_layercolor(
                layer=layerID, frame_color=col, fill_color=col, width=0
            )
            self.layer2color[color_name].append(dict(addcolor.iloc[0]))
            logger.info(
                "Added color to {}, {}".format(layerID, cfg.layername2LDT[layerID])
            )

        # Add polygon for each color:
        #   like in Klayout there can be multiple colors per layer.
        for colors in self.layer2color[color_name]:
            if not colors["visible"]:
                continue
            edgecolor = colors["frame_color"]
            facecolor = colors["fill_color"]
            alpha = colors["alpha"]
            lw = float(colors["width"])
            if colors["dither_pattern"] == "I1":
                fill = False
            else:
                fill = True
            # lw = 0 # matplotlib very slow with lw != 0
            polygon = matPolygon(
                points,
                closed=True,
                edgecolor=edgecolor,
                facecolor=facecolor,
                fill=fill,
                lw=lw,
                alpha=alpha,
            )

        if polygon is not None:
            self.patches.append(polygon)
            if self.minx > bbox[0]:
                self.minx = bbox[0]
            if self.maxx < bbox[2]:
                self.maxx = bbox[2]
            if self.miny > bbox[1]:
                self.miny = bbox[1]
            if self.maxy < bbox[3]:
                self.maxy = bbox[3]

            addbbox = False
            if addbbox:
                polygon = matPolygon(
                    [
                        (bbox[0], bbox[1]),
                        (bbox[2], bbox[1]),
                        (bbox[2], bbox[3]),
                        (bbox[0], bbox[3]),
                    ],
                    closed=True,
                    edgecolor="#00FF00",
                    facecolor=facecolor,
                    fill=False,
                    lw=1,
                    alpha=alpha,
                )
                self.patches.append(polygon)

    def add_polyline(self, layer, points, width, bbox):
        """Add a polyline to the plot.

        Not implemented.
        Polylines are transformed into polygons and send there.
        """

    def close(self):
        """Close plt and show plot.

        Returns:
            None
        """
        if not self.patches:
            print("Nothing to draw.")
            return None

        dx = self.maxx - self.minx
        dy = self.maxy - self.miny
        bufx = 0.05 * dx
        bufy = 0.05 * dy
        buf = max(bufx, bufy)
        domain = max(dx, dy)
        xtot = dx + 2 * buf
        ytot = dy + 2 * buf
        aspect = ytot / xtot
        factor = 1.0
        if aspect > factor:
            ysize = self.figsize
            xsize = ysize / aspect
        else:
            xsize = self.figsize
            ysize = 1.2 * xsize * aspect

        if cfg.plt_fig is None:
            fig, ax = mplt.subplots(
                figsize=(xsize, ysize), facecolor=cfg.plt_background_outside
            )
        else:
            fig, ax = cfg.plt_fig
        p = PatchCollection(self.patches, match_original=True)

        ax.add_collection(p)
        ax.set(
            xlim=[self.minx - buf, self.minx + dx + buf],
            ylim=[self.miny - buf, self.miny + dy + buf],
        )

        # format spines:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.spines["left"].set_visible(True)
        #ax.spines["left"].set_smart_bounds(True)

        ax.spines["bottom"].set_visible(True)
        #ax.spines["bottom"].set_smart_bounds(True)

        ax.yaxis.set_ticks_position("left")
        # ax.yaxis.set_ticklabels('')
        ax.xaxis.set_ticks_position("bottom")

        ax.set_title(self.title)

        # cell info available under self.cell:
        if self.show_pins:
            if not cfg.documentation_pins:
                # cfg.pinstyles['default']['layer'] != cfg.default_layers['docu_pin']:
                pin_ignore = cfg.pinstyle.get("pin_ignore", [])
                dt = 0.02 * domain
                for name, pin in self.cell.pin.items():
                    if (
                        name == "org"
                        or (pin.type == "bbox" and not cfg.bbox_stubs)
                        or name in pin_ignore
                    ):
                        continue
                    at = pin.pointer.a
                    xt = pin.pointer.x + dt * cos(radians(at))
                    yt = pin.pointer.y + dt * sin(radians(at))
                    ax.text(xt, yt, name, ha="center", va="center", rotation=0, wrap=True)
            else:  # docu-pins
                docu_patches = []
                d_text = 0.40
                for name, pin in self.cell.pin.items():
                    # customize appearance of each pin
                    if pin.xs is None:
                        style = "docu_default"
                    else:
                        style = get_xsection(pin.xs).pinstyle
                    if style is None:
                        style = "docu_default"
                    if pin.io == 1:
                        flip = True
                    else:
                        flip = False
                    pin_scale = cfg.pinstyles[style].get("scale", cfg.pinstyles['docu_default']['scale'])
                    pin_shape = cfg.pinstyles[style].get("shape", cfg.pinstyles['docu_default']['shape'])
                    pin_ignore = cfg.pinstyles[style].get("pin_ignore", cfg.pinstyles['docu_default']['pin_ignore'])
                    fontsize = 20 * 15 * pin_scale
                    scale = domain * pin_scale
                    shape = cfg.pinshapes[pin_shape]
                    dt = -d_text * pin_scale * domain

                    # not pin.show or
                    if (
                        name == "org"
                        or (pin.type == "bbox" and not cfg.bbox_stubs)
                        or name in pin_ignore
                    ):
                        continue
                    at = pin.pointer.a
                    xt = pin.pointer.x + dt * cos(radians(at))
                    yt = pin.pointer.y + dt * sin(radians(at))

                    ax.text(
                        xt,
                        yt,
                        name,
                        fontsize=fontsize,
                        rotation=0.0,
                        ha="center",
                        va="center",
                        weight="bold"
                        # bbox=dict(boxstyle="round", ec=(1, 1, 1), fc=(1, 1, 1))
                    )
                    points = util.transform_polygon(
                        points=shape,
                        da=at,
                        scale=scale,
                        flipy=flip,
                        x=pin.pointer.x,
                        y=pin.pointer.y
                    )
                    polygon = matPolygon(
                        xy=points,
                        closed=True,
                        edgecolor= cfg.pinstyles[style].get('edgecolor', "#000000"),
                        facecolor= cfg.pinstyles[style].get('fillcolor', "#FFFFFF"),
                        fill=True,
                        lw=3,
                        alpha=0.5,
                    )
                    docu_patches.append(polygon)
                if docu_patches:
                    p = PatchCollection(docu_patches, match_original=True)
                    ax.add_collection(p)

        mplt.axis("scaled")
        mplt.tight_layout()
        # ax.set_facecolor(cfg.plt_background_inside)

# ==============================================================================
#         #Watermark:
#         ax.set_facecolor((0.97, 0.97, 0.97))
#         import matplotlib.image as image
#         datafile = '/nazca/nazca-logo_illustration_beeldmerk.png'
#         im = image.imread(datafile)
#         mplt.imshow(im, extent=[self.minx-bufx, self.minx+dx+bufx, self.miny-bufy, self.miny+dy+bufy],
#             alpha=0.05, zorder=3)
#
# ==============================================================================
PLT = ClsMatplotlib()


# ==============================================================================
# plotly export
# ==============================================================================
class ClsPlotlylib:
    """Helper class to handle Plotly export of masks.

    Need to be run in a Jupyter notebook due to JavaScript output interface.
    """

    def __init__(self):
        """Construct an object handling Plotly exports.

        Returns:
            None
        """
        pass


    def open(self, cell=None, title=None, show_pins=True):
        """Inititialize Plotly mask output.

        Args:
            cell (Cell): Cell to draw
            title (str): title in Matplotlib output

        Returns:
            None
        """
        self.cell = cell

        if title is not None:
            # self.title = '{} - {}'.format(title, cell.cell_name)
            self.title = "{}".format(title)
        else:
            self.title = cell.cell_name
            self.title = "cell: {}".format(self.title)

        self.cellname = cell.cell_name
        self.show_pins = show_pins
        self.fig = go.Figure()

        # Create a layer-to-color map of type list of
        #   defined layer colors to speed up look-ups.
        self.layer2color = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ["color_name"]].iterrows():
                self.layer2color[row[0]].append(dict(cfg.colors.iloc[i]))

        # backup colormap for layers without color setting
        self.cmap = cfg.plt_cmap
        self.cmaplen = len(self.cmap)
        self.alpha = cfg.plt_alpha


    def add_polygon(self, layer, points, bbox):
        """Add a polygon to the plot.

        Args:
            layer (int, int): mask layer
            points (list of (float, float)): polygon points
            bbox (tuple of floats): bounding box of polygon (x1, y1, x2, y2)

        Returns:
            None
        """
        if layer == cfg.default_layers["docu_pin"]:
            return None

        layerID = get_layer(layer)
        L, D, T = cfg.layername2LDT[layerID]
        color_name = "{}/{}".format(L, D)
        if color_name not in self.layer2color.keys():  # add back-up color
            col = rgb2hex(self.cmap[L % self.cmaplen])
            addcolor = set_layercolor(
                layer=layerID, frame_color=col, fill_color=col, width=0
            )
            self.layer2color[color_name].append(dict(addcolor.iloc[0]))
            logger.info(
                "Added color to {}, {}".format(layerID, cfg.layername2LDT[layerID])
            )

        x, y  = [list(i) for i in zip(*points)]
        try:
            width = int(self.layer2color[color_name][0]['width'])
            # not type casting give a nan error on width in plotly
        except:
            width = 0
        self.fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                fill="toself",
                marker=dict(opacity=self.alpha),
                opacity=1 - self.alpha,
                fillcolor=self.layer2color[color_name][0]['fill_color'],
                line=dict(
                    color=self.layer2color[color_name][0]['frame_color'],
                    width=width,
                )
            )
        )


    def add_polyline(self, layer, points, width, bbox):
        """Add a polyline to the plot.

        Not implemented.
        Polylines are transformed into polygons and send there.
        """

    def close(self):
        """Close ploy and show plot.

        Returns:
            None
        """
        self.fig.update_layout(
            yaxis=dict(
                scaleanchor="x",
                scaleratio=1,
            ),
            showlegend=False,
            width=1000,
            height=1000,
            autosize=True
        )
        iplot(self.fig)
# end class ClsPlotlylib
PLOTLY = ClsPlotlylib()


# =============================================================================
# Nazca building block
# =============================================================================
class BBlock():
    """Helper class to export a layout as a Nazca GDS building block.

    The class create a .py module with a Cell loading the GDS file and placing the pins.
    Set bb=True in export_gds() to use it.
    """

    def __init__(self, filebasename, outpath="", bbpath=""):
        """Initialize BB file: open and add header.

        Args:
            filebasename (str): basename of the output file (without extension and no path).
            outpath (str): output directory
            bbpath (str): path to gds files on the gds py module that loads the BB.

        Returns:
            None
        """
        self.filename = os.path.join(outpath, filebasename)
        self.filebasename = filebasename
        self.outpath = outpath
        if bbpath == "":
            self.bbpath = 'localdir'
        else:
            self.bbpath = bbpath
        self.modulename = os.path.basename(filebasename).replace(".", "_") + "_gds"
        #self.path = os.path.dirname(filebasename)
        self.pyfilename = self.modulename + ".py"
        self.bblockfile = open(os.path.join(self.outpath, self.pyfilename), "w")
        bbfile_header = """# Nazca building block(s) file.
#
# Note1: Do *not* edit coordinates in the cell definition(s) below.
# Note2: Always use this file with the intended, accompanying gds library.
#
# In order to use the BB(s) in this file,
# include the following lines in your design file:
#
# import {0}
# {0}.<bb_name>.put()
#
# where '{0}' is the name of this file without .py extension and path, and
# where <bb_name> needs to be replaced with any cell variable name in this file.
#
# Executing this file stand-alone exports a gds with all BBs in this file.

from os.path import dirname
from nazca import Cell, Pin, load_gds, export_gds, put_stub
""".format(self.modulename, self.pyfilename)

        self.bblockfile.write(bbfile_header)
        self.bblockfile.write(f"localdir = dirname(__file__)")
        self.xpos = 0
        self.foot = ""
        self.BBcells = []

    def body(self, topcell):
        """Write BB file body of element.

        Args:
            topcell (Cell): topcell to export as BB

        Returns:
            None
        """
        indent = "    "
        # gdsfilename = os.path.join(self.bbpath, f"{self.filebasename}.gds")
        gdsfilename = f"localdir+'/{self.filebasename}.gds'"
        cellname = topcell.cell_paramsname
        cellvar = cellname.replace(".", "_")
        self.bblockfile.write(
            f"\nwith Cell(name='{cellname}', autobbox=False, instantiate=True) as {cellvar}:\n"
        )
        self.bblockfile.write(f"{indent}{cellvar}.default_pins('{topcell.default_in}', '{topcell.default_out}')\n")
        self.bblockfile.write(
            f"{indent}load_gds(filename={gdsfilename}, cellname='{cellname}', instantiate=False).put(0)\n"
        )
        self.foot += "{}{}.put('org', {:.1f})\n".format(indent, cellvar, self.xpos)
        self.BBcells.append(cellvar)

        try:
            length = topcell.length
        except:
            length = 0
        if length is not None:
            self.xpos += 1.1 * length
        for name, pin in sorted(topcell.pin.items()):
            if name == "org":
                continue
            x, y, a = pin.pointer.xya()
            if pin.xs is None:
                xs = None
            else:
                xs = "'{}'".format(pin.xs)
            self.bblockfile.write(
                f"{indent}Pin(name='{name}', width={pin.width}, xs={xs}).put({x:.5f}, {y:.5f}, {a:.5f})\n"
            )
        self.bblockfile.write(f"{indent}put_stub()")
        return None

    def footer(self):
        """Write BB file footer.

        Returns:
            None
        """
        self.bblockfile.write("\nBBcells = [{}]\n".format(', '.join(self.BBcells)))
        self.bblockfile.write("\nif __name__ == '__main__':\n")
        self.bblockfile.write(self.foot)
        self.bblockfile.write("    export_gds()\n")
        print("...wrote module file '{}'".format(self.pyfilename))
        return None

# =============================================================================
# uPDK building block
# =============================================================================
class UPDK():
    """Helper class to export a layout as a Nazca uPDK building block.
    """

    def __init__(self, filebasename):
        """Initialize uPDK file: open and add header.

        Args:
            filebasename (str): basename of the output file (without extension).

        Returns:
            None
        """
        self.filename = filebasename
        self.modulename = os.path.basename(filebasename).replace(".", "_") + "_gds"
        self.path = os.path.dirname(filebasename)
        self.uPDKfilename = self.modulename + ".yaml"
        self.uPDKfile = open(os.path.join(self.path, self.uPDKfilename), "w")
        self.uPDKfile_header = """#uPDK exported from Nazca Design"""
        self.blocks = {}


    def roundlist(self, L, digits=6):
        """Avoid very small fractions in floats."""
        return [round(l, digits) for l in L]


    def body(self, topcell):
        """Write uPDK file body of element.

        Args:
            topcell (Cell): topcell to export as uPDK

        Returns:
            None
        """
        blockdict = {}
        indent = "    "
        gdsfilename = self.filename + ".gds"
        cellname = topcell.cell_paramsname
        cellvar = cellname.replace(".", "_")

        blockdict["cellname"] = cellname

        blockdict["bbox"] = self.roundlist(list(topcell.bbox))
        blockdict["userbbox"] = topcell.userbbox

        # add pins:
        pindict = {}
        for name, pin in sorted(topcell.pin.items()):
            if name in ["org", "lb", "lc", "lt", "tl", "tc", "tr", "rt", "rc", "rb", "br", "bc", "bl", "cc"]:
                continue
            x, y, a = pin.pointer.xya()
            if pin.xs is None:
                xs = None
            else:
                xs = pin.xs
            pindict[name] = {
                'doc': pin.remark,
                'width': round(pin.width, 6),
                'xs': xs,
                'xya': f"({x}, {y}, {a})",
            }
        blockdict['pins'] = pindict

        # add parameters:
        parametersdict = {}
        for name, prange in sorted(topcell.ranges.items()):
            parametersdict[name] = prange
        blockdict['parameters'] = parametersdict

        self.blocks[cellname] = blockdict
        return None

    def footer(self):
        """Write uPDK file footer.

        Returns:
            None
        """
        #pprint(self.blocks)
        blocks = {'blocks': self.blocks}
        with open(self.uPDKfilename, 'w') as F:
            #print(yaml.dump(json.loads(json.dumps(blocks))))
            yaml.dump(json.loads(json.dumps(blocks)), F)

        print("...wrote uPDK file '{}'".format(self.uPDKfilename))
        return None


tab = "  "
class Export_layout:
    """Class to export a mask layout to various formats.

    Export types:
        gds
        plt
        svg
        nazca

    The Class makes use of one or more of the following class objects in
    the layout.py module to handle format specific output::

        GDS
        PLT
        SVG
        self.CELL

    """

    def __init__(
        self,
        layermap=None,
        layermapmode="all",
        prefix = '',
        suffix = '_N',
    ):
        """Initialize an Export_layout object.

        Args:
            layermap (dict): layermap oldlayr:newlayer
            layermapmode (str): 'all': copy all layers,
                'none': copy only mapped layers, remove others.

        Returns:
            None
        """
        self.gds_files = dict()  # existing gds files used in the layout
        self.layermapmodes = ["none", "all"]
        self.prefix = prefix
        self.suffix = suffix
        self.reset()

        # make sure hull layer is known
        if cfg.show_hull:
            if "hull" in cfg.default_layers.keys():
                get_layer(cfg.default_layers["hull"])

        self.setlayermap(layermap=layermap, mode=layermapmode)
        return None

    def reset(self):
        """Reset export settings to defaults.

        Returns:
            None
        """
        self.topcells = None
        self.path = ""
        self._filename = "nazca_export"
        self.bblock = False
        self.bbpath = ""
        self.infolevel = 0
        self.show_cells = False
        self.clear = cfg.export_clear
        self.title = None
        self.output = None
        self._ascii = False  # gds as ascii
        self._gds = False
        self._flat = False
        self._plt = False
        self._svg = False
        self._plotly = False
        self.nazca = False
        self.uPDK=False

        self.CELL = None # only utilized when self.nazca = True
        # CELL is an attribute so it would e possible to create multiple
        # new Nazca trees inside a single .

        self.layermap = {}
        self.layermapmode = "all"
        self.setlayermap()
        self.md5 = False
        self.submit = False
        return None

#    def newtree(self):
#        """Create a Nazca export object.
#
#        Returns:
#            ClsNazca:
#        """
#        self._newtree = ClsNazca()
#        return self._newtree

    def setlayermap(self, layermap=None, mode="all"):
        """Create the layermap for export.

        Layermap is set internally.

        Args:
            layermap (dict): layermap oldlayr:newlayer
            mode (str): 'all' (default): copy all layers,
                'none': copy only mapped layers, remove others.
                Values are case insensitive.

        Returns:
            dict: {layer: export_layer}, layermap as reference for the user
        """
        self.layermap = {}
        if layermap is None:
            layermap = {}

        if mode is None:
            mode = "all"
        elif mode.lower() in self.layermapmodes:
            self.layermapmode = mode.lower()
        else:
            logger.warning(
                "Trying to set a non-existing layermapmode '{}'."
                " Valid options as {}".format(mode, self.layermapmodes)
            )
            self.layermapmode = "all"

        if self.layermapmode == "all":
            for layername in cfg.layername2LDT.keys():
                self.layermap[layername] = layername
        elif self.layermapmode == "none":
            for layername in cfg.layername2LDT.keys():
                self.layermap[layername] = None

        if layermap:
            existing_layers = list(cfg.layername2LDT.keys()).copy()
            for layerold, layernew in layermap.items():
                layers_old = get_layer(layerold, aslist=True)
                for layer_old in layers_old:
                    if layer_old in existing_layers:
                        if layernew is None:
                            self.layermap[layer_old] = None
                        else:
                            layer_new = get_layer(layernew, aslist=True)
                            if len(layer_new) > 1:
                                print(f"WARNING: non unique mapping to layer {layernew}. Will take first entry and map to '{layer_new[0]}'")
                            elif layer_new[0] not in cfg.layername2LDT.keys():
                                logger.warning(
                                    "(rebuild): mapping to non-existing layer {0}."
                                    " Adding it for you, but honestly, you should add it:"
                                    " add_layer(layer={0})".format(layernew)
                                )
                                add_layer(name=str(layernew[0]), layer=layernew[0])
                            self.layermap[layer_old] = layer_new[0] #get_layer(layernew)
                    else:
                        logger.warning(f"Original layer {layerold} in layermap is not defined/unknown.")
        return self.layermap

    @property
    def filebasename(self):
        """Return filebasename.

        Returns:
            str: filebasename
        """
        return self._filebasename

    @property
    def filename(self):
        """Get filename for export.

        Returns:
            str: filename
        """
        return self._filename

    @filename.setter
    def filename(self, name):
        """Set filename and path for export.

        For names ending with .gds, .plt and .vsg extensions the corresponding
        export format is automatically selected.

        If a specific extension is missing in the provided filename when
        exporting to a specified format, like export_gds(), then the extension
        is appended to the filename automatically in during export.

        Args:
            name (str): filename. The name may include a path, in which case
                case the path becomes the base directory. Missing directories
                will be created.

        Returns:
            str: filename
        """
        if name is None:
            name = self._filename
        if isinstance(name, str):
            basename = os.path.basename(name)
            dirname = os.path.dirname(name)
            if dirname == "":
                if self.submit:
                    dirname = "submit/"
                else:
                    dirname = "./"
            if not os.path.exists(dirname):
                os.makedirs(dirname)
                logger.info("directory '{}' not existing. Creating it.".format(dirname))
            self.dirname = dirname
            last = basename[-4:].lower()
            base = basename[:-4]
            if last == ".gds":
                self.gds = True
                self._filebasename = base
            elif last == ".svg":
                self.svg = True
                self._filebasename = base
            elif last == ".plt":
                self.plt = True
                self._filebasename = base
            else:
                self._filebasename = basename
        else:
            logger.warning(
                "Filename provided is not a string but a {}."
                " Using '{}' instead.".format(type(name), self._filebasename)
            )
        return self._filebasename

    @property
    def svg(self):
        """Get svg export setting."""
        return self._svg
    @svg.setter
    def svg(self, value):
        """Set svg export setting."""
        if value:
            self.flat = True
        self._svg = value
        return self._svg

    @property
    def plt(self):
        """Get plt export setting."""
        return self._plt
    @plt.setter
    def plt(self, value):
        """Set plt export setting."""
        if value:
            self.flat = True
        self._plt = value
        return self._plt

    @property
    def plotly(self):
        """Get plotly export setting."""
        return self._plotly
    @plotly.setter
    def plotly(self, value):
        """Set plotly export setting."""
        if value:
            self.flat = True
        self._plotly = value
        return self._plotly

    @property
    def flat(self):
        """Get flat export setting.

        Returns:
            bool: flat state
        """
        return self._flat
    @flat.setter
    def flat(self, value):
        """Set flat export setting.

        Returns:
            bool: flat state"""
        if not value:
            self._svg = False
            self._plt = False
        self._flat = value
        return self._flat

    @property
    def ascii(self):
        """Get ascii export setting."""
        return self._ascii
    @ascii.setter
    def ascii(self, value):
        """Set flat export setting."""
        if value:
            self._gds = True
        self._ascii = value
        return self._ascii

    @property
    def gds(self):
        """Get gds export setting.

        Returns:
            bool: plt export state
        """

        return self._gds
    @gds.setter
    def gds(self, value):
        """Set gds export setting.

        Args:
            value (bool): gds export state

        Returns:
            bool: plt export state
        """
        if not value:
            self._ascii = False
        self._gds = value
        return self._gds


    # Functions to directly call self.CELL methods in Export_layout
    # for syntax reasons as user level.
    def open(self, params, new_cell_name=None):
        """Alias for self.CELL.open if self.nazca = True

        Returns:
            None
        """
        if self.nazca:
            self.CELL.open(params=params, new_cell_name=new_cell_name)
    cell_open = open

    def close(self, params):
        """Alias for self.CELL.open if self.nazca = True

        Returns:
            ClsNazca:
        """
        if self.nazca:
            NazcaCell = self.CELL.close(params)
            self.topcell = self.CELL.topcell
            return NazcaCell
    cell_close = close

    def put(self, *args, **kwargs):
        if self.nazca:
            self.CELL.topcell.put(*args, **kwargs)
    def topcell(self):
        if self.nazca:
            return(self.CELL.topcell)
    def add_polygon(self, layer, points):
        if self.nazca:
            self.CELL.add_polygon(xy=points, layer=layer)
    def add_polyline(self, layer, points, width):
        if self.nazca:
            self.CELL.add_polyline(layer=layer, xy=points, width=width)
    def add_annotation(self, *args, **kwargs):
        if self.nazca:
            self.CELL.add_annotation(*args, **kwargs)
    def add_instance(self, inode, trans, flip):
        if self.nazca:
            self.CELL.add_instance(inode=inode, xya=trans, flip=flip)
    def add_gds(self, *args, **kwargs):
        if self.nazca:
            self.CELL.add_gds(*args, **kwargs)


    def add_polygons(self, params):
        """Add polygons content to cell.

        Args:
            pgon_iter (iterator): polygon iterator
            level (int): hierarchy level
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        pgon_iter = params.iters["polygon"]
        level = params.parent_level

        for pgon, xy, bbox in pgon_iter:
            layermap = self.layermap[pgon.layer]
            if layermap is None:
                continue
            L, D = cfg.layername2LDT[layermap][0:2]
            if self.gds:
                pols = limit_polygon_points(xy)
                for p in pols:
                    GDS.content[level].append(gbase.gds_polygon(p, lay=L, datatype=D))
            if self.plt:
                PLT.add_polygon(layermap, xy, bbox)
            if self.svg:
                SVG.add_polygon(layermap, xy, bbox)
            if self.plotly:
                PLOTLY.add_polygon(layermap, xy, bbox)
            if self.nazca:
                self.CELL.add_polygon(layermap, xy)
        return None


    def add_polylines(self, params):
        """Add polylines content to cell.

        In Matplotlib and svg output the polylines are converted to polygons.

        Args:
            pline_iter (iterator): polyline iterator
            level (int): hierarchy level
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        pline_iter = params.iters["polyline"]
        level = params.parent_level

        for pline, xy, bbox in pline_iter:
            layermap = self.layermap[pline.layer]
            if layermap is None:
                continue
            L, D = cfg.layername2LDT[layermap][0:2]
            if self.gds:
                if not cfg.export_polyline_as_polygon:
                    GDS.content[level].append(
                        gbase.gds_polyline(
                            xy,
                            pline.width,
                            lay=L,
                            datatype=D,
                            pathtype=pline.pathtype,
                        )
                    )
                else:
                    pols = limit_polygon_points(xy)
                    for p in pols:
                        GDS.content[level].append(
                            gbase.gds_polygon(p, lay=L, datatype=D)
                        )
            if self._plt or self._svg:
                # No need to split: use polyline2polygon (not polyline2polygons).
                polygon_xy = util.polyline2polygon(
                    xy, width=pline.width, miter=pline.miter
                )
                # TODO: polygon_xy can also be supplied by the iterator
                if self.plt:
                    PLT.add_polygon(layermap, polygon_xy, bbox)
                if self.svg:
                    SVG.add_polygon(layermap, polygon_xy, bbox)
                if self.plotly:
                    PLOTLY.add_polygon(layermap, polygon_xy, bbox)
            if self.nazca:
                self.CELL.add_polyline(layer=layermap, xy=xy, width=pline.width)
        return None

    def add_annotations(self, params):
        """Add annotation content to cell.

        Args:
            instantiate (bool): the instantiation level of the cell where
               the annotation are from: To check for black boxes, which
               can't be flattened.
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        anno_iter = params.iters["annotation"]
        level = params.parent_level
        create = params.cell_create
        cell = params.cell

        for anno, xy in anno_iter:
            layermap = self.layermap[anno.layer]
            if layermap is None:
                continue
            L, D = cfg.layername2LDT[layermap][0:2]
            if layermap == "bb_name":
                if not create:
                    pass
                    # print("ERROR: it's not allowed to flatten a black-box cell '{}'".format(cell.cell_name))
                name = anno.text.split("\n")
                if not cell.cell_name.startswith(name[0]):
                    pass
                    # print("ERROR: bb_name and cell_name of a black-box must be the same: '{}' != '{}'".\
                    #    format(name, cell.cell_name))
            if self.gds:
                GDS.content[level].append(
                    gbase.gds_annotation(xy=xy, string=anno.text, lay=L, texttype=D)
                )
            if self.nazca:
                # TODO: use 'layermap' instead of 'anno.layer'?:
                self.CELL.add_annotation(anno.text, anno.layer, xy)
        return None


    def add_instances(self, params):
        """Add instances to cell.

        Args:
            instance_iter (iterator): instance iterator
            level (int): hierarchy level
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        instance_iter = params.iters["instance"]
        level = params.parent_level
        hierarchy = params.hierarchy
        for inode, [x, y, a], flip in instance_iter:
            cell = inode.cell
            # only add references for instantiated cells
            if self.gds and (cell.instantiate or hierarchy == 'full'): # TODO: how to address hierarchy setting = 'full'
                GDS.content[level].append(
                    gdsmod.cell_reference(
                        [x, y],
                        cell.cell_name,
                        a,
                        flip=flip ^ inode.flip,
                        array=inode.array,
                        mag=inode.scale,
                    )
                )
                # TODO: check: for external gds files in cells named "load_gds",
                # the array setting needs to be transferred to the gds in this instance
                # if instantiate is False on "load_gds".
                if self.infolevel > 2:
                    print(
                        "{}add instance '{}' @ ({:.3f}, {:.3f}, {:.3f}), instance flip={}".format(
                            (level+1) * '. ',
                            cell.cell_name,
                            x, y, a,
                            inode.flip,
                            inode.array,
                        )
                    )
            if self.nazca and (cell.instantiate or hierarchy in ['self', 'full']):
                self.CELL.add_instance(inode, [x, y, a], flip)
        return None


    def add_gdsfiles(self, params):
        """Add gds instances to cell.

        Args:
            gdsfile_iter (iter): iterator over gdsfiles in the cell
            level (int): depth in the hierarchy

        Returns:
            None
        """
        for gdsinfo, [x, y, a], flip in params.iters['gdsfile']:
            if self.gds:
                filename, cellname, newcellname, layermap, cellmap, scale, strm = (
                    gdsinfo
                )
                GDS.content[params.level].append(
                    gdsmod.cell_reference([x, y], newcellname, a, mag=scale, flip=flip)
                )
                if newcellname not in self.gds_files:
                    self.gds_files[newcellname] = (
                        filename,
                        cellname,
                        layermap,
                        cellmap,
                        scale,
                        strm,
                    )
            if self.nazca:
                self.CELL.add_gds()
        return None


    def export_topcells(self, topcells):
        """Export topcells.

        Loop over each topcell in the topcells list.

        Args:
            topcells (list of Cell): list of topcells

        Returns:
            None
        """
        # initialize export formats with single output file across topcells.
        bbstr = ""
        if self.bblock:
            bbstr = "" # "_lib"
            bblock = BBlock(filebasename=self.filebasename + bbstr, outpath=self.dirname, bbpath=self.bbpath)
        if self.uPDK:
            uPDK = UPDK(filebasename=os.path.join(self.dirname, self.filebasename + bbstr))
        if self.gds:
            GDS.open(filebasename=os.path.join(self.dirname, self.filebasename + bbstr))
        if self.output == "file":
            mplt.ioff()
            # print('Switched of Matplotlib interactive output.')
        if self.nazca:
            self.CELL = ClsNazca(
                instantiate=self.instantiate,
                hierarchy=self.hierarchy,
                prefix=self.prefix,
                suffix=self.suffix,
            )

        cells_visited = set()
        for topcell in topcells:
            # initialize formats with a new output file per topcell
            if self.plt:  # separate plot for each topcell.
                PLT.open(topcell, title=self.title, show_pins=self.show_pins)
            if self.plotly:  # separate plot for each topcell.
                PLOTLY.open(topcell, title=self.title, show_pins=self.show_pins)
            if self.svg:  # separate plot for each topcell.
                SVG.open(self.filename)

            NL = Netlist()
            cell_iter = NL.celltree_iter(
                topcell,
                hierarchy=self.hierarchy,
                cells_visited=cells_visited,
                infolevel=self.infolevel,
            )
            for params in cell_iter:
                if params.cell_start:
                    if params.cell_create:
                        # take care of formats with cell hierarchy opening
                        if self.gds:
                            GDS.content[params.level].append(
                                gdsmod.cell_open(params.cell.cell_name)
                            )
                        if self.nazca:
                            self.CELL.open(params)
                    self.add_polygons(params=params)
                    self.add_polylines(params=params)
                    self.add_annotations(params=params)
                    self.add_instances(params=params)
                    self.add_gdsfiles(params=params)
                elif params.cell_close:
                    # take care of formats with cell hierarchy closing
                    if self.gds:
                        GDS.content[params.level].append(gdsmod.cell_close())
                        GDS.write(params.level)
                        del GDS.content[params.level]
                    if self.nazca:
                        self.CELL.close(params=params)

            cells_visited |= NL.cells_visited

            # Close exports that generate a file per topcell.
            if self.plt:
                PLT.close()
                if self.info:
                    print("...plotting plt")
                if self.output != "file" and cfg.plt_fig is None:
                    mplt.show()
            if self.svg:
                SVG.close()
                print("...Wrote file '{}'".format(SVG.filename))
            if self.plotly:
                if self.info:
                    print("...exorted to plotly")
                PLOTLY.close()
            if self.output is not None:
                try:
                    file = os.path.join(self.path, topcell.cell_paramsname) + ".png"
                except:
                    logger.warning("BB has no basename '{}'".format(topcell.cell_name))
                    file = os.path.join(self.path, topcell.cell_name) + ".png"
                print("export BB file: '{}'".format(file))
                mplt.savefig(file)
            if self.bblock:
                bblock.body(topcell)
            if self.uPDK:
                uPDK.body(topcell)
            if self.show_cells:
                print("Cells processed:", cells_visited)

        # Close exports that generate one file with all topcells.
        if self.bblock:
            bblock.footer()

        if self.uPDK:
            uPDK.footer()

        if self.gds:  
            # Save external GDS file based instances
            if self.infolevel > 0:
                print("----\nsave external gds instances")
            for newcellname, mapping in self.gds_files.items():
                filename, cellname, layermap, cellmap, scale, strm = mapping
                g = gstream.GDSII_stream(filename, cellmap=cellmap, layermap=layermap)
                stream = g.GDSII_stream_cell(newcellname)
                GDS.outfile.write(stream)
                if self.infolevel > 0:
                    print("{}'{}'".format(tab * 2, filename))
            GDS.outfile.write(gdsmod.layout_close())
            GDS.outfile.close()
            os.replace(GDS.tempname, GDS.filename)
            print("...Wrote file '{}'".format(GDS.filename))

            # Copy libraries of whitebox files and make whitebox mapping list:
            if self.submit:
                whitelibs = {}
                for Cname in cells_visited:
                    try:
                        C = cfg.cellnames[Cname]
                    except:
                        continue  # TODO: cell visited not on cellnames. Probably an export_gds loop without clear=False
                    if hasattr(C, "whitebox"):
                        whitebox = C.whitebox
                        pgpfile = whitebox.get("file", "pgp file entry not found")
                        md5 = whitebox.get("md5", "md5sum entry not found")
                        whitename = whitebox.get("cellname", "white cellname entry not found")
                        whitelibs[Cname] = (pgpfile, md5, whitename)
                whitelog = GDS.filename + ".ipblocks.csv"

                if whitelibs != {} or os.path.isfile(whitelog):
                    # Only create a whitelibs files if needed
                    with open(whitelog, "w") as F:
                        F.write("#cellname,white_cellname,encrypted_filename,md5_gds\n")
                        for name, (pgpfile, md5, whitename) in whitelibs.items():
                            F.write(f'"{name}","{whitename}","{os.path.basename(pgpfile)}",{md5}\n')
                            if not os.path.isfile(pgpfile):
                                logger.critical(
                                    f"missing ipblock library file as referenced in cell '{name}': {pgpfile}"
                                )
                            else:
                                self.md5 = True
                                copy2(pgpfile, self.dirname)
                    logger.info("exported: {}".format(whitelog))

            if self.ascii:
                ga = gstream.GDSII_stream(GDS.filename)
                ga.ASCII_write(GDS.filename + ".asc")

            if self.md5:
                md5sum = util.file2md5(os.path.join(GDS.filename), save=False)
                logger.info("md5sum of gds output: {}".format(md5sum))
                with open(GDS.filename + ".md5", "w") as F:
                    F.write(
                        "{}  {}".format(
                            md5sum, os.path.join(os.path.basename(GDS.filename))
                        )
                    )  # do not use a \n

            if util.isnotebook():
                display(
                    HTML(
                        '<pre>...<a href="{0}" target="_blank">{0}</a></pre>'.format(
                            GDS.filename
                        )
                    )
                )


    def generate_layout(self, topcells=None):
        """Internal wrapper function before exporting the layout.

        Create final topcells list and set final flags for export.

        Args:
            topcells (Cell | list of Cells): Cell(s) to export
                (default = None, which, exports the 'nazca' default gds cell)
            filename (str): gds output filename (default = 'nazca_export.gds')
            ascii (bool): export ascii version of gds (default = False)
            show_cells (bool): print exported cell names to stdout (default = False)
            gds (bool): export gds (default = True)
            flat (bool): export flat gds, i.e. no hierarchy (default = False)
            plt (bool): generate matplotlib based layout (default = False)
            clear (bool): clear mask layout between consecutive exports (default = True)
            title (str): title for the layout if the format allows for a title
            output (str): type of output stream (screen, file ...)
            path (str): output dir for saved Matplotlib plots (default = '')

        Returns:
            None
        """
        if isinstance(topcells, str):
            logger.ERROR(
                "You are trying to export a string instead of a cell. "
                "Did you mean: export_gds(filename='" + topcells + "')?"
            )
            return 0

        if self.infolevel > 0:
            print("gds:{}, nazca:{}, plt:{}, flat:{}".format(self.gds, self.nazca, self.plt, self.flat))
            print(f"hierarchy: '{self.hierarchy}'")
            print("Generate layout:")

        # close cells in reverse order.
        if not self.clear:  # keep default topcell open
            for cell in cfg.cells[1:-1:-1]:
                cell.close()
                if not cfg.solve_direct:
                    cfg.cells[0]._solve()  # solve (still open) default topcell
        else:
            for cell in cfg.cells[::-1]:
                cell.close()

        # construct a list of topcells
        if topcells is None:
            topcells = []
        if not topcells:  # 0 topcells
            topcells = [cfg.defaultcell]
        elif isinstance(topcells, Cell):  # 1 topcell
            topcells = [topcells]
        else:  # >1 topcells
            topcells = set(topcells + cfg.topcells)

        self.export_topcells(topcells)

        if self.clear:  # Create a new default cell.
            clear_layout()
            if self.infolevel > 0:
                print("Recreated topcell '{}'.".format(cfg.defaultcellname))

        return None


def clear_layout():
    """Remove all cell references to start a brand new layout.

    A new topcell 'nazca' will be created.

    Returns:
        None
    """
    for name, cell in cfg.cellnames.items():
        del cell
    cfg.cellnames = {}
    cfg.cells = []

    # cfg.cellnames.pop(cfg.defaultcellname) #avoid triggering cell reuse warning
    cfg.defaultcell = Cell(name=cfg.defaultcellname)
    cfg.topcell = cfg.defaultcell
    cfg.cp = cfg.defaultcell.org
    return None


def export_clear():
    """Clear the default topcell.

    Note that export_clear does not clear or delete any other cells.

    Returns:
        None
    """
    if cfg.cells:
        for cell in cfg.cells[::-1]:
            cell.close()
        cfg.cellnames.pop(cfg.defaultcellname)
        cfg.defaultcell = Cell(name=cfg.defaultcellname)
        cfg.cp = cfg.defaultcell.org
        # if self.infolevel:
        #    print("Recreated topcell '{}'.".format(cfg.defaultcellname))
    return None


# ==============================================================================
#
# ==============================================================================
def verify_topcells(topcells):
    """Verify if <topcells> has the correct format.

    Args:
        topcells (Cell | list of Cell): variable to check

    Returns:
        None: if topcells okay

    Exceptions:
        ValueError: topcells not valid
    """
    if topcells is None:
        return None
    elif isinstance(topcells, Cell):
        return None
    elif isinstance(topcells, list):
        for cell in topcells:
            if not isinstance(cell, Cell):
                raise ValueError("You need to provide a Cell or list of Cells for topcells in export.")
        return None
    raise ValueError(
        "You need to provide a Cell object or list of Cell objects"
        " for keyword 'topcells'. Instead got an object of type {}.".format(
            type(topcells)
        )
    )


def rebuild(
        cell,
        instantiate=True,
        flat=False,
        layermap=None,
        layermapmode=None,
        hierarchy='self',
        infolevel=0,
        clear=False,
    ):
    """Flatten, rename, relayer, reshape and/or filter a Nazca Cell object.

    The original cell remains unchanged.

    Args:
        cell (Cell): input cell(tree)
        instantiate (bool): instantiate the toplevel (True)
        flat (bool): flatten whole rebuild cell
        layermap (dict): layermap oldlayer:newlayer
        layermapmode (str): 'all': copy all layers,
            'none': copy only mapped layers, remove others.
        clear (bool): clear all cells after rebuild (default=False)

    Returns:
        Cell: rebuilt input cell
    """
    export = Export_layout()
    export.nazca = True
    export.flat = flat
    export.infolevel = infolevel
    export.setlayermap(layermap, layermapmode)
    export.instantiate = instantiate
    export.hierarchy = hierarchy
    export.clear = clear
    export.generate_layout(cell)
    return export.CELL.topcell


def celltree_iter(
    cell, level=0, position=None, flat=False, cells_visited=None, infolevel=0
):
    """Alias for a Nelist celltree_iter2 function."""
    return Netlist().celltree_iter2(
        cell=cell,
        position=position,
        flat=flat,
        cells_visited=cells_visited,
        infolevel=infolevel,
    )


def export(
    topcells=None,
    filename=None,
    gds=False,
    ascii=False,
    plt=False,
    svg=False,
    plotly=False,
    flat=False,
    infolevel=0,
    show_cells=False,
    info=True,  # progress info
    clear=None,
    title=None,
    output=None,
    path="",
    bb=False,
    bbpath="",
    uPDK=False,
    cellmap=None,
    layermap=None,
    layermapmode=None,
    md5=True,
    submit=False,
    show_pins=True,
    hierarchy='apply'
):
    """Export layout to gds file for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which, exports the 'nazca' default gds cell)
        filename (str): gds output filename (default = 'nazca_export.gds')
            The file name may include a path, in which case
            case the path becomes the base directory. Missing directories
            will be created.
        clear (bool): clear mask layout between consecutive exports (default = True)
        gds (bool): export gds (default = True)
        ascii (bool): export ascii version of gds (default = False)
        svg (bool): export svg file (default = False)
        plt (bool): generate matplotlib based layout (default = False)
        plotly (bool): generate plotly based layout (default = False)
        flat (bool): export flat gds, i.e. no hierarchy (default = False)
        infolevel (int): amount of debug info to stdout (default = 0)
        show_cells (bool): (default = False)
        info (bool): (default = True)
        title (str): title for the layout if the format allows for a title
        output (str): type of output stream (screen, file ...)
        path (str): output dir for saved Matplotlib plots (default = '')
        bb (bool): export as a library bb (default = False)
        bbpath (str): path to use to load a gds for bb=True generated modules (default="")
        cellmap (dict):
        show_pins (bool): show pins in the output (default=True)

    Returns:
        None
    """
    verify_topcells(topcells)
    if cfg.validation_layermapmode is not None:
        layermapmode = cfg.validation_layermapmode
    if cfg.validation_layermap is not None:
        layermap = cfg.validation_layermap    
    export = Export_layout(layermap=layermap, layermapmode=layermapmode)
    
    export.submit = submit  # keep this line before export.filename in this method!
    export.filename = filename
    export.ascii = ascii
    export.infolevel = infolevel
    export.show_cells = show_cells
    export.gds = gds
    hierarchy_allowed = ['flat', 'apply', 'self', 'full']
    if flat:
        hierarchy = 'flat'
    elif hierarchy in ['apply', 'full', 'flat']:
        hierarchy = hierarchy
    if gds and hierarchy == 'self':
        print(f"Can not use hierarchy '{hierarchy}' for gds output. Setting hierarchy = 'apply'.")
        hierarchy = 'apply'
    if plt and hierarchy in ['self', 'full', 'apply']:
        print(f"Can not use hierarchy '{hierarchy}' for plt output. Setting hierarchy = 'flat'.")
        hierarchy = 'flat'
    if svg and hierarchy in ['self', 'full', 'apply']:
        print(f"Can not use hierarchy '{hierarchy}' for svg output. Setting hierarchy = 'flat'.")
        hierarchy = 'flat'
    if plotly and hierarchy in ['self', 'full', 'apply']:
        print(f"Can not use hierarchy '{hierarchy}' for plotly output. Setting hierarchy = 'flat'.")
        hierarchy = 'flat'
    elif hierarchy not in hierarchy_allowed:
        raise ValueError(f"Value '{hierarchy}' for hierarchy unknown. Allowed values: {hierarchy_allowed}")
    export.hierarchy = hierarchy

    export.plt = plt
    if svg and not cfg.SVGWRITE:
        export.svg = False
        logger.error("Could not load module 'svgwrite'. Skipping svg export.")
        print("Could not load module 'svgwrite'. Skipping svg export.")
    else:
        export.svg = svg
    if plotly and not cfg.PLOTLY:
        export.plotly = False
        logger.error("Could not load module 'plotly'. Skipping plotly export.")
        print("Could not load module 'plotly'. Skipping plotly export.")
    else:
        export.plotly = plotly
    export.info = info
    export.title = title
    export.output = output
    export.path = path
    export.bblock = bb
    export.bbpath = bbpath
    export.uPDK = uPDK
    export.md5 = md5
    export.show_pins = show_pins
    if clear is None:
        export.clear = cfg.export_clear
    else:
        export.clear = clear
    if export.info:
        print("Starting layout export...")

    exported = False
    if export.gds:
        exported = True
        if info:
            print("...gds generation")
    if export.plt:
        exported = True
        if export.info:
            print("...matplotlib generation")
        if cfg.plt_fig is None:
            mplt.show()
    if export.svg:
        exported = True
        if info:
            print("...svg generation")
    if export.plotly:
        exported = True
        if info:
            print("...plotly generation")

    if not exported:
        msg = "No export type defined: plt, gds, svg, or plotly."
        logger.error(msg)
    export.generate_layout(topcells)
    return None


def export_svg(topcells=None, title=None, path="", **kwargs):
    """Export layout to svg format (scalable vector graphics).

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        title (str): title for the layout if the format allows for a title
        path (str): output dir for saved Matplotlib plots (default = '')

    Returns:
        None"""
    export(
        topcells,
        plt=False,
        plotly=False,
        gds=False,
        svg=True,
        info=False,
        title=title,
        **kwargs)


def export_plt(
    topcells=None,
    clear=None,
    title=None,
    output=None,
    path="",
    show_pins=True,
    hierarchy='flat',
    **kwargs,
):
    """Export layout with Matplotlib for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        clear (bool): clear mask layout between consecutive exports (default = True)
        title (str): title for the layout if the format allows for a title
        output (str): Matplotlib output stream (screen, file ...) (default = None -> screen)
        path (str): output dir for saved Matplotlib plots (default = '')

    Returns:
        None
    """
    export(
        topcells,
        plt=True,
        plotly=False,
        gds=False,
        info=False,
        clear=clear,
        title=title,
        output=output,
        path=path,
        show_pins=show_pins,
        hierarchy='flat',
        **kwargs
    )
export_png = export_plt


def export_plotly(
    topcells=None,
    clear=None,
    title=None,
    output=None,
    path="",
    show_pins=True,
    hierarchy='flat',
    **kwargs,
):
    """Export layout with Matplotlib for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        clear (bool): clear mask layout between consecutive exports (default = True)
        title (str): title for the layout if the format allows for a title
        output (str): Matplotlib output stream (screen, file ...) (default = None -> screen)
        path (str): output dir for saved Matplotlib plots (default = '')

    Returns:
        None
    """
    export(
        topcells,
        plt=False,
        plotly=True,
        gds=False,
        info=False,
        clear=clear,
        title=title,
        output=output,
        path=path,
        show_pins=show_pins,
        hierarchy='flat',
        **kwargs
    )


def export_gds(
    topcells=None,
    filename=None,
    flat=False,
    clear=None,
    bb=False,
    uPDK=False,
    md5=False,
    submit=False,
    bbpath="",
    **kwargs,
):
    """Export layout to gds file for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        filename (str): gds output filename (default = 'nazca_export.gds')
            The filename may include a path, in which case
            case the path becomes the base directory. Missing directories
            will be created.
        flat (bool): export flat gds, i.e. no hierarchy (default = False)
        clear (bool): clear mask layout between consecutive exports (default = True)
        bb (bool): Export design as a building block (default = False)'
        submit (bool): create a complete fileset for foundry submission (default=False)
        md5 (bool): create md5sum (default=False)
        bbpath (str): path to use to load a gds for bb=True generated modules (default="")

    Returns:
        None
    """
    verify_topcells(topcells)
    export(
        topcells,
        filename=filename,
        plt=False,
        gds=True,
        info=True,
        flat=flat,
        clear=clear,
        bb=bb,
        uPDK=uPDK,
        md5=md5,
        submit=submit,
        bbpath=bbpath,
        **kwargs
    )


def layout(
        layermap=None,
        layermapmode="all",
        infolevel=0,
        prefix='',
        suffix='_N',
    ):
    """Create a layout object for rebuilding cells.'

    Args:
        layermap (dict): layermap oldlayer: newlayer
        layermapmode (str): 'all': copy all layers,
            'none': copy only mapped layers, remove others.

    Returns:
        Export_layout: layout object
    """
    ly = Export_layout(
        layermap=layermap,
        layermapmode=layermapmode,
        prefix=prefix,
        suffix=suffix,
    )
    ly.nazca = True
    ly.CELL = ClsNazca(infolevel=infolevel, prefix=ly.prefix, suffix=ly.suffix,)
    return ly
