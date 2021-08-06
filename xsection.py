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
# @author: Ronald Broeke (c) 2016-2017
# @email: ronald.broeke@brightphotonics.eu
#-----------------------------------------------------------------------


"""Module defining xsection (cross section) functionality.

This includes the Material class.
"""


from itertools import count
import warnings

from functools import partial
import inspect
import numpy as np
import pandas as pd
from scipy import interpolate

# for 2D plots
from matplotlib.collections import PolyCollection
# for 3D plots
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import matplotlib.pyplot as plt

import nazca as nd
from nazca.simglobal import sim
from nazca import cfg
from nazca.solver_multilayerc import MultiLayerSlabSolverC


default_solverclass = MultiLayerSlabSolverC


nazca_folder = nd.__path__[0]
data_folder = nazca_folder + "/index_data"


class Material_util():
    """Class with utilities to define the effective index on material data.
    """

    def __init__():
        # TODO: calss must have init and selfin methods!
        pass


    def wl_check(self, wl, wlmin, wlmax, name):
        """Check if wavelength (wl) is in the allowed domain.

        Validates wl with the range of wlmin to wlmax.

        Args:
            wl (float): the wavelength to be checked in um.
            wlmin (float): the minimum wavelength of the data for interpolation.
            wlmax (float): the maximum wavelength of the data for interpolation.
            name (str): the name of the material data.

        Returns:
            bool: True if the wavelength check is passed.
        """
        # If wl is not in data range, put warning in log-file
        if (wl < wlmin) or (wl > wlmax):
            msg = f"Wavelength out of range in material '{name}' for wavelength {wl}: {wlmin} < wl < {wlmax} um."
            nd.main_logger(msg, 'warning')
        return True



    # TODO: not clear what the relation is between spline and data_spline
    def spline(self, wl, name, dat):
        """Generate an interpolated effective index based on material data.

        Generate an cubic spline if the input wavelength is in the data
        range. Evaluates the cubic spline at the input wavelength.

        Args:
            wl (float): wavelength to be interpolated in um.
            name (str): the name of the material data.
            dat (pd.DataFrame): data on known wavelength-Neff relation

        Returns:
            float:
        """
        wlmin = float(min(dat["wl"]))
        wlmax = float(max(dat["wl"]))
        if Material_util.wl_check(self, wl, wlmin, wlmax, name):
            Ipl = interpolate.CubicSpline(dat["wl"], dat["n"])
            return Ipl


    def data_spline(self, wl, name, n_type="real"):
        """Returns the effective index based on material model <name>.

        Reads the <name>.csv and interpolates the effective index for the
        given wavelength.

        Args:
            wl (float): wavelength to be interpolated in um.
            name (str): name of the material model .csv file.
            n_type (str): determines if the output is "real" or "complex" (if available).
                        Default value is "real", complex values slow calculations down.

        Returns:
            float: Interpolated effective index at wavelength wl.
        """
        filename = str(name + ".csv")
        path = data_folder + "/" + filename
        # name = mat_index
        # if not path.exists(): # TODO: Find alternative check.
        #     raise Exception(f"Material data for {name} not found.")
        # else:
        data = pd.read_csv(path)
        # Check if there is a string in the DataFrame
        strcheck = list(data.isin(["wl"]).any()[data.isin(["wl"]).any() == True].index)
        if not strcheck:
            # If there is no string, return the data
            return Material_util.spline(self, wl=wl, name=name, dat=data)(wl)
        else:
            # If there is a string, cut the data before the string and return the data
            row = int(list(data.isin(["wl"])["wl"][data.isin(["wl"])["wl"] == True].index)[0])
            new_dat = data[:row]
            if n_type == "complex":
                # If the Neff should be complex, gather the data from the
                # bottom of the file and add them.
                # Since this takes a long time, it is advised only if a complex
                # value for Neff is needed.
                comp_dat = [
                    float(data['n'][m])+float(data['n'][row+1+m])*1j for m in range(row)
                    ]
                new_dat['n'] = comp_dat
            return Material_util.spline(self, wl=wl, name=name, dat=new_dat)(wl)



class Material():
    """Define a material.

    The material can be a physical material or a virtual material.

    A virtual material is e.g. a slab layer stack and its index is function
    of wavelength and polarization.

    The material contains the Neff and color.
    """

    def __init__(self, Nmat=None, name='', rgb=()):
        """Initialize the material.

        Args:
            Nmat (float): refractive index of the material
            name (str): name of the material e.g. 'air'
            rgb (): color to display material. Value None give a random color

        Returns:
            None
        """
        if rgb == ():
            rgb = tuple(np.random.random(3))
        self.rgb = rgb
        self.name = name
        if callable(Nmat):
            self._Neff = Nmat
        elif isinstance(Nmat, float):
            self._Neff = Nmat
        elif Nmat is None:
            self._Neff = Nmat
            self._Neff = self.Neff()
        else:
            self._Neff = Nmat


    def Nmat(self, **kwargs):
        """
        Returns:
            float: materials its refractive index
        """
        return self.Neff(**kwargs)


    def Neff(self, **kwargs):
        """Get the index of the material.

        Returns:
            float: effective refractive index
        """
        self.wl = kwargs.pop('wl', sim.wl)
        self.pol = kwargs.pop('pol', sim.pol)
        if isinstance(self._Neff, float):
            return self._Neff
        elif callable(self._Neff):
            return self._Neff(wl=self.wl, pol=self.pol)
        else:
            try:
                self._Neff = Material_util.data_spline(self, wl=self.wl, name=self.name)
            except FileNotFoundError:
                raise Exception(
                    f"No effective index (Nmat) given for material {self.name}.\n"\
                    f"If index data should be used, check if '{self.name}.csv' is present\n"\
                    "in the {data_folder} directory."
                )
            return self._Neff

    def color(self, rgb):
        """Assign a color to the material.

        Args:
            rgb (tuple): (red, green, blue), values in  [0, 1]

        Returns:
            None
        """
        self.rgb = rgb


class Stack():
    """Class to define layer stacks.

    A Stack is a list of layers as [(material, width), ...], bottom up with
    the width (or height) in um.

    Note that the material act like a "virtual material". The Stack
    can be assigned a slab solver having method Neff.

    A stack has an etch function to optionally remove layers from the
    top depending on a specified etch depth.
    """
    _ids = count(0)  # stack id

    def __init__(
        self,
        layers=None,
        etchdepth=0,
        background=None,
        name='',
        solver=None,
        view='topview',
    ):
        """Initialize a stack object.

        Args:
            layers (list): list of tuples as (material, width)
            solver (class | function | float): slab solver class derived from class SolverBase.

        Returns:
            None
        """
        self.id = next(self._ids)
        self.background = background
        self.set(layers, etchdepth, solver, name, view)


    def simplifylayers(self, layers):
        """Merge neighbouring layers with the same material.

        Returns:
            list: layers, a list of tuples as (material, width)
        """
        mat = []
        width = []
        for i, (m, w) in enumerate(layers):
            if i > 0 and m == mat[-1]:
                width[-1] += w
            else:
                mat.append(m)
                width.append(w)
        #print (mat, width)
        return list(zip(mat, width))


    def applyetch(self, layers, etchdepth=0):
        """Apply etch to the layer stack.

        This reduces layer thickness or even removes layers if the etchdepth is
        larger than the layer thickness. The etch material is replaced with
        background material.

        Args:
            layers (list): list of layer tuples (material, width)
            etchdepth (float): etch depth measured from top of material

        Returns:
            list: layers after etch, list of tuples as (material, width)
        """
        if self.view != 'sideview' and etchdepth != 0:
            raise Exception(f"Etching requires a sideview Stack. Instead view='{self.view}'.")
        mat, width = map(list, zip(*layers))
        stack = reversed(width)
        drop = 0
        n = len(width)-1
        for i, w in enumerate(stack):
            #print (i, w, etchdepth)
            if w > etchdepth:
                width[n-i] = w-etchdepth
                break
            else:
                drop += 1
                etchdepth -= w
        if drop > 0 :
            width = width[:-drop]
            mat = mat[:-drop]
        if etchdepth > 0: # fill etch with background
            width.append(etchdepth)
            mat.append(self.background)
        return list(zip(mat, width))


    def set(
        self,
        layers=None,
        etchdepth=0,
        name=None,
        solver=None,
        view='',
    ):
        """Set the layer stack parameters.

        Args:
            layers (list): list of tuples (material, thickness)
            etchdepth (float): etchdepth into the stack measure from the top layer
            function (function): index function of the stack
            name (str): name of the stack
            solver (class, function, float): optional solver to obtain Neff for the stack.
            view (str): 'topview' or 'sideview'

        Returns:
            None
        """
        self.view = view
        layers = self.applyetch(layers, etchdepth)
        layers = self.simplifylayers(layers)
        self.layers = layers
        self.materials, self.widths = map(list, zip(*layers))  # map: need a list a later, not a tuple
        self.sumwidth = sum(self.widths)
        self.etchdepth = etchdepth
        self.solver = solver

        if name is None:
            self.name = 'stack-'.join(str(self._ids))
        else:
            self.name = name


    def reset(self):
        self._ids = count(0)


    def count(self):
        return self._ids


    def __str__(self):
        output = "name: '{}'\n".format(self.name)
        for mat, width in self.layers[::-1]:
            output += '  {}\t{:.3f}\n'.format(mat.name, width)
        return output


    def Neff(self, wl=None, pol=None, mode=None):
        """Add Neff to the stack to act like a virtual material."""
        if self.solver is None:
            nd.main_logger("No solver assigned to the stack to calculate Neff.", "error")
            return 1.0

        if hasattr(self.solver, "Neff"):
            self.solver.geometry(layers=self.layers)
            return self.solver.Neff(wl=wl, pol=pol, mode=mode)
        elif callable(self.solver):
            return self.solver(wl=wl, pol=pol, mode=mode)
        elif isinstance(self.solver, float):
            return self.solver
        else:
            nd.main_logger("Solver type not recognized {solver}.", "error")



class Xsection():
    """Define a cross sectional waveguide structure (XS).

    A structure refers to a type of optical waveguide or metal line.

    The xsection can hold multiple properties::

       - A table of (gds) mask layers needed to draw the xsection.
       - A set of Stack objects that describe the xsection's physical geometry in 1D or 2D.

    Note that for e.g. mode simulation gds layers are not needed and vise versa.

    The XS xsection waveguide geometry in this class is defined as a horizontal
    stack (hstack) of one or more 'vertical' layer stacks (vstack). This description
    can be used for 1D and 2D EIM.

    Example::

        three vstacks in a hstack

    Each vstack has one or more layers, and each layer has a material assigned to it.
    A vstack stack has a width attribute to express its with in the hstack.

    A function or mode solver can be assigned to materials, vstacks and hstacks
    to provide or calculate properties like the refractive index, mode profiles, etc.
    """

    def __init__(self, background=None, layers=None, name='', origin='', description=''):
        """Construct a xsection.

        Args:
            background ():
            layers ():
            name (str): xsection name
        """
        #waveguide simulation properties
        self.vstack = []
        self.hstack = []
        self.name = name
        self._stack_ids = count(0)
        self.layers = layers #epi
        self.background = background

        # waveguiding layout properties:
        self.os = 0 # straight-bend offset
        self.width = None # waveguide-width
        self.radius = None # radius of curvature
        self.minimum_radius = None # minimum radius in um
        self.mask_layers = pd.DataFrame() # for storing mask layer information
        self.symmetry = True # if waveguide in the xs are asymmetric
        self.pinstyle = None # visualistion of the xs pin in a layout
        self.origin = origin # origin of the Xsection for documentation
        self.description = description # description of the Xsection for documentation


    def background(self, mat):
        """Set the background material of the xsection.

        Returns:
            None
        """
        self.background = mat


    def reset(self):
        """Empty the xsection structure vstack and hstack.

        Returns:
            None
        """
        self.vstack = []
        self.hstack = []


    def stack_iter(self):
        """Returns:
            Interator over the stacks
        """
        for stack in self.stack:
            yield stack


    def levelhstack(self):
        """Make all vertical layers the same height.

        Fill the difference with background material.

        Returns:
            None
        """
        vwidths = []
        for materials, hwidth in self.hstack.layers:
            if hasattr(materials, 'sumwidth'):
                vwidths.append(materials.sumwidth)
        maxvwidth = max(vwidths)

        for vmat, hwidth in self.hstack.layers:
            if hasattr(vmat, 'sumwidth'):
                if vmat.sumwidth < maxvwidth:
                    if vmat.layers[-1][0] == self.background:
                        vmat.layers[:-1].append(
                            (vmat.layers[-1][0], maxvwidth - vmat.sumwidth - vmat.layers[-1][1])
                        )
                    else:
                        vmat.layers[:-1].append(
                            (self.background, maxwidth - vmat.layers.sumwidth)
                        )
        # TODO: reset the solver vstacks


    def add_vstack(
        self,
        layers=None,
        etchdepth=0,
        name='',
        solver=None,
    ):
        """Add a 'sideview' layer stack to the xsection.

        Args:
            layers (list): list of tuples (material, thickness)
            etchdepth (float): etchdepth into the stack measure from the top layer
            name (str): name of the stack
            function (function): index function of the stack returning a float
            solver (class): slab solver class to use for the stack

        Returns:
            Stack: object describing stack layers
        """
        #materials, widths = zip(*layers)
        if layers is None:
           layers = self.layers

        #TODO: check is name already exists.
        vstack = Stack(
            layers = layers,
            etchdepth = etchdepth,
            background = self.background,
            name = name,
            view='sideview',
            solver=solver,
        )
        self.vstack.append(vstack)
        vstack_id = len(self.vstack)
        return vstack


    def add_hstack(
        self,
        layers,
        name='',
        function=None,
        solver=None,
    ):
        """Define the hstack of the xsection.

        Args:
            layers (list): list of tuples (material, thickness)
            name (str): name of the stack
            function (function): index function of the stack returning a float
            solver (class): slab solver class to use for the stack

        Returns:
            Stack: horizontal stack
        """
        if layers is not None:
            self.hlayers = layers

        self.hstack = Stack(
            layers=layers,
            etchdepth=0,
            background=None,
            name=name,
            solver=solver,
            view='topview',
        )
        self.levelhstack()
        return self.hstack


    def Neff(self, wl=sim.wl, pol=sim.pol, radius=0.0, **kwargs):
        """Return the effective index of the xsection.

        Args:
            wl (float): Wavelength in um.
            pol (int): Polarization. 0 for "TE" or 1 for "TM".
            radius (float): Bending radius at the center of the waveguide.
                Default is 0.0 (straight waveguide).

        Returns:
            float: effective refractive index
        """
        R = kwargs.pop("R", None)
        if R is not None and radius == 0.0:
            radius = R
        return self.solver.Neff(wl=wl, pol=pol, radius=radius, **kwargs)


    def maxNeff(self):
        """Obtain maximum refractive index of all material in the xsection.

        Returns:
            float: maximum refractive index of all material in the xsection
        """
        Neffs = []
        stack_ids = [a.id for a in self.hstack['vstacks']]
        for v in self.stack:
             if v.id in stack_ids:
                 for mat in v.materials:
                     Neffs.append(mat.Neff(R=sim.R, wl=sim.wl))
        return max(Neffs)


    def minNeff(self):
        """Obtain minimum refractive index of all material in the xsection.

        Returns:
            float: minimum refractive index of all material in the xsection
        """
        Neffs = []
        stack_ids = [a.id for a in self.hstack['vstacks']]
        for v in self.stack:
             if v.id in stack_ids:
                 for mat in v.materials:
                     Neffs.append(mat.Neff(R=sim.R, wl=sim.wl))
        return min(Neffs)


    #def __str__(self):


    def info(self):
        """Create a string with information on the xsection.

        For printing the settings of all stacks and hstacks.

        Returns:
            str: information on the xsection
        """
        info = ''

        stack_defined = self.layers is not None
        if stack_defined:
            info = "XS name: '{}'\n".format(self.name)
            info += '* epi:\n'
            for mat, w in self.layers[::-1]:
                info += '  {}\t{:.3f}\n'.format(mat.name, w)
        else:
            return None
            #warnings.warn("No layers defined for xsection '{}':".\
            #    format(self.name), stacklevel=2)

        if self.hstack:
            info += '\n* hstack:\n'
            for i, (vstack, w) in enumerate(zip(self.hstack['vstacks'], self.hstack['widths'])):
                info += "'{}'({:.2f})".format(vstack.name, w)
                if i < len(self.hstack['vstacks'])-1:
                    info += ' | '
            info += '\n'
            stack_ids = [a.id for a in self.hstack['vstacks']]
        else:
            info += 'No hstacks defined yet. Use add_hstack().\n'

        return info


    # @property
    # def index(self):
    #     return self._IM
    # @index.setter
    # def index(self, IM):
    #     # TODO: make standard functions from index model visible as Structure class level?
    #     # RB: maybe just run all index function through XS.index.funcs()
    #     IM.xs = self  # register the xs to the solver

    #     self._IM = IM
    #     # self.Nridge = self._IM.Nridge
    #     # self.Nr = self._IM.Nridge
    #     # self.Nbg = self._IM.Nbg
    #     # self.Nsub = self._IM.Nsub
    #     # self.dNradius = self._IM.dNradius
    #     # self.Neff = partial(self._IM.Neff, width=self.width)
    #     return self._IM


    def showGuide2D(self, draw_index=False, ax=None):
        """Display a 2D representation of xsection.

        Returns:
            None
        """

        no_stack_defined = not self.hstack
        if no_stack_defined:
            warnings.warn("Warning: No layerstack defined to 2D draw xsection '{}'.\n Use help(nazca.Xsection) for options".\
                format(self.name), stacklevel=2)
            return None

        cfg.formatplot()
        if ax is None:
            fig, ax = plt.subplots()
            plotHere = True
        else:
            plotHere = False

        x0 = 0.0
        ymax = []
        w = self.hstack['widths']
        vertex = []
        colors = []
        for i, (width,stack) in enumerate(zip(self.hstack['widths'],self.hstack['vstacks'])):
            y0 = 0
            for j, (h, mat) in enumerate(zip(stack.widths, stack.materials)):
                #print (t, thick, mat)
                if draw_index:
                     intens = (mat.Neff()-self.minNeff()) / (self.maxNeff()-self.minNeff())
                     col = (intens, intens, intens)
                     #if col==(0,0,0): #do not draw white materials (air)
                     #     break
                else:
                    col = stack.materials[j].rgb
                    #if col==(1,1,1): #do not draw white materials (air)
                    #     break

                y1 = y0+h
                x1 = x0+w[i]
                x = np.array([x0, x1, x1, x0])
                y = np.array([y0, y0, y1, y1])
                vertex.append(np.swapaxes([x,y], 0, 1))
                colors.append(col)
                y0 = y1
                ymax.append(y0)
            x0 += width
        coll = PolyCollection(vertex, facecolors=colors, edgecolors='none')
        ax.add_collection(coll)
        #ax.autoscale_view()
        ax.set_xlim([0,x0])
        ax.set_ylim([0,max(ymax)])
        ax.set_xlabel('x [um]')
        ax.set_ylabel('y [um]')
        ax.set_title(self.name)
        #fig.colorbar(coll, ax=ax)
        if plotHere:
            plt.show()


    def showGuide3D(self, draw_index=False):
        """Display a 2D representation of xsection.

        Returns:
            None
        """
        cfg.formatplot()
        ax = a3.Axes3D(plt.figure())
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_zlabel('y')
        ax.set_title(self.name)

        z0 = 1.0
        x0 = 0.0
        ymax = []
        #w = [1.0] + self.hstack['widths'] + [1.0]
        try:
            w = self.hstack['widths']
        except:
            print("Warning: No layerstack defined to 3D draw xsection '{}'. Use 'help(nazca.Xsection)' for options.".format(self.name))
            return None
        for i, (width, stack) in enumerate(zip(self.hstack['widths'], self.hstack['vstacks'])):
            y0 = 0
            for j, (h, mat) in enumerate(zip(stack.widths, stack.materials)):
                if draw_index:
                    intens = (mat.Neff() - self._minNeff()) / (self._maxNeff()-self._minNeff())
                    col = (intens, intens, intens)
                    if col==(0, 0, 0):  # do not draw black materials
                        break
                else:
                    col = stack.materials[j].rgb
                    if col==(1, 1, 1):  # do not draw white materials (air)
                        break
                vertex = []
                y1 = y0+h
                x1 = x0+w[i]
                #front
                x = [x0, x1, x1, x0]
                y = [0,  0,  0,  0]
                z = [y0, y0, y1, y1]
                vertex.append(list(zip(x,y,z)))
                #back
                x = [x0, x1, x1, x0]
                y = [z0, z0, z0, z0]
                z = [y0, y0, y1, y1]
                #right
                vertex.append(list(zip(x,y,z)))
                x = [x1, x1, x1, x1]
                y = [z0,  0,  0,  z0]
                z = [y0, y0, y1, y1]
                #left
                vertex.append(list(zip(x,y,z)))
                x = [x0, x0, x0, x0]
                y = [z0,  0,  0, z0]
                z = [y0, y0, y1, y1]
                vertex.append(list(zip(x,y,z)))
                #top
                x = [x0, x1, x1, x0]
                y = [z0,  z0,  0,  0]
                z = [y1, y1, y1, y1]
                vertex.append(list(zip(x,y,z)))
                #bot
                x = [x0, x1, x1, x0]
                y = [z0,  z0,  0,  0]
                z = [y0, y0, y0, y0]
                vertex.append(list(zip(x,y,z)))

                y0 = y1
                ymax.append(y0)
                q = a3.art3d.Poly3DCollection(vertex)
                q.set_color(colors.rgb2hex(col))
                q.set_edgecolor('k')
                ax.add_collection3d(q)
            x0 += width
        ax.set_xlim([0,x0])
        ax.set_ylim([0,z0])
        ax.set_zlim([0,max(ymax)])
        print (x0, z0, max(ymax))
        plt.show()

# end class Xsection


def modalplot(buildXS, width, modes=None):
    XS = buildXS(max(width))
    N2 = []
    if modes is None:
        modes = XS.solver.maxModes
    for m in range(modes):
        N = []
        for w in width:
            XS = buildXS(w)
            XS.solver.mode = m
            n = XS.solver.Neff()
            if n<=0:
                n = np.nan
            N.append(n)
        N2.append(N)
    #return N2
    for i, n in enumerate(N2):
        plt.plot(w, n, '-')
    plt.xlabel('width [um]')
    plt.ylabel('Neff')


def modaldata(XSfunc, width, modes=None):
    """Calculate Neff against waveguide <width> for <modes> number of modes.

    Args:
        XSfunction (function): function F(width) to calculate Neff
        width (list of float): waveguide width in um
        modes (int): number of modes to calculate (if existing)

    Returns:
        DataFrame: index = width, columns = modes, data = Neff
    """
    N2 = [width]
    names = ['width']
    if modes is None:
        XS = XSfunc(max(width))
        modes = XS.solver.maxModes
    for m in range(modes):
        N = []
        names.append('mode'+ str(m))
        for w in width:
            XS = XSfunc(w)
            XS.solver.mode = m
            n = XS.solver.Neff()
            if n <= 0:
                n = np.nan
            N.append(n)
        N2.append(N)
        df = pd.DataFrame(np.array(N2).T, columns = names)
    return df


if __name__=="__main__":
    """
    TODO:
    Example code (Xsection):
        Retrieve a stack from demofab.
        Plot the xsection with materials and colors.
        Calculate the effective index of the xsection for various wavelengths.
    """
    from nazca.demofab import pdk_10_materials as mat
    from nazca.demofab.pdk_15_xsections import xsShallow


    # Define a wavelength range on which to study the xsection
    wl_range = np.linspace(1.3, 1.6, 11)
    sh_te_Neff = []
    bg_Neff = []
    rid_Neff = []
    inp_Neff = []
    q125_Neff = []
    print("--Nazca: xsection.py--")
    print("Calculate Neff as a function of wavelength")
    # Determine Neff (TE) for all defined wavelengths
    for n, wl in enumerate(wl_range):
        if n % 10 == 0:
            print(f"Step {n}/{len(wl_range)-1}")
        sh_te_Neff.append(xsShallow.Neff(wl=wl, pol=0))
        bg_Neff.append(xsShallow.stack[0].Neff(wl=wl, pol=0))
        rid_Neff.append(xsShallow.stack[1].Neff(wl=wl, pol=0))
        inp_Neff.append(mat.InP.Neff(wl=wl, pol=0))
        q125_Neff.append(mat.Q125.Neff(wl=wl, pol=0))

    # Plot the effective index as a function of wavelength.
    fig, ax = plt.subplots()
    sh_plot = ax.plot(wl_range, sh_te_Neff, label='Xsection: Shallow')
    rid_plot = ax.plot(wl_range, rid_Neff, label='Stack: Ridge')
    q125_plot = ax.plot(wl_range, q125_Neff, label='Material: Q125')
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Effective index (TE)")
    ax.set_title("Demofab Shallow xsection")
    plt.legend()
    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    bg_plot = ax.plot(wl_range, bg_Neff, label='Stack: Background')
    inp_plot = ax.plot(wl_range, inp_Neff, label='Material: InP')
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Effective index (TE)")
    ax.set_title("Demofab Shallow xsection")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Print out the information about the xsection
    print("\n" + "Stack information:")
    print(xsShallow.info())