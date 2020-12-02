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
"""


from itertools import count
import warnings

import numpy as np
import pandas as pd
#for 3D plots
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import matplotlib.pyplot as plt

# 2D
from matplotlib.collections import PolyCollection

from .simglobal import sim
from . import slabsolver as slab
from . import cfg


class Material():
    """Define a material.
    
    The material can be a physical material or a virtual material. 
    
    A virtual material is e.g. a slab layer stack and its index is function
    of wavelength and polarization. 
    
    The material contains the Neff and color.
    """

    def __init__(self, Nmat=None, name='', rgb=None):
        """Initialize the material.
        
        Args:
            Nmat (float): refractive index of the material
            name (str): name of the material e.g. 'air'
            rgb (): color to display material. Value None give a random color    
        
        Returns:
            None
        """
        if rgb == None:
            rgb = tuple(np.random.random(3))
        self.rgb = rgb
        self.name = name
        if Nmat is None:
            self._Neff = 1.0
        elif isinstance(Nmat, float):
            self._Neff = Nmat
        else:
            self.Nmat = Nmat
            self._Neff = Nmat


    def Nmat(self, **kwargs):
        """
        Returns:
            float: materials its refractive index
        """
        return self._Neff


    def Neff(self, **kwargs):
        """Get the index of the material.
        
        Returns:
            float: effective refractive index
        """
        if isinstance(self._Neff, float):
            return self._Neff
        self._wl = kwargs.get('wl', sim.wl)
        self._pol = kwargs.get('pol', sim.pol)
        kwargs['wl'] = self._wl
        kwargs['pol'] = self._pol
        #print (kwargs)
        return self._Neff(**kwargs)


    def color(self, rgb):
        """Assign a color to the material.
        
        Args:
            rgb (tuple): (red, green, blue), values in  [0, 1]

        Returns:
            None
        """
        self.rgb = rgb  
        

class Stack():
    """Define slab layer stacks.

    The layers (optional) in stack are materials,
    each having a specific thickness in um.
    A refractive index function or solver (optional) can be assigned
    to the stack via a function or solver pointer.
    """
    _ids = count(0)

    def __init__(
        self, layers=None, 
        etchdepth=0, 
        background=None, 
        name='', 
        function=None, 
        solver=None, 
        view='topview'
    ):
        """"""""
        self.id = next(self._ids)
        self.background = background
        self.set(layers, etchdepth, function, solver, name, view)

        #self.name = function.__name__

    def simplifylayers(self, layers):
        """Merge neighbouring layers with the same material.

        Returns:
            list: list of tuples as (material, width)
        """
        mat = []
        width = []
        for i, (m, w) in enumerate(layers):
            if i > 0 and m == mat[-1]:
                width[-1] += w
            else:
                mat.append(m)
                width.append(w)
        #print (mat,width)
        return list(zip(mat, width))


    def applyetch(self, layers, etchdepth):
        """Apply etch to the layer stack.

        This reduces layer thickness or even removes layers.

        Args:
            layers (list): list of layer tuples (material, width)
            etchdepth (float): etch depth measured from top of material

        Returns:
            list: list of tuples as (material, width)
        """
        etchdepth
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


    def set(self, layers=None, etchdepth=0, function=None, solver=None, name=None, view=''):
        """Set the layer stack parameters.

        Args:
            layers (list): list of tuples (material, thickness)
            etchdepth (float): etchdepth into the stack measure from the top layer
            function (function): index function of the stack
            solver (function): solver object for the stack
            name (str): name of the stack
            view (str): 'topview' or 'sideview'

        Returns:
            None
        """

        layers = self.applyetch(layers, etchdepth)
        #if  self.background is not None:
        #    layers.append((self.background, 0.0))
        layers = self.simplifylayers(layers)
        self.layers = layers
        self.materials, self.widths = map(list, zip(*layers)) # map: need a list a later, not a tuple
        self.etchdepth = etchdepth
        self.function = function

        if name is None:
            self.name = 'stack-'.join(str(self._ids))
        else:
            self.name = name

        if solver is None:
            if len(layers) == 1:
                self._Neff = layers[0][0].Neff
                self.solver = None
            else:
                self.solver = slab.Slabsolver(name=self.name,
                    layers=layers, view=view)
                self._Neff = self.solver.Neff
        else:
            self.solver = solver
            solver.view = view
            self._Neff = solver.Neff


    def reset(self):
        self._ids = count(0)


    def count(self):
        return self._ids


    def Neff(self, **kwargs):
        """

        Returns:
            float: effective refractive index
        """
        return self._Neff(**kwargs)


    def __str__(self):
        output = "name: '{}'\n".format(self.name)
        for mat, w in self.layers[::-1]:
            output += '  {}\t{:.3f}\n'.format(mat.name, w)
        return output


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
        self.stack = []
        self.name = name
        self.hstack = {}
        self._stack_ids = count(0)
        self.layers = layers #epi
        self.background = background

        # waveguiding layout properties:
        self.os = 0 # straight-bend offset
        self.width = None # straight-bend offset
        self.radius = None # straight-bend offset
        self.minimum_radius = None # straight-bend offset
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
        self.stack = []
        self.hstack = {}


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
        h = []
        for s in self.hstack['vstacks']:
            h.append(sum(s.widths))
        m = max(h)
        #print ('h', h)
        for i, s in enumerate(self.hstack['vstacks']):
            if h[i]<m:
                if s.materials[-1] is self.background:
                    s.widths[-1] = m-sum(s.widths[:-1])
                else:
                    s.materials.append(self.background)
                    s.widths.append(m-sum(s.widths[:]))


    def add_vstack(
        self, 
        layers=None, 
        etchdepth=0, 
        name='', 
        function=None,
        solver=None
    ):
        """Add a 'sideview' layer stack to the xsection.

        Args:
            layers (list): list of tuples (material, thickness)
            etchdepth (float): etchdepth into the stack measure from the top layer
            name (str): name of the stack
            function (function): index function of the stack returning a float
            solver (function): solver object for the stack

        Returns:
            Stack: object describing stack layers
        """
        #materials, widths = zip(*layers)
        if layers is None:
           layers = self.layers

        #TODO: check is name already exists.
        stack = Stack(
            layers, 
            etchdepth, 
            self.background, 
            name, 
            function, 
            solver,
            view='sideview'
        )
        self.stack.append(stack)
        stack_id = len(self.stack)
        #return stack_id-1
        return stack


    def add_hstack(self, layers, name='', function=None, solver=None):
        """Define the hstack of the xsection.

        Args:
            layers (list): list of tuples (material, thickness)
            name (str): name of the stack
            function (function): index function of the stack returning a float
            solver (function): solver object for the stack

        Returns:
            None
        """
        if layers is not None:
            self.hlayers = layers
        sts, widths = zip(*self.hlayers)
        self.hstack = {
            'widths': widths,
            'vstacks': sts,
            'function': function,
            'solver': solver
        }

        self.solver = solver
        if solver is not None:
            self.solver.view = 'topview'
            self._Neff = self.solver.Neff
        else:
            self.solver = slab.Slabsolver(layers=self.hlayers, name=name, view = 'topview')
            self._Neff = self.solver.Neff

        self.levelhstack()


    def Neff(self, **kwargs):
        """Return the effective index of the xsection.

        Args:
            kwarga: wl, R, w
            
        Returns:
            float: effective refractive index
        """
        R = kwargs.pop('R', 0)
        w = kwargs.pop('w', 0)
        wl = kwargs.pop('wl', sim.wl)
        pol = kwargs.pop('pol', sim.pol)
        return self._Neff(R=R, wl=wl, pol=pol)

    def set_Neff(self,func):
        self._Neff=func


    def showSolvers(self):
        """Show assigned solvers.

        Returns:
            None
        """
        print ('Solvers:')
        #waveguide
        stack_ids = [a.id for a in self.hstack['vstacks']]
        print ('rib-order', stack_ids, ', function:', self.hstack['function'], ', solver:', self.hstack['solver'])
        #stack
        for v in self.stack:
             if v.id in stack_ids:
                 print ('rib', v.id, ', function:', v.function, ', solver:', v.solver)


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
        """Obtain minnimum refractive index of all material in the xsection.

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


#    def __levelStacks(self):
#        Heights = []
#        stack_ids = [a.id for a in self.hstack['vstacks']]
#        for v in self.stack:
#             if v.id in stack_ids:
#                h0 = 0
#                for h in v.widths:
#                    h0 += h
#                heights.append(h0)
#        maxHeight = max(height)
#        return max(Neffs)


#    def showStacks(self):
#        for i, v in enumerate(self.stack):
#            print ('rib', i, ', function:', v.function, ', solver:', v.solver)


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

        if self.stack:
            for num, stack in enumerate(self.stack):
                if stack.solver is None:
                    solvername = 'None'
                else:
                    solvername =  stack.solver.name
                if stack.function is None:
                    functionname = 'None'
                else:
                    functionname = 'yes'

                info = info + "\n* vstack no.{}:\n".format(num)
                info += stack.__str__()
                info += "  solver: '{}'\n  function: {}\n".format(solvername, functionname)
        else:
            info += 'No vstacks defined yet. Use add_vstack().\n'

        return info


    @property
    def index_model(self):
        return self._IM
    @index_model.setter
    def index_model(self, IM):
        # make standard functions from index model visible as Structure class level
        self._IM = IM
        self.Nrte = self._IM.Nrte
        self.Nrtm = self._IM.Nrtm
        self.Nbte = self._IM.Nbte
        self.Nbtm = self._IM.Nbtm
        self.Nste = self._IM.Nste
        self.Nstm = self._IM.Nstm
        self.dNRte = self._IM.dNRte
        self.dNRtm = self._IM.dNRtm
        self.Neff2Dte = self._IM.Neff2Dte
        self.Neff2Dtm = self._IM.Neff2Dtm
        self.oste = self._IM.SBOte
        self.ostm = self._IM.SBOtm
        return self._IM


    def Nrte(self, wl):
        """

        Returns:
            float: effective refractive index of the 'ridge'
        """
        return self._IM.Nrte(wl)


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
            if n<=0:
                n = np.nan
            N.append(n)
        N2.append(N)
        df = pd.DataFrame(np.array(N2).T, columns = names)
    return df
