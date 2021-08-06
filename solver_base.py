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
# (c) 2020-2021 Bright Photonics B.V.
# Author(s): Ronald Broeke
#
"""Interface class to access mode-solvers in a uniform way.

Use the SolverBase class as parent class to define placeholder/fallback
function for derived classes for solver interfaces.

Solver interfaces may require external solvers to be installed before usage.

Waveguide geometry is to be defined/stored in the Xsection module,
or provided by other means to the solver.
"""

import os
from sys import exit
from functools import wraps
import subprocess
from math import sin, cos, atan, exp, sqrt, pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import nazca as nd
from nazca.simglobal import sim
from nazca import cfg


class SolverBase():
    """Base class for adding and accessing Nazca (and external) mode solvers.

    This class provides a predefined set of template and fallback functions
    which are common for various solvers.

    Note that geometry information can be stored in Nazca's Xsection module.
    It is possible that new geometry approaches are needed for new solvers
    and these geometries can be also stored in Xsection.

    By design, the Xsection should be able to drive the solver interfaces. 
    Therefore, any object that has access to a Xsection object can do the same, 
    e.g. a pin or node.
    
    Electrical fields are stored in dict containers to combine data and
    metadate is a single object. The expected Efield labels are::
    
        solver_name (str)
        solver_type (?)
        layers (1D slab layer stack)
        Nsub (float), substrate index for the optical field
        view (str), 'sideview' or 'topview'
        wavelength (float)
        polarization (int), TE=0, TM=1
        power (float): field power
        data: modal field profile data in a DataFrame of (x, 0, 1, ...), 
            where columns 0, 1 contain the Efield of mode 0, 1, ...
        
    A 1D slab layer stack is described/built as a list of tuples::
        
        [(material1, thickness1), (material2, thicknes2), ...] 
        
    Materials can be of come in three flavors:
        - float
        - a callable returning a float
        - an object having a class method Neff() returning a float.         
    """
    def __init__(self, wl=None, pol=None, mode=0, view='topview', points=101):
        """Initialize the solver interface with geometry info from Xsection."""
        self.solves1D = False   # Ability to solve 1D applications
        self.solves2D = False   # Ability to solve 2D applications
        self.solvesStacks = False  # Ability to solve stacks
        self.add_params = False  # Additional parameters
        self.doRecalc = True
        self.type = "none"  # ?
        
        self.wl = None
        self.pol = None
        self.mode = None
        self._wl0 = None
        self._pol0 = None
        self._mode0 = None
        self.set(wl=wl, pol=pol, mode=mode, hard=True)
        self.points = points
        self.x0 = 0  # x-offset to apply to the data (solvers may have different x=0 locations).
        
        self._view0 = view
        if view is not None:
            self._view = view
        else:
            self._view = None

        self.modesflagcnt = 0  # mode counter flag. Only show message for 0, but this flag counter will count all occurances.
        self.fieldlabels = [  # expected Efield labels for reference
            'solver_name', 
            'solver_type', 
            'layers', 
            'Nsub', 
            'view',
            'wavelength',
            'polarization',
            'data',
        ]
                 
        
    def simsettings(func):
        """Decorator to set and reset wl, pol and mode settings for a specific "local" calculation.

        This decorator centralizes the handling of the solver parameters wl, pol and mode.
        
        The specific user settings in a call only need to be applied in
        the first call to the solver. Any sequential calls inside the solver
        use the same values.

        Some methods can be called by the user directly as well as internally.
        Such a method needs set wl, pol and mode utilizing this decorator.

        Returns:
            function: wrapper that executes func with the updated wl, pol and mode.
        """
        #@wraps(func)
        stack = []  # keep track of settings per call level for restoring
        def wrapper(self, *args, **kwargs):
            """Set wl, pol and mode based on inputs are fallback rules.
            
            Returns:
                output: result of wrapped function
            """
            wl = kwargs.pop('wl', self.wl)
            pol = kwargs.pop('pol', self.pol)
            mode = kwargs.pop('mode', self.mode)
            #update = (wl, pol, mode) != (self.wl, self.pol, self.mode) 
            #if update: 
            stack.append((self.wl, self.pol, self.mode))
            self.set(wl=wl, pol=pol, mode=mode, hard=False)
            result = func(self, wl=self.wl, pol=self.pol, mode=self.mode, *args, **kwargs)
            #if update:
            self.wl, self.pol, self.mode = stack.pop()
            return result
        return wrapper
    
    
    def set(self, wl=None, pol=None, mode=None, hard=True):
        """Set and update default modesolver settings.

        The default setting is to reuse the set values.        

        Args:
            wl (float): wavelength
            pol (int): polarization
            mode (int): mode
            hard (bool): If True define new reset values for wl, pol and mode
                if they are not None. Note that a set(hard=True) can be changed
                back to the Nazca system default with a reset(hard=True).
                
        Returns:
            self: updated solver object
        """
        if pol in ["TE", 'te', 0]:
            pol = 0
        elif pol in ["TM", "tm", 1]:
            pol = 1        
        
        if hard:  # store sticky reset values
            if wl is not None:
                self._wl0 = wl
            if pol is not None:
                self._pol0 = pol
            if mode is not None:
                self._mode0 = mode

        # fill parameters, with fallback options if needed:
        if wl is None:
            if self._wl0 is None:
                wl = sim.wl
            else: 
                wl = self._wl0
        if pol is None:
            if self._pol0 is None:
                pol = sim.pol
            else: 
                pol = self._pol0  
        if mode is None:
            if self._mode0 is None:
                mode = 0
            else:      
                mode = self._mode0
          
        if (wl, pol, mode) != (self.wl, self.pol, self.mode):
            self.wl, self.pol, self.mode = wl, pol, mode
            self.doRecalc = True  
        return self

    
    def geometry(self, **kwargs):
        """Placeholder for setting the geometry to solve.
        
        kwargs are whatever geometry info the solver needs.
        
        Returns:
            self: solver updated with the new geometry setting.
        """
        nd.main_logger(f"geometry() not defined in '{self.solverfile}'.", "warning")
        return self
        
    
    def _settings(self):
        """Obtain the current state of the solver settings.
        
        Can be added to output data of other functions, for example the field.
        
        Args:
            None.
        
        Returns:
            dict: Current settings for the calculated system.
        """
        # TODO: Do not do this here/ Use _localsimsettings! 
        if self.wl is None:
            wl = sim.wl
        else:
            wl = self.wl
        if self.pol is None:
            pol = sim.pol
        else:
            pol = self.pol
        
        layers = []
        for Nmat, width in self.layers:
            if callable(Nmat):
                Nval = Nmat(wl=wl, pol=pol)
            else:
                Nval = Nmat
            layers.append((Nval, width))
        Nsub = [self.Nsub(wl=wl, pol=pol) if callable(self.Nsub) else self.Nsub][0]
        settings = {}
        settings['solver_name'] = self.name
        settings['solver_type'] = self.type
        settings['layers'] = layers
        settings['Nsub'] = Nsub
        settings['view'] = self.view
        settings['wavelength'] = wl
        settings['polarization'] = "TM" if pol == 1 else "TE"
        return settings


    def setlayers(self, layers):
        """Register layers.

        Generate also self..materials list and self.widths lists for easy
        access to these values.

        Args:
            layers (list)

        Returns:
            list: layers
        """
        if layers is not None:
            self.layers = layers
            self.materials, self.widths = zip(*self.layers)
        else:
            self.layers = None
            self.materials, self.widths = None, None
        return self.layers


    def reset(self, hard=True):
        """Reset and update default modesolver settings.

        A soft/default reset will reset to the hard-set values.

        In a hard reset the wavelength and polarization are set to None.
        They will then default to the simglobal settings in functions.
        The mode is reset to 0.

        Args:
            hard (bool): hard reset (default=False) to Nazca system default.

        Returns:
            self: updated solver object
        """
        if hard:  # store absolute reset values.
            self._wl0 = None
            self._pol0 = None
            self._mode0 = 0
        self.wl = self._wl0
        self.pol = self._pol0
        self.mode = self._mode0
        return self


    @simsettings
    def modes(self, wl=None, pol=None):
        """Placeholder to return the number of guided modes.
        """
        if self.modesflagcnt == 0:
            nd.main_logger(f"Method modes() needs to defined in solver '{self.name}'. Returning 0 (this messsage only appears once).", "error")
            self.modesflagcnt += 1
        return 0


    def _polRot(self, pol):
        """Adjust polarization depending on the slab orientation

        Polarization TE is nominally defined as having the main E-vector in
        the horizontal direction (plane of the substrate) when looking at a 2D xsection.

        Translating a nominal 2D E-field to vertical and horizontal slabs
        results in the following effective physical polarization for the solver::

          - a TE solver for a vertically stacked slab; 'sideview'
          - a TM solver for a horizontally stacked slab; 'topview'

        Returns:
            int: polarization state for the slab simulation
        """
        if self._view == 'topview':
            pol = 1 if pol == 0 else 0
        return pol


    @simsettings
    def plotfield1D(self, wl=None, pol=None, mode=None, field=None, **kwargs):
        """Plot 1D field (x, E).

        If no field is provided, the field will be calculated via field1D().

        Args:
            field (list of DataFrame): optional set of fields.

        Returns:
            axis: Matplotlib plot object
        """
        if field is None:
            field = self.field1D(wl=self.wl, pol=self.pol, mode=self.mode)
            if field is False:
                return False 
        if not isinstance(field, list):
            field = [field]
        ax = None
        legend = []
        for fld in field:
            data = fld['data']
            legend.append(
                f"{fld['polarization']}{fld['mode']}, wl:{fld['wavelength']}, {fld['view'][:3]}\n"
            )
            if fld['polarization'] == 'TM':
                style = '--'
            else:
                style = '-'
            if ax is None:
                ax = data.plot(ax, style=style)
            else:
                data.plot(ax=ax, style=style)
        # TODO: title needs to adapt when multiple fields in "field" are provided
        ax.set_title(
            f"Efield 1D-slab solver:{fld['solver_name']}, norm:{fld['normalization']}"
        )
        ax.set_ylabel("Efield")
        ax.legend(legend, bbox_to_anchor=(1, 1), loc="upper left")
        plt.show()
        return ax
    field1Dplot = plotfield1D 


    def plotfield2D(self, **kwargs):
        """Contour plot of 2D mode solution.

        Replace this function in a child class when required for the solver.

        Returns:
            axis: Matplotlib plot object
        """
        X, Y, E = self.field2D(**kwargs)
        ax = plt.contourf(X, Y, E.real)
        plt.show()
        return ax
    showMode2D = plotfield2D  # TODO: let the Neff run via a solver interface, not directly in Xsections


    def field1Dpower(self, field):
        """Calculate power in modes in <field>
        
        Args:
            field (dict): Efield container.
            
        Returns:
             dict: {mode: power} power in the mode profiles in <field>
        """
        Edatain = field['data']
        x = Edatain.index.values
        powers = {}
        for mode in Edatain.columns:
            E = Edatain[mode].values           
            power = 0.5* sum([
                (E[i-1] + E[i])**2 * (x[i] - x[i-1]) 
                for i in range(1, len(x))
            ])            
            powers[mode] = power
        return powers
        
    
    def field1D_normpower(self, field, power=1.0):
        """Normalize all Efields in <field> to have power <power>.
                
        Args:
            field (dict): Efield container.
            power (float): value to scale the power. Default is 1.0.
        
        Returns:
            dict: normalized Efield container output data.
        """
        Edatain = field['data']
        x = Edatain.index.values
        Edataout = pd.DataFrame()
        powers = self.field1Dpower(field)
        for mode in Edatain.columns:
            if powers[mode] != 0:
                Edataout[mode] = Edatain[mode] * sqrt(power/powers[mode])
        Edataout.index = x
    
        fieldout = {}
        for label in field.keys():
            if label == 'data':
                fieldout[label] = Edataout
            else:
                fieldout[label] = field[label]  
        fieldout['normalization'] = power     
        return fieldout
   
     
    def field1D_normfield(self, field, fieldmax=1.0):
        """Normalize all Efields in <field> to have Efield max at <fieldmax>.
                
        Args:
            field (dict): Efield container.
            fieldmax (float): value to scale the fieldmax. Default is 1.0.
        
        Returns:
            dict: normalized Efield container output data.
        """
        Edatain = field['data']
        x = Edatain.index.values
        Edataout = pd.DataFrame()
        for mode in Edatain.columns:
            E = Edatain[mode].values           
            fieldinmax = max(E)
            Edataout[mode] = E * fieldmax/fieldinmax
            normalization = 'field'
        Edataout.index = x
    
        fieldout = {}
        for label in field.keys():
            if label == 'data':
                fieldout[label] = Edataout      
            else:
                fieldout[label] = field[label] 
        fieldout['normalization'] = 'field'
        return fieldout
   

    def field1D(self, wl=None, pol=None, mode=None):
        """Placeholder for getting 1D field and refractive index."""
        raise Exception(f"Method field1D() not defined yet in solver '{self.solverfile}'.")
     

    def field2D(self, wl=None, pol=None, mode=None):
        """Placeholder for getting 2D field and refractive index."""
        raise Exception(f"Method field2D() not defined yet in solver '{self.solverfile}'.")
    

    # def confinement1D(self):  # TODO: missing keywords
    #     """Placeholder for getting 1D field confinement."""
    #     raise Exception("confinment1D not defined yet")
      

    # def confinement2D(self):  # TODO: missing keywords
    #     """Placeholder for getting 2D field confinement."""
    #     raise Exception("confinment1D not defined yet")
      

    def Neff(self):  # TODO: missing keywords
        """Placeholder Neff"""
        raise Exception("Neff not defined yet")
      

    def PrintModeInfo(self, **kwargs):
        """Placeholder Neff"""
        nd.main_logger(f"PrintModeInfo() not defined in '{self.solverfile}'.", "warning")




