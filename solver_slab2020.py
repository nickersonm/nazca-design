#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
# 2017-2020 (c) Bright Photonics B.V.
# Author(s): Ronald Broeke
#-----------------------------------------------------------------------

"""
Simple slab solver for symmetrical 3 layer waveguides.

The solver can be topview or sideview, which determines the polarization
handling of the slabsolver. By definition: in a 2D view
-- with x pointing right and y pointing upward --
the component with the main E-field pointing in +x is TE2D
and the main E-field in the +y direction is TM2D.

"sideview":
    A stack with interfaces perpendicular to the y-direction
(boundaries parallel to the x-axis) has TE2D parallel to the boundaries,
hence, uses a TE solver for nominal TE polarization

"topview":
    A stack with interfaces perpendicular to the x-direction
(boundaries parallel to the y-axis) has TM2D parallel to the boundaries,
hence, uses a TM solver for nominal TE polarization.

"""

import os
from math import sin, cos, atan, exp, sqrt, pi
import numpy as np
import pandas as pd

import nazca as nd
from nazca.solver_base import SolverBase

thisfile = os.path.basename(__file__)

allowed_views = ['topview', 'sideview']


# TODO: add doRecalc concept to avoid repeating calculations.

class SlabMode(SolverBase):
    """Simple slab mode solver with sin, cos and exp elements.

    TE polarization for a 2D waveguide is defined as having the main E-field in
    the horizontal x-direction (horizontal). The y-directiony (vertical) is
    along the thickness of the guide. The z-direction (longitudinal) is the
    propagation direction.

    The default slab orientation is the 'topview' of a waveguide, where in the
    above definition a TM solver will be used for a nominal TE mode.
    For the "sideview" the nominal polarization is the same as the solver
    polarization.
    """

    def __init__(
        self,
        wl=None,
        pol=None,
        mode=0,
        name="SlabMode",
        layers=None,
        Nsub=None,
        view="topview",
        points=None,
    ):
        """Initialize a slab mode.

        Args:
            d (float):
            wl (float):
            mode (int): mode number, starting at 0.
            pol (int): polarization, TE->0, TM->1.
            view (str): 'sideview' for a ridge, or 'topview' for an EIM guide.
            points (int): number of point used in plotting fields.

        Returns:
            None
        """
        SolverBase.__init__(self)
        self.name = name
        self.type = "3-layer, symmetrical slab solver"
        self.solverfile = __file__

        self.slabFlag = 1
        self.solves1D = True
        self.solvesStacks = True
        self.properties = {}
        if points is not None:
            self.points = points

        self.dbmin = 1e-10  # accuracy for solving b
        self.set(wl=wl, pol=pol, mode=mode)
        if layers is not None:
            self.geometry(layers=layers, Nsub=Nsub)
        self.Nsub = Nsub  # substrate index

        if view is not None:
            self.view = view

        self.doRecalc = True
        self.flags = set()  # register warnings to issue them once where needed.

        # Field parameters used in mode profile:
        self.xs = ''
        self.V = None  # V-parameter
        self.b = None  # propagation contstant
        self.uu = None  # mode lateral wavevector under ridge
        self.cc = None  # mode lateral wavevector outside ridge
        self.maxModes = None  # max. number of modes based on V-parameter
        self.xclip = 0  # x value beyond which to clip the mode to 0
        self.Ge = 1.0  # factor for even slab modes under the ridge
        self.Go = 1.0  # factor for oneven slab modes under the ridge
        self.Pe = 0  # even mode constant
        self.Po = 0  # oneven mode constant
        self.sqrtPe = 0  # even mode constant
        self.qrtPo = 0  # oneven mode constant
        self.Je = 0  # even mode constant
        self.Jo = 0  # oneven mode constant

      

    def geometry(self, layers, Nsub=None):
        """Set geometry for the solver.

        Args:
            layers:

        Returns:
            self: solver with update geometry settings
        """
        self.doRecalc = True
        self.setlayers(layers=layers)
        if layers is None and self.layers is None:
            raise Exception(f"Layer stack 'layers' should be a symmetrical 3-layer stack for solver '{thisfile}'. Given: layers=None.")
        elif len(layers) == 1:
            if hasattr(self.Mat[1], "Neff"):
                self.Neff = self.materials[0].Neff
            else:
                self.Neff = self.materials[0]
        elif (
            len(layers) == 3
            and self.materials[0] == self.materials[2]  # symmetrical slab:  same material in layer 0 and 2
        ):
            self.d = self.widths[1]
            self.Nr = self.materials[1]
            self.Nbg = self.materials[0]
            if Nsub is not None:
                self.Nsub = Nsub
            else:
                self.Nsub = self.Nbg
        else:
            raise Exception(f"Layer stack 'layers' should be a symmetrical 3-layer stack for solver '{thisfile}'.")
        return self
    

    @property
    def view(self):
        return self._view
    @view.setter
    def view(self, val):
        if val is None:
            val = 'topview'
        if val in allowed_views:
            self._view = val
        else:
            raise Exception(f"Invalid view: '{val}'. Allowed views: {allowed_views}.")


    def setd(self, d):
        self.d = d
        self.doRecalc = True


    def _calcPe(self):
        """Get power normalization factor of an even mode.

        TODO: TE/TM?

        Pe = power in mode. Note: E/sqrt(Pe) gives normalized E-field
        Normalization can be clipped to [X0/2,-X0/2]

        Returns:
            float:
        """
        He = (self.Ge * cos(self.cc * self.d / 2.0))**2
        Pet = He / self.uu  # tail from -inf to -d/2 (edge slab)
        if self.xclip / 2.0 > self.d / 2.0:  # tail correction when clipping domain at xclip
            Pet = Pet - He*exp(self.uu * (self.d - self.xclip)) / self.uu
        Per = self.Ge**2 * (self.d / 2.0 + sin(self.cc * self.d) / (2.0 * self.cc)) # part under ridge
        Pe = Per + Pet
        return Pe


    def _calcPo(self):
        """Get power normalization of an odd mode.

        #TODO: TE/TM?

        Po = power in mode. Note: E/sqrt(Pe) gives normalized E-field
        Normalization can be clipped to [X0/2,-X0/2]

        Returns:
            float
        """
        Ho = (self.Go * sin(self.cc * self.d / 2.0))**2.0
        Pot = Ho / self.uu  # tail from -inf to -d/2 (edge slab)
        if self.xclip/ 2.0 > self.d/2.0: # tail correction when clipping domain at xclip
            Pot = Pot - Ho * exp(self.uu * (self.d - self.xclip)) / self.uu
        Por = self.Go**2 * (self.d/2.0 - sin(self.cc * self.d) / (2.0 * self.cc)) # part under ridge
        Po = Por + Pot
        return Po


    @SolverBase.simsettings
    def calculate_mode(self, wl=None, pol=None, mode=None):
        """Calculate all mode properties and store in class attributes.

        Set vertical stack mode numbers  to 0

        Args:
            wl (float):
            pol (int):    
            mode (int): mode number for horizontal stack

        Returns:
            bool: True if a mode does exist, False if not.
        """
        self.doRecalc = False
        if hasattr(self.Nr, "Neff"):
            Nr = self.Nr.Neff(wl=self.wl, pol=self.pol, mode=0)
        elif callable(self.Nr):
            Nr = self.Nr(wl=self.wl, pol=self.pol, mode=0)
        else:
            Nr = self.Nr

        if hasattr(self.Nbg, "Neff"):
            Nbg = self.Nbg.Neff(wl=self.wl, pol=self.pol, mode=0)
        elif callable(self.Nbg):
            Nbg = self.Nbg(wl=self.wl, pol=self.pol, mode=0)
        else:
            Nbg = self.Nbg
        assert Nbg > 0

        if hasattr(self.Nsub, "Neff"):
            Nsub = self.Nsub.Neff(wl=self.wl, pol=self.pol)
        elif callable(self.Nsub):
            Nsub = self.Nsub(wl=self.wl, pol=self.pol)
        else:
            Nsub = self.Nsub

        if self.slabFlag == 0:
            nd.logger.error("(Initialize slab properties before setting the mode.\n")

        modeinfo = f"(width={self.d}, wl={self.wl}, pol={self.pol})"
        try:
            V = (2 * pi * self.d / self.wl) * sqrt(Nr**2 - Nbg**2)
        except TypeError as E:
            #raise Exception(f"V-parameter values: {self.d=}, {self.wl=}, {Nr=}, {Nbg=}. {E}.") only python >= 3.8?
            raise Exception(f"V-parameter values: self.d={self.d}, self.wl={self.wl}, Nr={Nr}, Nbg={Nbg}. {E}.")
         
        self.maxModes = int(V / pi) + 1  # one polarization, not taking into account Nsub
       
        # TODO: is maxModes guided modes?
        if self.mode >= self.maxModes:
            flagid = (self.wl, self.pol, self.mode, self.d)
            if flagid not in self.flags:
                nd.main_logger(f"No guided mode #{self.mode} exists (max mode #{self.maxModes-1}), {modeinfo}", "warning")
                self.flags.add(flagid)
            self._Neff = 0
            self.b = 0
            self.V = V
            self.cc = 0
            self.uu = 0
            self.Pe = 0
            self.Po = 0
            self.sqrtPe = sqrt(self.Pe)
            self.sqrtPo = sqrt(self.Po)
            return False # mode does not exist.

        b = 0.8  # start value
        dbmin = self.dbmin
        db = (1.0 - b) - dbmin
        count = 0  # iteration counter
        countMax = 40  # maximum allowed iterations
        s = 1.0
        if self._polRot(self.pol) == 0:
            while db > dbmin and count < countMax: # abs to prevent tan from changing sign accidently
                count += 1
                value = V * sqrt(1.0 - b) - 2.0 * atan(sqrt(b / (1.0 - b)))
                #print(f"{value}, {db}, {b}")
                if value > self.mode * pi:
                    db /= 2.0
                    b += db
                    s = 2.0
                else:
                    db /= s
                    b -= db
                    if b < 0:
                        b = 0
            if count >= countMax:
                nd.main_logger(f"Reached max steps in '{thisfile}'. {modeinfo}", "error")
        elif self._polRot(self.pol) == 1:
            n2 = Nr**2 / Nbg**2
            while db > dbmin and count < countMax:  # abs to prevent tan from changing sign accidently
                count += 1
                value = V * sqrt(1.0 - b) - 2.0 * atan(n2 * sqrt(b / (1.0 - b)))
                if value > self.mode * pi:
                    db /= 2.0
                    b += db
                    s = 2.0
                else:
                    db /= s
                    b -= db
                    if b < 0:
                        b = 0
            if count >= countMax:
                nd.main_logger(f"Reached max steps in '{thisfile}'. {modeinfo}", "error")
        else:
            nd.main_logger(
                f"provided unknown polarization to InitMode() in '{thisfile}' {modeinfo}.\n",
                "error",
            )
        if count==countMax:
            nd.main_logger(
                f"Too many iterations: something went wrong in InitMode() in '{thisfile}'. {modeinfo}.\n",
                "error",
            )
        self._Neff = sqrt(b * (Nr**2 - Nbg**2) + Nbg**2)
        if self.view == "topview" and self._Neff < Nsub:
            #print(f"Mode cut-off: Neff < Nsub ({self._Neff:0.3f} < {Nsub:0.3f})")
            self._Neff = 0
            return False
        self.cc = (2 * pi / self.wl) * sqrt(Nr**2 - self._Neff**2)
        self.uu = (2 * pi / self.wl) * sqrt(self._Neff**2 - Nbg**2)
        self.b = b
        self.V = V
        # mode properties not related to x can be calculated once in advance:
        self.Pe = self._calcPe()
        self.Po = self._calcPo()
        self.sqrtPe = sqrt(self.Pe)
        self.sqrtPo = sqrt(self.Po)
        if self._polRot(self.pol) == 0:
            self.Je = self.Ge * cos(self.cc * self.d/2.0) / exp(-self.uu * self.d/2.0)  # continuous at waveguide edge TE
            self.Jo = self.Go * sin(self.cc * self.d/2.0) / exp(-self.uu * self.d/2.0)  # continuous at waveguide edge TE
        else:
            self.Je = self.Ge * cos(self.cc * self.d/2.0) / exp(-self.uu * self.d/2.0) * n2  # jump of n2 at waveguide edge TM 
            self.Jo = self.Go * sin(self.cc * self.d/2.0) / exp(-self.uu * self.d/2.0) * n2  # jump of n2 at waveguide edge TM
        #print ("(Nbg**2/Nridge**2)=" , (Nbg**2/Nridge**2))
        #print ("Ge*0.5*cc*sin(0.5*cc*d)=" , Ge*0.5*cc*sin(0.5*cc*d))
        #print ("Je*0.5*uu*exp(-0.5*uu*d)=", self.Je*0.5*self.uu*exp(-0.5*self.uu*self.d))

        self.lastmode = self.mode

        self.properties = {
            'view': self._view,
            'wl': self.wl,
            'pol': self.pol,
            'mode': self.mode,
            'width': self.d,
            'Neff': self._Neff,
            'Nridge': Nr,
            'Nbg': Nbg,
            'Nsub': Nsub,
            'Vparameter': self.V,
            'beta': self.b,
            'xs': self.xs,
            'cc': self.cc,
            'uu': self.uu,
            'Pe': self.Pe,
            'Po': self.Po,
            'sqrtPe': self.sqrtPe,
            'sqrtPo': self.sqrtPo,
            'Je': self.Je,
            'Jo': self.Jo,
        }

        return True


    def setXclip(self, x):
        """Set the clipping of the field.

        The field will be made zero outside the clipping.

        Args:
            x (float): clipping x position (half clipping width).

        Returns:
            None
        """
        self.xclip = x
        self.doRecalc = True


    def _getEe(self, x):
        """Return even mode E-field.

        Clips domain at x0 if x0>d

        Returns:
            float:
        """
        if self.doRecalc:
            self.calculate_mode()
        if x < -self.d/2.0:
            return self.Je * exp(self.uu*x) / self.sqrtPe**2
        elif x > self.d/2.0:
            return self.Je * exp(-self.uu*x) / self.sqrtPe**2
        else:
            return self.Ge * cos(self.cc*x) / self.sqrtPe**2


    def _getEo(self, x):
        """Return odd mode E-field.

        No clipping is applied.

        Args:
            x (float): x-position of the field.

        Returns:
            float:
        """
        if self.doRecalc:
            self.calculate_mode()
        if x < -self.d/2.0:
            return -self.Jo * exp(self.uu*x) / self.sqrtPo**2
        elif x > self.d/2.0:
            return self.Jo *  exp(-self.uu*x) / self.sqrtPo**2
        else:
            return self.Go * sin(self.cc*x) / self.sqrtPo**2


    @SolverBase.simsettings
    def getE(self, x, wl=None, pol=None, mode=None):
        """Return E field strength at x. Normalized per um.

        No clipping is applied.

        Args:
            x (float): x-position of the field.

        Returns:
           float: E field at x
        """
        if self.doRecalc:
            self.calculate_mode(wl=self.wl, pol=self.pol, mode=self.mode)
        if self.Pe == 0:
            return 0
        if self.xclip != 0 and abs(x) > self.xclip:
            E = 0
        else:
            self.absx01 = abs(self.xclip)
            if self.mode % 2 == 0:
                E = self._getEe(x)
            else:
                E = self._getEo(x)
        return E


    @SolverBase.simsettings
    def field1D(self, wl=None, pol=None, mode=[0]):
        """Calculate the the Efield of the slabmode.

        The field is power nomalized as sum(dx * E**2) = 1 .

        Args:
            wl (float): wavelength in um.
            pol (int): 0 for TE polarization, 1 for TM polarization.
            mode (list): List of mode numbers. [0, 1, ...]

        Returns:
            output (dict): Dictionary with parameter information and field data.
                Contains DataFrame 'data': field strength, index: "x" in um, columns: "E"
        """
        # TODO: allow for an array of modes
        if not isinstance(mode, list):
            mode = [mode]
        if self.layers is None:
            nd.main_logger("Can not calculate field1D: No layers set in the solver. Use geometry() method the set the layers first.")
            return False
        if self.points is None:
            self.points = 101
            
        xmin = -self.widths[0] - 0.5 * self.widths[1]
        xmax = self.widths[2] + 0.5 * self.widths[1]
        x = np.linspace(xmin, xmax, self.points)
        
        df = pd.DataFrame()
        for m in mode:
            df[m] = [
                self.getE(xi, wl=self.wl, pol=self.pol, mode=m) for xi in x
            ]
            if m == mode[0]:
                fieldinfo = self._settings()
        df.index = x
        df.index.name = "x"
        fieldinfo['mode'] = mode
        fieldinfo['normalization'] = "power"
        fieldinfo['data'] = df
        return self.field1D_normpower(fieldinfo)


# TODO: move to solver base
    def getEwidth(self, Efrac=0.5):
        """Find the waveguide width where field strenggth has fraction Efrac left w.r.t. center of guide.

        Args:
            Efrac (float): Fraction of the field relative to the field maximum. Efrac in [0, 1].

        Returns:
            float: width (this is a half width)
        """
        if self.doRecalc:
            self.calculate_mode()
        E0 = self.getE(0)
        Elimit = Efrac*E0
        x=0.0
        dx=0.3
        dxmin = 0.001
        t = 0
        tmax = 100
        while t < tmax and dxmin < dx:
            E = self.getE(x)
            if E > Elimit:
                x += dx
            else:
                x -= dx
                x = abs(x)
                dx /= 2.0
            t += 1
            #print(t, " ", x, " ", E, "\n")
            
        if t == tmax:
            nd.logger.warning ("Too many loops in getEwidth(). Returning 0 width.\n")
            return 0
        else:
            return x


    @SolverBase.simsettings
    def Neff(self, wl=None, pol=None, mode=None):
        """Get effective index of the slabmode.

        Note that geometry parameters must be initialed before a call to this method.

        Args:
            mode (int | list): mode number

        Returns:
            float: Neff
        """
        if self.doRecalc:
            neffs = []
            if not isinstance(self.mode, list):
                modes = [self.mode]
            else:
                modes = self.mode
            for m in modes:    
                self.calculate_mode(wl=self.wl, pol=self.pol, mode=m)
                neffs.append(self._Neff) 
            if len(neffs) == 1:
                return neffs[0]
            else:
                return neffs
        else:
            return self._Neff


    @SolverBase.simsettings
    def modes(self, wl=None, pol=None, **kwargs):
        """Calculate the number of guided modes in the slab.

        Note that geometry parameters must be initialed before a call to this method.

        Args:

        Returns:
            int: number of guided modes
        """
        if self.Nsub is None:
            # TODO: only show this message once per solver
            nd.main_logger(f"Nsub not set in '{thisfile}', setting it to 1.0", "error")
            Nsub = 1.0
        elif callable(self.Nsub):
            Nsub = self.Nsub(wl=self.wl, pol=self.pol)
        elif hasattr(self.Nsub, "Neff"):
            Nsub = self.Nsub.Neff(wl=self.wl, pol=self.pol)
        else:
            Nsub = self.Nsub

        if self.doRecalc:
            self.calculate_mode(wl=self.wl, pol=self.pol, mode=0)
            #print(f"{self.d}, {self.wl}, {self.pol}, {self.maxModes}")

        gm = 0
        for m in range(self.maxModes):
            if self.Neff(wl=self.wl, pol=self.pol, mode=m) > Nsub:
                gm += 1
            else:
                m = self.maxModes  # jump out of loop
        return gm


    @SolverBase.simsettings
    def PrintModeInfo(self, wl=None, pol=None, mode=0, width=None):
        """Print mode info to standard output.

        Args:
            mode (int): mode number, default=-1, which uses the mode set in the class attribute "mode".

        Returns:
            None
        """
        if isinstance(mode, list):
            nd.main_logger(f"Method PrintModeInfo can't take a list of modes {mode}, setting it to the first ({mode[0]}).")
            mode = mode[0]       

        if self.layers is None:
            nd.main_logger("Can not print mode info. No layers defined. Set a geometry first.", "error")
            return None
        if width is not None:
            self.d = width
        if self.doRecalc:
            self.calculate_mode(wl=self.wl, pol=self.pol, mode=mode)
        modes = self.modes(wl=self.wl, pol=self.pol)
        if mode > modes-1:
            print(f"Error mode > maxmodenumber: {mode} > {modes-1}")
            return None

        if mode != self.lastmode:
            self.doRecalc = True
        if self.doRecalc:
            self.calculate_mode(wl=self.wl, pol=self.pol, mode=mode)
        for varname, value in self.properties.items():
            if isinstance(value, float):
                print(f"{varname:12} = {value:6.5f}")
            elif isinstance(value, int):
                print(f"{varname:12} = {value:d}")
            else:
                print(f"{varname:12} = {value}")
        return None


    def PrintAllModes(self):
        """Print Neff for all guided modes in the slab.

        Returns:
            None
        """
        print("mode# Neff guided")
        for m in range(self.maxModes):
            if self._Neff > self.Nsub:
                guided = "yes"
            else:
                guided = "no"
            print(f"{m}: {self.getNeff(m)} {guided}")


    def calculate_overlap(self, mode, E, xclip=None):

        """Calculate the power overlap of the slabmode with a flat E field.

        The E-field is assumed constant across the slab mode.

        [-d/2, d/2] is domain of ridge.
        [-xclip/2, xclip/2] is domain where overlap is calculated, if x0>d
        [-inf, inf] is domain where overlap is calculated, if xclip < d/2 (waveguide width)

        Args:
            mode (int): mode number in arm, only used to determine even or odd mode.
            E (complex): incoming E-field density (in V/m / um ).

        Returns:
            float: field amplitude of mode due to overlap with E
        """
        if self.doRecalc:
            self.CalculateMode()

        if xclip is None:
            xclip = self.xclip # clip is half width
        a = -self.d / 2.0
        b = +self.d / 2.0

        if abs(xclip) < self.d / 2.0:
            nd.logger.warning (f"Wrong domain in def Overlap().\n...x0<d: d={self.d}, xclip={self.xclip}")
            return 0

        if mode % 2 == 0:
            # Field integral of ridge and tail sections for even modes:
            Ore = 2 * self.Ge * sin(self.cc * b) / self.cc
            Ote = 2 * self.Je * (exp(self.uu * a) - exp(-self.uu * xclip)) / self.uu
            #O = sqrt(2*xclip) * E * (Ore + Ote) / self.sqrtPe  # all power E*sqrt(x0)
            O =  E * (Ore + Ote) / self.sqrtPe  # all power E*sqrt(x0)
        else:
            Oro = 0
            Oto = 0
            O =  sqrt(2*xclip) * E * (Oro + Oto) / self.sqrtPo  # all power E*sqrt(x0)
        return O

# end class slabMode


if __name__ == "__main__":
    slab = SlabMode(points=401)
    slab.geometry(layers=[(1.5, 4.0), (3.2, 4.0), (1.5, 4.0)], Nsub=3.17)
    
    slab.PrintModeInfo(mode=0)
    slab.PrintModeInfo(mode=1)

    slab.plotfield1D(mode=1)
    slab.plotfield1D()
    slab.plotfield1D(wl=2.0)
    slab.plotfield1D(wl=5.0)
    slab.plotfield1D()
    slab.plotfield1D(wl=5.0, pol=1)
    slab.plotfield1D()
    slab.set(wl=5.0).plotfield1D()
    slab.plotfield1D(pol=1)
    slab.reset().plotfield1D()

    fields = []
    fields.append(slab.field1D(mode=0))
    fields.append(slab.field1D(mode=1))
    fields.append(slab.field1D(mode=2))
    fields.append(slab.set(pol=1).field1D(mode=0))
    fields.append(slab.field1D(mode=1))
    fields.append(slab.field1D(mode=2))
    slab.plotfield1D(field=fields)
   
    
    # slab.set(wl=1.3, pol=0, hard=True)
    
    # modelist = list(range(slab.modes()))
    # fieldsTM = slab.field1D(wl=1.55, pol=1, mode=modelist)
    # slab.plotfield1D(field=fieldsTM)
    # Neffs = slab.Neff(mode=modelist)
    # print(Neffs)
    
    # fieldsTE = slab.field1D(mode=[0, 1, 2])
    # slab.plotfield1D(field=fieldsTE)
    # # slab.plotfield1D(field=[fields[0], fields[1]])


    slab.PrintModeInfo()


