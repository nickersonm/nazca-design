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
# @author: Ronald Broeke (c) 2016-2017
# @email: ronald.broeke@brightphotonics.eu
#
"""Multi-layer solver interface.

This module interfaces with an external slab solver program in C.
IO communication is via an input file and an output file.

Other solvers can be connected using the same concept. Depending on the solver
this may be extended to sockets or other communication channels.

"""

from sys import exit
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import nazca as nd
from nazca import cfg
from nazca.solver_base import SolverBase

# TODO: add doRecalc concept to avoid repeating calculations.

# path of the solver relative to the work dir or absolute
slabpath = os.path.dirname(__file__)
thisfile = os.path.basename(__file__)


class MultiLayerSlabSolverC(SolverBase):
    """Multi-layer slabmode c-solver interface.

    The last used values for wl and pol are stored in self.wl and self.pol.
    """

    def __init__(
        self,
        wl=None,
        pol=None,
        mode=0,
        name='Multi_Layer',
        layers=None,
        Nsub=None,
        points=100,
        view='topview',
    ):
        """Initialize the Nazca to slabc multilayer solver interface.

        Args:
          wl (float): wavelength in um
          pol (int): polarization of the mode (0 for TE or 1 for TM)
          mode (): optical mode
          name (): name of the solver
          layers (lsit): list of layer tuples (material, width)
          points (): point to calculate for a mode field
          view (): 'topview' or 'sideview'

        Returns:
            None
        """
        SolverBase.__init__(self)

        self.name = name
        self.type = "MLS (multi-layer-solver)"
        self.solverfile = __file__

        self.set(wl=wl, pol=pol, mode=mode)
        self.point = points
        self.geometry(layers=layers, Nsub=Nsub, view=view)
        self.plotNum = 0
        if name is None:
            self.solvername = "multi-layer-slab-solver"
        else:
            self.name = name


    def geometry(self, layers, Nsub=None, view='topview'):
        """Set the geomtery to solve.

        Returns:
            None
        """
        if layers is not None:
            self.setlayers(layers=layers)
            #raise Exception(f"Provide a valid layer list [<material>, <width>), ...]. given instead: '{layers}'")
            x0 = 0
            for n in range(len(self.widths) // 2):
                x0 += self.widths[n]    
            self.x0 = x0 + 0.5* self.widths[n]
            self.Nsub = Nsub
        if view in ['topview', 'sideview']:
            self.view = view
        else:
            raise Exception(f"View '{view}' not existing")


    @SolverBase.simsettings
    def _calculate_mode(self, wl=None, pol=None, mode=None, **kwargs):
        """Calculate the E-field and Neff of the mode by calling external solver.

        Args:
            wl (float):
            pol (int):

        Returns:
            DataFrame, float: field as (x: E), Neff
        """
        if self.view == 'topview':
            pol = self._polRot(self.pol)
        else:
            pol = self.pol
        if pol == 0:
            polstr = "TE"
        elif pol == 1:
            polstr = "TM"

        if self.materials is None:
            exit("No materials defined. Can not calculate field.")

        N = []
        for material in self.materials:
            if hasattr(material, "Neff"):
                N.append(material.Neff(wl=self.wl, pol=self.pol))
            elif callable(material):
                N.append(material(wl=self.wl, pol=self.pol))
            else:  # assuming a float index
                N.append(material)

        if len(self.materials) == 1:  # single layer stack
            field, Neff = None, N[0]

        else:  # create slabin.dat file for use by the solver
            slabin = \
f"""Wavelength (um) = {self.wl}
Polarisation = {polstr}
Mode order = {self.mode}
Number of layers = {len(N)}
Refractive indices = {" ".join(map(str, N))}
Thicknesses (um) = {" ".join(map(str, self.widths))}
Number of plot intervals = {self.points}
"""
            slabin_filename = "slabin.dat"
            with open(slabin_filename, "w") as text_file:
                text_file.write(slabin)

            # create slabout.dat, read data and use subprocess to get return value
            solver_exec = os.path.join(slabpath, "slab")
            proc = subprocess.Popen(
                os.path.join(slabpath, "slab"),
                stdout=subprocess.PIPE,
                shell=True
            )
            stdout_bytestr, err = proc.communicate()
            if err is not None:
                print(f"Slabsolver err: {err}")
            stdout = stdout_bytestr.decode("UTF-8")
            stdout = stdout.split("=")[-1]
            stdout = stdout.rstrip("\n")
            try:
                self._Neff = float(stdout)
                if self._Neff < self.Nsub:
                    self._Neff = 0
                    self._field = pd.DataFrame(index=['x'], columns=['E'])
                else:
                    self._field = pd.read_csv("slabout.dat", delimiter=" ", names=["x", "E"])
                    self._field['x'] -= self.x0
                    self._field.set_index("x", inplace=True)
            except:
                self._Neff = 0
                self._field = pd.DataFrame()
                nd.main_logger(
                    f"Solver module '{thisfile}' received message: {stdout}\nSee for more details input file '{slabin_filename}':\n{slabin}",
                    "warning"
                )
            field, Neff = self._field, self._Neff

        self.doRecalc = False
        return field, Neff


    def modes(self, wl=None, pol=None, **kwargs):
        """Calculate the number of guided modes in the slab.

        Note that geometry parameters must be initialed before a call to this method.

        Args:
            wl (float): wavelength in um.
            pol (int): 0 for TE polarization, 1 for TM polarization.

        Returns:
            int: number of guided modes
        """
        maxmodecnt = 20  # do not try to find more modes.
        if self.Nsub is None:
            # TODO: only show this message once per solver
            nd.main_logger(f"Nsub not set in '{thisfile}', setting it to 1.0", "error")
            Nsub = 1.0
        elif type(self.Nsub) in [float, int, np.int64, np.float64]:
            Nsub = self.Nsub
        elif callable(self.Nsub):
            Nsub = self.Nsub(wl=self.wl, pol=self.pol)
        elif hasattr(self.Nsub, "Neff"):
            Nsub = self.Nsub.Neff(wl=self.wl, pol=self.pol)
        else:
            nd.main_logger("Nsub not defined. Setting it to 1.0", "warning")
            Nsub = 1.0

        N = self.Neff(wl=wl, pol=pol, mode=0)
        modecnt = 0
        while N > Nsub and modecnt < maxmodecnt:
            modecnt += 1
            N = self.Neff(wl=wl, pol=pol, mode=modecnt)
        return modecnt
    

# TODO: move to /integrate with solver base plotfield1D:
    def _plot(self, title=None, rotate=False, **kwargs):
        """Plot the E-field.

        Args:
            title (str): plot title
            rotate (bool): rotate the plot by 90 degrees
            **kwargs ():

        Returns:
            MatPlotLib axis: axis object of the plot produced.
        """
        cfg.formatplot()

        if title is None:
            title = self.name

        field, Neff = self.getField()
        x = kwargs.pop("x", "x")
        y = kwargs.pop("y", "E")
        xlabel = "x [um]"
        ylabel = "E [field normalized]"
        clip_on = kwargs.get("clip_on", False)
        if rotate == True:
            x, y = y, x
            ylabel, xlabel = "y [um]", ylabel
        ax = kwargs.pop("ax", None)
        if ax is None:
            fig, ax = plt.subplots(figsize=cfg.modeplotsize)

        # ax = self.ax
        col = "b" if self.pol == 0 else "g"
        linestyles = ["-", "-", "-", "-"] * 8
        ls = linestyles[self.mode]

        Dc = np.cumsum(self.D)
        Dc = np.insert(Dc, 0, [0])

        if self.plotNum >= 0:
            ax.set_title(title)
            if rotate == False:
                ax.set_xlim([Dc[0], Dc[-1]])
            else:
                ax.set_ylim([Dc[0], Dc[-1]])
            for i, d in enumerate(Dc[:-1]):
                if rotate == False:
                    ax.axvspan(
                        Dc[i],
                        Dc[i + 1],
                        facecolor="m",
                        alpha=0.01 + 0.05 * (self.Mat[i].Neff() - 1),
                    )
                else:
                    ax.axhspan(
                        Dc[i],
                        Dc[i + 1],
                        facecolor="m",
                        alpha=0.01 + 0.05 * (self.Mat[i].Neff() - 1),
                    )
            [ax.spines[key].set_linewidth(2.0) for key in ax.spines.keys()]

        field.plot(x=x, y=y, ls=ls, color=col, clip_on=clip_on, ax=ax, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.tight_layout()

        self.plotNum += 1
        return ax


    def __str__(self):
        out = "Slab-solver settings:"
        out += "\n- slab name                            :'" + self.name + "'"
        out += "\n- type                                 : MLS (multi-layer-solver)"
        out += "\n- view                                 : " + self.view
        out += "\n- wavelength                 - wl [um] : " + str(self.wl)
        out += "\n- polarization               - pol     : " + str(self.pol)
        out += "\n- materials                  - M       : " + ", ".join(
            [m.name for m in self.Mat]
        )
        out += "\n- layer thickeness           - D [um]  : " + ", ".join(
            map(str, self.D)
        )
        out += "\n- refractive indices         - N       : " + ", ".join(
            map(str, [m.Neff() for m in self.Mat])
        )
        out += "\n- number of points for Field - points  : " + str(self.points)
        out += "\n"
        return out


    @SolverBase.simsettings
    def Neff(self, wl=None, pol=None, mode=None, **kwargs):
        """Calculated Neff of the slabmode.

        Args:
            wl (float): wavelength in um
            pol (int): polarization of the mode (0 for TE or 1 for TM)
            mode (list): list of mode numbers [0, 1, ...].

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
                field, neff = self._calculate_mode(wl=self.wl, pol=self.pol, mode=m, **kwargs)
                neffs.append(neff) 
            self._Neff = neffs
            if len(neffs) == 1:
                return neffs[0]
            else:
                return neffs
        else:
            return self._Neff


    @SolverBase.simsettings
    def field1D(self, wl=None, pol=None, mode=[0], **kwargs):
        """Calculate the the Efield of the slabmode.

        The field is normalized to Emax.

        Args:
            wl (float): wavelength in um.
            pol (int): 0 for TE polarization or 1 for TM polarization.
            mode (list): list of mode numbers [0, 1, ...].
        
        Returns:
            output (dict): Dictionary with parameter information and field data.
                Contains DataFrame 'data': field strength, index: "x" in um, columns: "mode"
        """
        if type(self.mode) is not list:
            mode = [mode]
        else:
            mode = self.mode
        df = pd.DataFrame()
        for m in mode:
            field, neff = self._calculate_mode(wl=self.wl, pol=self.pol, mode=m)
            df[m] = field['E']
            if m == mode[0]:
                df.index = [m for m in field.index]
                df.index.name = "x"
                fieldinfo = self._settings()            
        fieldinfo['mode'] = mode
        fieldinfo['normalization'] = "peak"
        fieldinfo['data'] = df
        return self.field1D_normpower(fieldinfo)        
      

if __name__ == "__main__":
    from nazca import xsection
    XS = nd.add_xsection('testXS')

    #Define materials. Use value or functions to set the refractive index:
    M1 = xsection.Material(Nmat=1.5, name='mat1', rgb=(0.0, 0.4, 0.9))
    M2 = xsection.Material(Nmat=3.2, name='mat2', rgb=(0.0, 0.8, 0.3))
    Mair = xsection.Material(Nmat=1.0, name='air', rgb=(0.95, 0.95, 1.0))

    wguide = 3.0
    hfilm = 0.6
    hsub = 1.0
    hclad = 1.0
    epi = [(M1, hsub), (M2, hfilm), (M1, hclad)]

    # test waveguide
    XS.layers = epi
    XS.background = Mair
    vs1 = XS.add_vstack(name='backgrnd', etchdepth=hclad + 0.2,)
    vs2 = XS.add_vstack(name='ridge',)
    hs = XS.add_hstack(
        name='guide',
        layers=[(vs1, 1.0), (vs2, wguide), (vs1, 1.0)],
    )
    
    # Calculate modes in waveguide + parasitic waveguide
    layers = [(2.5, 4.0), (3.5, 1.0), (2.5, 8.0), (3.5, 2.0), (2.5, 4.0),]
    mulsol = MultiLayerSlabSolverC(view="topview", points=1900)
    mulsol.geometry(layers=layers, Nsub=1.0)
    print("# of modes:", mulsol.modes(wl=1.55, pol=0))
    
    modes = mulsol.modes()
    fields = mulsol.set(wl=1.55, pol=0).field1D(mode=list(range(modes)))
    
    # TODO: hits missing mode:
    fieldsTM = mulsol.set(pol=1).field1D(mode=[0, 1, 2])
    
    pnorm = mulsol.field1D_normpower(field=fields)
    mulsol.plotfield1D(field=fields)
    mulsol.plotfield1D(field=fieldsTM)
    mulsol.plotfield1D(field=pnorm)
    
    # Calculate modes in parasitic waveguide
    multi2 = MultiLayerSlabSolverC(view="topview", points=1300)
    multi2.geometry(layers=layers[:3], Nsub = 1.0)
    fields2 = multi2.field1D(wl=1.55, pol=0, mode=list(range(multi2.modes())))
   
    
    # Calculate overlap between both systems.
    # TODO: x coordinated not the same for fld1 and fld2: need spline!
    # OVL_mat = np.zeros((len(fields['data'].keys()),len(fields2['data'].keys())))
    # dx = fields['data'].index.values[1] - fields['data'].index.values[0]
    # for fld1 in fields['data'].columns:
    #     for fld2 in fields2['data'].columns:
    #         OVL_mat[fld1, fld2] = dx * sum([abs(fields['data'][fld1][x] * fields2['data'][fld2][x]) for x in fields2['data'].index])
    #         if OVL_mat[fld1, fld2] > 1e-3:
    #             print(fld1, fld2, 10*np.log10(OVL_mat[fld1, fld2]))
    # print(OVL_mat)