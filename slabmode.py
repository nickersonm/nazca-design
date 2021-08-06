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
#
# @author: Ronald Broeke (c) 2016-2017
# @email: ronald.broeke@brightphotonics.eu
# -----------------------------------------------------------------------

"""Simple slab mode solver."""

from sys import exit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
from .simglobal import sim
from scipy.optimize import fsolve
from . import cfg


class EIM(object):
    """
    Simple 3 layer symmetric slab solver.
    """

    def __init__(
        self, Mr=None, Mb=None, w=2.0, wl=sim.wl, name=None, points=100, view=None
    ):
        if view is not None:
            self.view = view
        else:
            self.view = "topview"

        if name is None:
            self.name = "symmetric-slab-layer-solver"
        else:
            self.name = name

        self._w = w  # waveguide width
        self._wl = wl  # wavelength
        self._pol = 0  # polarization
        self._mode = 0  # modes
        self.points = 100  # number of points in plot of field

        self.Mr = Mr
        self.Mb = Mb

        # self._Nridge = Nridge # ridge index
        # self._Nbg = Nbg    # background index
        self._Nsub = 0.0  # substrate index
        self._Neff = 0.0  # effective index

        self._V = 0.0    # V-parameter
        self._b = 0.0    # propagation constant
        self._dbmin = 1e-12 #needs detail in case of very wide slabs, e.g. Wslab = 1000*lambda
        self.__u = 0.0   # mode lateral wavevector under ridge
        self.__c = 0.0   # mode lateral wavevector outside ridge

        self._maxmodes = 1  # max. number of modes based on V-parameter
        self.__slabFlag = False
        self.__doRecalc = True

        # Parameters used in mode profile:
        self._clip = 0.0
        self.__Ge = 1.0  # factor for even slab modes under the ridge
        self.__Go = self.__Ge
        self.__Pe = 0.0  # even mode constant
        self.__Po = 0.0  # oneven mode constant
        self.__sqrtPe = 0.0  # even mode constant
        self.__sqrtPo = 0.0  # oneven mode constant
        self.__Je = 0.0  # even mode constant
        self.__Jo = 0.0  # oneven mode constant

        self.EFieldResolution = (
            0.05  # in [um] resolution  when plotting a field profile
        )

    def __polRot(self, pol):
        if self.view == "sideview":
            return pol
        elif self.view == "topview":
            if pol == 0:
                return 1
            else:
                return 0
        else:
            return None

    def __getPe(self):
        # TODO: TE/TM?
        # Get power normalization of even mode.
        # Pe = power in mode. Note: E/sqrt(Pe) gives normalized E-field
        # Normalization can be clipped to [X0/2,-X0/2]
        # x0 = abs(x0)
        He = (self.__Ge * m.cos(self.__c * self._w / 2)) ** 2
        Pet = He / self.__u  # tail from -inf to -self._w/2 (edge slab)
        if self._clip / 2 > self._w / 2:
            Pet = Pet - He * exp(self.__u * (self._w - self._clip)) / self.__u  # tail
        Per = (
            self.__Ge
            * self.__Ge
            * (self._w / 2 + m.sin(self.__c * self._w) / (2 * self.__c))
        )  # ridge
        self.__Pe = Per + Pet
        return self.__Pe

    def __getPo(self):
        # Get power normalization of odd mode.
        # Po = power in mode. Note: E/sqrt(Pe) gives normalized E-field
        # Normalization can be clipped to [X0/2,-X0/2]"""
        # x0 = abs(x0);
        Ho = (self.__Go * m.sin(self.__c * self._w / 2)) ** 2
        Pot = Ho / self.__u  # tail from -inf to -d/2 (edge slab)
        if self._clip / 2 > self._w / 2:
            Pot = (
                Pot - Ho * m.exp(self.__u * (self._w - self._clip)) / self.__u
            )  # tail correction when clipping domain at x0
        Por = (
            self.__Go
            * self.__Go
            * (self._w / 2 - m.sin(self.__c * self._w) / (2 * self.__c))
        )  # part under ridge
        self.Po = Por + Pot
        return self.Po

    def __calculateMode(self, **kwargs):
        # if not self.__slabFlag:
        #    print 'bp slabmode: initialize slab properties before setting the mode.'
        #    return None

        self._wl = kwargs.pop('wl', self.wl)
                              # sim.wl)
        self._pol = kwargs.pop('pol', self.pol)
                               # sim.pol)

        # print ('pol=', self._pol)
        # print ('view=', self.view)
        # print ('polrot=', self.__polRot(self._pol))

        self._Nridge = self.Mr.Neff(wl=self._wl, pol=self._pol)
        self._Nbg = self.Mb.Neff(wl=self._wl, pol=self._pol)

        if self._Nridge <= self._Nbg:
            exit("(slabmode): Nr<Nb error.")
        if self._Nbg < 1:
            exit("(slabmode): background index Nb<1 error.")

        self._V = (2 * m.pi * self._w / self._wl) * m.sqrt(
            self._Nridge ** 2 - self._Nbg ** 2
        )
        self._maxmodes = int(self._V / m.pi) + 1
        # one polarization

        # TODO: is maxModes guided modes?
        # printf (mode,"\n");
        if self._mode >= self._maxmodes:
            # print ('bpMessage.Warning ("No guided mode #"+ mode+ " exists (max mode #"+ (self._maxmodes-1)+").\n")')
            self._Neff = 0
            self.__c = 0
            self.__u = 0
            self._b = 0
            self.__Pe = 1
            self.__Po = 1
            self.__sqrtPe = 1
            self.__sqrtPo = 1
            return -1  # mode does not exist.

        b = 0.8 # propagation start value
        dbmin = self._dbmin
        db = (1.0-b)-dbmin
        value = 0
        count = 0  # iteration counter
        countMax = 50  # maximum allowed iterations
        s = 1.0

        if self.__polRot(self._pol) == 0:  # TE, polRot for topview
            while (
                db > dbmin and count < countMax
            ):  # abs to prevent tan from changing sign accidently
                count += 1
                value = self._V * m.sqrt(1 - b) - 2 * m.atan(m.sqrt(b / (1 - b)))
                if value > self._mode * m.pi:
                    db /= 2.0
                    b += db
                    s = 2.0
                else:
                    db /= s
                    b -= db
                    if b < 0:
                        b = 0

        elif self.__polRot(self._pol) == 1:  # TM, polRot for topview
            n2 = self._Nridge ** 2 / self._Nbg ** 2
            while (
                db > dbmin and count < countMax
            ):  # abs to prevent tan from changing sign accidently
                count += 1
                value = self._V * m.sqrt(1 - b) - 2 * m.atan(n2 * m.sqrt(b / (1 - b)))
                if value > (self._mode) * m.pi:
                    db /= 2.0
                    b += db
                    s = 2.0
                else:
                    db /= s
                    b -= db
                    if b < 0:
                        b = 0

        else:
            print(
                'bpMessage.Error ("(slabmode): provided unknown polarization to function InitMode().")'
            )
        if count >= countMax:
            print("(slabmode): Too many iterations: Maybe the ")

        self._b = b
        self._Neff = m.sqrt(
            self._b * (self._Nridge ** 2 - self._Nbg ** 2) + self._Nbg ** 2
        )
        self.__c = (2 * m.pi / self._wl) * m.sqrt(self._Nridge ** 2 - self._Neff ** 2)
        self.__u = (2 * m.pi / self._wl) * m.sqrt(self._Neff ** 2 - self._Nbg ** 2)

        self.__Pe = self.__getPe()
        self.__Po = self.__getPo()
        self.__sqrtPe = m.sqrt(self.__Pe)
        self.__sqrtPo = m.sqrt(self.__Po)

        # continuous at waveguide edge TE and TM:
        self.__Je = self.__Ge * m.cos(self.__c * self._w / 2)
        self.__Jo = self.__Go * m.sin(self.__c * self._w / 2)

        self.__doRecalc = False
        return 0

    @property
    def b(self):
        if self.__doRecalc:
            self.__calculateMode()
        return self._b

    @b.setter
    def b(self, val):
        print("b can not be set")

    @property
    def dbmin(self):
        return self._dbmin

    @dbmin.setter
    def dbmin(self, val):
        print(
            "dbmin can not be set"
        )  # TODO: not yet, first enable readout (18-08-20, SVU)

    @property
    def dbmin(self):
        return self._dbmin

    @dbmin.setter
    def dbmin(self, val):
        print("dbmin can not be set") #TODO: not yet, first enable readout (18-08-20, SVU)

    @property
    def V(self):
        if self.__doRecalc:
            self.__calculateMode()
        return self._V

    @V.setter
    def V(self, val):
        print("V can not be set")

    def __c(self):
        if self.__doRecalc:
            self.__calculateMode()
        return self.__c

    def __u(self):
        if self.__doRecalc:
            self.__calculateMode()
        return self.__u

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, modeNumber):
        self._mode = modeNumber
        self.__doRecalc = True

    @property
    def Nr(self):
        return self.Nr

    @Nr.setter
    def Nr(self, val):
        self._Nr = val
        self.__doRecalc = True

    @property
    def Nbg(self):
        return self._Nbg

    @Nbg.setter
    def Nbg(self, val):
        self._Nbg = val
        self.__doRecalc = True

    @property
    def Nsub(self):
        return self._Nsub

    @Nsub.setter
    def Nsub(self, val):
        self._Nsub = val
        self.__doRecalc = True

    @property
    def clip(self):
        return self.clip

    @clip.setter
    def clip(self, x):
        self._clip = x
        self.__doRecalc = True

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, width):
        self._w = width
        self.__doRecalc = True

    @property
    def wl(self):
        return self._wl

    @wl.setter
    def wl(self, val):
        self._wl = val
        self.__doRecalc = True

    @property
    def pol(self):
        return self._pol

    @pol.setter
    def pol(self, pol):
        self._pol = pol
        self.__doRecalc = True

    def initSlab(
        self,
        Nr,  # Ntv ridge */
        Nbg,  # Ntv background */
        w,  # waveguide width */
        wl,  # wavelength */
        pol=0,  # polarization */
        mode=0,  # mode number
        clip=0.0,
    ):  # where to force the modeprofile to 0 (clipping) */

        """
        Initialize a slab waveguide.

        Args:
             Nr (float): ridge index
             Nbg (float): background index
             w (float): waveguide width in um
             wl (float): wavelength in um
             pol (int): polarization (TE=0, TM=1)
             mode (int): mode number

        Returns:
            int: maximum number of guides modes in the slab
        """
        if Nbg > Nr:
            print("bp: Error: Nbg>=Nr ", Nbg, Nr, "\n")
        self.__slabFlag = True
        self._Nridge = Nr
        self._Nbg = Nbg
        self._w = w
        self._wl = wl
        self._pol = pol
        self._mode = mode
        self._clip = abs(clip)

        # self.setMode(self._mode) # default initialize mode 0
        self.__doRecalc = True
        return self._maxmodes

    # -----------------------------------------------------------------------
    # field profiles
    # -----------------------------------------------------------------------

    def __getEe(self, xx):
        # Return even mode E-field.
        # Clips domain at x0 if x0>d
        X = []
        for x in xx:
            if x < -self._w / 2.0:
                X.append(
                    self.__Je * m.exp(self.__u * (x + self._w / 2.0)) / self.__sqrtPe
                )
            elif x > self._w / 2.0:
                X.append(
                    self.__Je * m.exp(-self.__u * (x - self._w / 2.0)) / self.__sqrtPe
                )
            else:
                X.append(self.__Ge * m.cos(self.__c * x) / self.__sqrtPe)
        return np.array(X)

    def __getEo(self, xx):
        # Return odd mode E-field.
        # Clips domain at x0 if x0>d/2
        X = []
        for x in xx:
            if x < -self._w / 2.0:
                X.append(
                    -self.__Jo * m.exp(self.__u * (x + self._w / 2.0)) / self.__sqrtPo
                )
            elif x > self._w / 2.0:
                X.append(
                    self.__Jo * m.exp(-self.__u * (x - self._w / 2.0)) / self.__sqrtPo
                )
            else:
                X.append(self.__Go * m.sin(self.__c * x) / self.__sqrtPo)
        return np.array(X)

    def E(self, x, wl=sim.wl, pol=sim.pol):
        """
        Calculate Efield at position <x> for wavelength <wl> and polarization <pol>.

        Args:
            x (float): spatial postional along
            wl (float): optional, wavelength in um
            pol (int): optional, polarization (TE=0, TM=1)

        Returns:
            float: E field strength at <x>,  normalized per um
        """
        if self._wl != wl:
            self._wl = wl
            self.__doRecalc = 1
        if self._pol != pol:
            self._pol = pol
            self.__doRecalc = 1
        if self.__doRecalc:
            self.__calculateMode()
        if self._clip != 0 and abs(x) > self._clip:
            return 0
        absx01 = abs(self._clip)
        if self._mode % 2 == 0:
            Ef = self.__getEe(x)
        else:
            Ef = self.__getEo(x)
        return Ef

    def getEwidth(self, Efrac=0.5):
        """
        Find width where field has fraction Efrac left w.r.t. center of guide.

        Args:
            Efrac (float): fraction between 0 and 1 (defailt = 0.5, HWHM)

        Returns:
            float: distance of mode center where filds strength drops to Efrac
        """
        if self.__doRecalc:
            self.__calculateMode()
        E0 = self.getE(0)
        Elimit = Efrac * E0
        x = 0.0
        dx = 0.3
        dxmin = 0.001
        t = 0
        tmax = 100
        E = 0.0
        while t < tmax and dxmin < dx:
            E = self.getE(x)
            if E > Elimit:
                x += dx
            else:
                x -= dx
                x = abs(x)
                dx /= 2.0

            t += 1
            # printf(t, " ", x, " ", E, "\n");

        if t == tmax:
            print('bpMessage.Warning ("Too many loops in getEwidth()\n")')
            return 0
        else:
            return x

    def Neff(self, **kwargs):
        self._wl = kwargs.pop("wl", sim.wl)
        self._pol = kwargs.pop("pol", sim.pol)
        if self._wl is None:
            print("in Neff call, no wl parameter found.")
            self._Neff = 0
        else:
            if self.__doRecalc:
                self.__calculateMode(**kwargs)
        return self._Neff

    # def Neff_bak(self, mode=-1):
    #    oldMode = self._mode
    #    if mode==-1:
    #        mode = self._mode
    #    else:
    #        self._mode = mode
    #    if self.__doRecalc:
    #        self.__calculateMode()
    #    self._mode = oldMode
    #    return self.Neff

    @property
    def maxModes(self):
        if self.__doRecalc:
            self.__calculateMode()
        return self._maxmodes  # // calculated in slab init

    @maxModes.setter
    def maxModes(self, val):
        self._maxmodes = val

    #    #print "Can not set maxModes. skipping."

    @property
    def numberOfGuidedModes(self):
        oldMode = self._mode
        self.__calculateMode()
        gm = 0
        mm = self._maxmodes
        # print ('_maxmodes=', mm)
        for m in range(mm):
            #    #if self.__doRecalc:
            self._mode = m
            self.__calculateMode()
            # print ('_mode=', self._mode)
            # print ('_Neff=', self._Neff)
            # print ('_Nsub=', self._Nsub)
            if self._Neff > self._Nsub:
                gm += 1
            else:
                break
        self._mode = oldMode
        self.__calculateMode()
        # print ('gm-1=', gm-1)
        return gm

    @numberOfGuidedModes.setter
    def numberOfGuidedModes(self):
        print("Can not set numberOfGuidedModes. skipping.")

    def printModeInfo(self, mode=-1):
        if mode == -1:
            mode = self._mode
        else:
            self.setMode(mode_)

        oldMode = self._mode
        if self.__doRecalc:
            self.__calculateMode()

        guided = "no"
        if self._Neff > self._Nsub:
            guided = "yes"

        print(
            "Slab-modeinfo:",
            "\n- mode#         = ",
            self._mode,
            "\n- Neff          = ",
            self._Neff,
            "\n- w             = ",
            self._w,
            "\n- wl            = ",
            self._wl,
            "\n- pol           = ",
            self._pol,
            "\n- guided-modes  = ",
            self.numberOfGuidedModes,
            "\n- Nr            = ",
            self._Nridge,
            "\n- Nbg           = ",
            self._Nbg,
            "\n- Nsub          = ",
            self._Nsub,
            "\n- beta          = ",
            self.beta,
            "\n- V             = ",
            self._V,
            "\n- u             = ",
            self.__u,
            "\n- c             = ",
            self.__c,
        )
        self._mode = oldMode

    def printAllModes(self):
        print("mode# Neff guided")
        for m in range(self._maxmodes):
            print(
                m, " ", self.getNeff(m),
            )
            if self._Neff > self._Nsub:
                print(" yes")
            else:
                print(" no")

    def __str__(self):
        out = "Slab-solver settings:"
        out += "\n- slab name                            :'" + self.name + "'"
        out += (
            "\n- type                                 : S3S (symmetric-3-layer-solver)"
        )
        out += "\n- view                                 : " + self.view
        out += "\n- wavelength                 - wl [um] : " + str(self._wl)
        out += "\n- polarization               - pol     : " + str(self._pol)
        out += "\n- materials                  - M       : " + ", ".join(
            [self.Mb.name, self.Mr.name, self.Mb.name]
        )
        out += "\n- layer thickeness           - D [um]  : " + str(self.w)
        out += "\n- refractive indices         - N       : " + ", ".join(
            map(str, [self.Mb.Neff(), self.Mr.Neff(), self.Mb.Neff()])
        )
        out += "\n- number of points for Field - points  : " + str(self.points)
        out += "\n"
        return out

    def getField(self):
        x = np.linspace(-1.5 * self.w, 1.5 * self.w, self.points)
        E = self.E(x)
        field = pd.DataFrame({"x": x, "E": E})
        return field

    def Efield(self):
        output = pd.DataFrame({})
        output["x"] = np.linspace(-1.5 * self.w, 1.5 * self.w, self.points)
        output["E"] = self.getField()
        output["E"] = [m ** 2 for m in output["E"]]
        return output

    def plot(self, title=None, ax=None, **kwargs):
        F = self.getField()
        if title is None:
            title = self.name

        field = self.getField()
        x = kwargs.pop("x", "x")
        y = kwargs.pop("y", "E")
        clip_on = kwargs.get("clip_on", False)
        xlabel = "x [um]"
        ylabel = "E [power normalized]"

        if self.pol == 0:
            col = "b"
        else:
            col = "g"

        if ax is None:
            fig, ax = plt.subplots(figsize=cfg.modeplotsize)

        [ax.spines[key].set_linewidth(2.0) for key in ax.spines.keys()]
        # if self.plotNum >= 0:
        ax.set_title(title)

        F.plot(x=x, y=y, color=col, clip_on=clip_on, ax=ax, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-1.5 * self.w, 1.5 * self.w)
        plt.tight_layout()

        # self.plotNum += 1
        return ax
