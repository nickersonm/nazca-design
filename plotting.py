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
from .simglobal import sim
import matplotlib.pyplot as plt


def __plotProperty1D(F, R, w, wl, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
        plotHere = False
    else:
        plotHere = True

    linewidth = 3
    if len(R) > 1:
        Y = [F(R=Ri, width=w[0], wl=wl[0]) for Ri in R]
        ax.plot(R, Y, linewidth=linewidth)
        ax.set_xlabel('radius [um]')

    if len(w) > 1:
        Y = [F(R=R[0], width=wi, wl=wl[0]) for wi in w]
        ax.plot(w, Y, linewidth=linewidth)
        ax.set_xlabel('width [um]')

    if len(wl) > 1:
        Y = [F(R=R[0], width=w[0], wl=wli) for wli in wl]
        ax.plot(wl, Y, linewidth=linewidth)
        ax.set_xlabel('wavelength [um]')

    # plt.ylabel('offset [um]')
    ax.grid()
    if plotHere:
        plt.show()
    return None


def __plotProperty2D(F, R, w, wl, ax=None):
    print(ax)
    if ax is None:
        fig, ax = plt.subplots()
        plotHere = False
    else:
        plotHere = True

    if len(R) == 1:
        Y = [[F(R=R[0], width=wi, wl=wli) for wi in w] for wli in wl]
        img = ax.contourf(w, wl, Y)
        ax.set_xlabel('width [um]')
        ax.set_ylabel('wavelength [um]')
    if len(w) == 1:
        Y = [[F(R=Ri, width=w[0], wl=wli) for Ri in R] for wli in wl]
        img = ax.contourf(R, wl, Y)
        ax.set_xlabel('radius[um]')
        ax.set_ylabel('wavelength [um]')
    if len(wl) == 1:
        Y = [[F(R=Ri, width=wi, wl=wl[0]) for Ri in R] for wi in w]
        img = ax.contourf(R, w, Y)
        ax.set_xlabel('radius [um]')
        ax.set_ylabel('width [um]')

    plt.colorbar(img, ax=ax)
    ax.grid()
    if plotHere:
        plt.show()
    return None


def plotProperty(F, ax=None, **kwargs):
    R = kwargs.pop('R', sim.R)
    w = kwargs.pop('width', 0)
    wl = kwargs.pop('wl', sim.wl)

    try: len(R)
    except TypeError: R = [R]
    try: len(w)
    except TypeError: w = [w]
    try: len(wl)
    except TypeError: wl = [wl]

    t = 0
    if len(R) > 1: t += 1
    if len(w) > 1: t += 1
    if len(wl) > 1: t += 1
    if t < 1 or t > 2:
        print('Please provide 1 or 2 list arguments, not', t)
        return 0
    elif t == 1:
        __plotProperty1D(F, R, w, wl, ax=ax)
    elif t == 2:
        __plotProperty2D(F, R, w, wl, ax=ax)
    return None
