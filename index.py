#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generic index model base class and utilities.

Class IndexBase() is to be used as a parent class that defines all common 
index model functions.  Many of them are fallback functions that should be 
defined in the derived child index model, i.e. the child should
defines methods as Nridge, Nsub, end Neff, be it as constants, 
functions (like fits) or 1D or 2D solvers. The signatures of the respective
methods shall match 

Example:

    If a new index model is defined, include IndexBase as the argument:

        from index import IndexBase
        class NewIndexModel(IndexBase):
            ...

It also contains a generic plot function to do a visual check for the data (plot_index).
Furthermore, the Separator function can isolate a set of values from a list of tuples.

Bright Photonics B.V. (c) 2020-2021
"""

import numpy as np
import matplotlib.pyplot as plt

import nazca as nd
from .version import __version__

nd.logger.info("Index models version: {}".format(__version__))
from nazca.simglobal import sim


class IndexBase:
    """Base class to unify implementations of index models.

    The class defines two types of methods:
    1. Placeholder/fallback index methods that should be defined in a derived index class.
    2. Model independent methods.

    Anticipated index model types to build upon IndexBase::

        'float': constant number (independent of input parameters).
        'slab' : 1D slab-solver / 1D multi-layer solver.
        '2Dslab' vertical 1D index slab-solver + horizontal 1D index slab-solver
        'poly' : vertical index as polynomial as function of N(wl, pol).
        'data' : Interpolated value generated from a data set.
        'diff' : Diffusion model N(wl, pol, wid
    """

    def __init__(self):
        """Initialize an IndexBase object.

        Returns:
            None
        """
        # self.dwl = 2e-2  # step size for group index calculation in um. (Reproduce results)
        self.dwl = 1e-3  # step size for group index calculation in um. (Best results)
        self.flagNridge = False  # store if missing Nridge() as been reported.
        self.flagNbg = False  # store if missing Nbg() as been reported.
        self.flagNsub = False  # store if missing Nsub() as been reported.
        self.flagNeff = False  # store if missing Neff() as been reported.
        self.flagdNradius = False  # store if missing dNraidus() as been reported.
        self.flagMaxMode = False  # store if missing maxMode function has been reported.
        self.flagModesetter = False  # store if missing modesetter() has been reported.
        self.flagsboffset = False  # store if missing sboffset() has been reported.

        self.properties = {
            "name": "Index Base",
            "modelNr": None,
            "modelNbg": None,
            "modelNeff": None,
            "Neff2D": None,
        }

        self.mode = 0
        self.max_mode = None

    def __str__(self):
        """Index model information.

        Retrieves and sums up all properties of the index model.
        """
        out = "Index model properties:\n"
        long_key = max([len(key) for key in self.properties.keys()])
        for n, key in enumerate(self.properties.keys()):
            out += key + (long_key + 1 - len(key)) * " " + f": {self.properties[key]}"
            if n != len(self.properties.keys()) - 1:
                out += "\n"
        return out

    def setMode(self, mode):
        """Set the desired mode to study.

        NOTE: do not directly overwrite this function. This will remove the
        check in place.

        First checks if the mode is available for setting.
        If so, execute the 'modesetter()' function, which can be overwritten
        locally.

        Args:
            mode (int): Waveguide mode number.

        Returns:
            None
        """
        if mode != self.mode:
            if self.max_mode is None:
                self.getMaxMode()
            if mode <= self.max_mode and mode >= 0:
                self.mode = mode
                self.modesetter(mode=mode)
            else:
                nd.main_logger(
                    f"Mode could not be set to {mode} (Maximum mode: {self.max_mode}). "
                    f"Using the old mode setting (mode = {self.mode}).",
                    "error",
                )

    def validateModes(self, mode):
        """Validate which modes in <modes> are existing modes.

        Filter out the modes below the maximum mode.
        All deviating entries will be printed to the log file.

        Args:
            modes (list): List of modes to be validated.
                If 'all' is given, this will calculate all available modes.

        Raises:
            Exception:
                If none of the values in 'modes' can be calculated.

        Returns:
            list: filtered out guided modes in <modes>
        """
        if self.max_mode is None:
            self.getMaxMode()
        if isinstance(mode, list):
            if mode == []:
                nd.main_logger(
                    "No modes were given to 'multimode_Neff'. Returning Neff for mode = 0.",
                    "error",
                )
                mode = [0]
            elif min(mode) > self.max_mode:
                raise Exception(
                    f"Neff can not be determined for modes {mode}\n"
                    f"Maximum mode: {self.max_mode}"
                )
            elif max(mode) > self.max_mode:
                rejects = []
                new_mode = []
                for m in mode:
                    if m > self.max_mode:
                        rejects.append(m)
                    else:
                        new_mode.append(m)
                mode = new_mode
                nd.main_logger(
                    f"{self.properties['name']}: Neff can not be calculated for modes {rejects}"
                    f" (maximum mode = {self.max_mode}). Continuing with mode(s) {new_mode}.",
                    "warning",
                )
        elif mode == "all":
            mode = [m for m in range(self.max_mode + 1)]
        else:
            raise Exception("Input for mode={modes} not recognized.")
        return mode

    def Ngrp(self, wl=None, pol=None, width=None, radius=0, mode=None, **kwargs):
        """Calculate the group index.

        Calculate the group index of a straight or bend waveguide based on Neff().

        Args:
            width (float): width of the ridge in um.
            wl (float): wavelength of the light in um.
            pol (int): polarization of the light. Use 0 for TE and 1 for TM.
            radius (float): radius of the bend at the center of the waveguide
                in um. The default is 0.
            modes (list): list of modes to calculate Ngrp for. Default is None.

        Returns:
            float | dict: Ngrp as float for modes is None, else dict {modenumber: Ngrp}
        """
        if pol is None:
            pol = sim.pol
        if wl is None:
            wl = sim.wl

        Neff = self.Neff(wl=wl, pol=pol, width=width, radius=radius, mode=mode, **kwargs)
        Neff2 = self.Neff(wl=wl + 0.5 * self.dwl, pol=pol, radius=radius, width=width, mode=mode, **kwargs)
        Neff1 = self.Neff(wl=wl - 0.5 * self.dwl, pol=pol, radius=radius, width=width, mode=mode, **kwargs)
       
        if mode is not None:
            Ngrp = {}
            for m in Neff.keys():
                Ngrp[m] = Neff[m] - wl * ((Neff2[m] - Neff1[m]) / self.dwl)
            return Ngrp
        else:
            return Neff - wl * ((Neff2 - Neff1) / self.dwl)

# =============================================================================
# PLACEHOLDER FUNCTIONS
# =============================================================================
    def getMaxMode(self):  # rename function
        """Placeholder to obtain the maximum guided mode number.

        Args:
            None

        Returns:
            int: Maximum allowed mode number.
        """
        if self.max_mode is None:
            nd.main_logger(
                f"Function 'getMaxMode' needs to be overwritten for {self.properties['name']}. "
                "Setting 'max_mode' to 0.",
                "error",
            )
            self.max_mode = 0
        return self.max_mode


    def modesetter(self, mode):  # rename function
        """Placeholder function for local mode script.

        Should be overwritten in the index model.

        Args:
            mode (int): Waveguide mode number.

        Returns:
            None
        """
        if not self.flagModesetter:
            nd.main_logger(
                f"Function 'modesetter' is not defined for {self.properties['name']}. "
                "Mode is not transferred to underlying functions.",
                "warning",
            )
            self.flagModesetter = True


    def Neff(self, wl, pol, width, radius=0.0, mode=0):
        """Placeholder Neff function.

        Obtain the effective index of the waveguide.

        Writes a log error if Neff() is not defined in a derived class.
        To be defined in the derived index class.

        Args:
            width (float): width of the ridge in um.
            wl (float): wavelength of the light in um.
            pol (int): polarization of the light. Use 0 for TE and 1 for TM.
            radius (float): Radius at the center of the waveguide in um.
                Default is 0.0 (straight waveguide).
            modes (list): list of modes to calculate Neff for. Default is None.

        Raises:
            Exception:
                If this function is not overwritten in the index model.

        Returns:
            None
        """

        if not self.flagNeff:
            nd.main_logger(
                "Method Neff(self, width, wl, pol) needs to be defined for ({self.properties['name']}). Returning None.",
                "error",
            )
        self.flagNeff = True
        return None


    def dNradius(self, wl, pol, mode=0, width=None, radius=None):
        """Placeholder dNradius function.

        Obtain the wavegudie effective index change due to curvature
        w.r.t. to a straight waveguide.
        To be defined in the derived index class.

        Writes a log error if dNradius() is not defined in a derived class.

        Args:
            width (float): width of the ridge in um.
            wl (float): wavelength of the light in um.
            pol (int): polarization of the light. Use 0 for TE and 1 for TM.
            radius (float): radius of the bend at the center of the waveguide
                in um.

        Returns:
            float: correction on Neff for bend curvature index (default = 0.0).
        """
        # TODO: check for higher order modes separately
        if not self.flagdNradius:
            nd.main_logger(
                f"No dNradius defined for index model '{self.properties['name']}', mode={mode}. Returning 0.0.",
                "error",
            )
        self.flagdNradius = True
        return 0.0


    def Nridge(self, wl, pol):
        """Placeholder Nridge function.

        Obtain the ridge effective index.
        To be defined in the derived index class.

        Writes a log error if Nridge() is not defined in a derived class.

        Args:
            wl (float): wavelength of the light in um.
            pol (int): polarization of the light. Use 0 for TE and 1 for TM.

        Returns:
            None
        """
        if not self.flagNridge:
            nd.main_logger(
                "Method Nridge(self, wl, pol) needs to be defined for {self.properties['name']}. Returning None.",
                "error",
            )
        self.flagNridge = True
        return None


    def Nbg(self, wl, pol):
        """Placeholder Nbg function.

        Obtain the background effective index.
        To be defined in the derived index class.

        Writes a log error if Nbg() is not defined in a derived class.

        Args:
            wl (float): wavelength of the light in um.
            pol (int): polarization of the light. Use 0 for TE and 1 for TM.

        Returns:
            None
        """
        if not self.flagNbg:
            nd.main_logger(
                "Method Nbg(self, wl, pol) needs to be defined for {self.properties['name']}. Returning None.",
                "error",
            )
        self.flagNbg = True
        return None


    def Nsub(self, wl, pol):
        """Placeholder Nsub function.

        Obtain the substrate index. Note this is the substrate the optical
        mode experiences. This could be different from the physical/mechanical
        substrate below it.

        To be defined in the derived index class.

        Writes a log error if Nsub() is not defined in a derived class.

        Args:
            wl (float): wavelength of the light in um.
            pol (int): polarization of the light. Use 0 for TE and 1 for TM.

        Returns:
            None
        """
        if not self.flagNsub:
            nd.main_logger(
                "Method Nsub(self, wl, pol) needs to be defined for {self.properties['name']}. Returning None.",
                "error",
            )
        return None


    def sboffset(self, width, radius):
        """Placeholder for the offset from straight to bend waveguides.

        To be defined in the derived index class.

        Args:
            wl (float): wavelength
            pol (int): polarisation

        Returns:
            float: straight-to-bend offset value
        """
        if not self.flagsboffset:
            nd.main_logger(
                f"Function 'modesetter' is not defined for {self.properties['name']}. "
                "Mode is not transferred to underlying functions.",
                "warning",
            )
            self.flagsboffset = True
        return 0.0

    @property
    def properties(self):
        self._propDict["max mode"] = self.max_mode
        self._propDict["mode"] = self.mode
        return self._propDict

    @properties.setter
    def properties(self, newDict):
        self._propDict = newDict


# =============================================================================
# Auxiliary functions
# =============================================================================
def legacyspline2nazca(splinedata):
    """Reorganize legacy spline data [wl, Neff, ...] into [wl], [Neff].

    legacy: [wl1, Neff1, wl2, Neff2, ...]
    out: [wl1, wl2, ...], [Neff1, Neff2, ...]

    Args:
        data (list): spline data in legacy format.

    Returns:
        list, list: wl, Neff
    """
    wls = []
    Neffs = []
    for i, value in enumerate(splinedata):
        if i % 2 == 0:
            wls.append(value)
        else:
            Neffs.append(value)
    return wls, Neffs


def plot_index(name, dat, xdat="Lambda", ydat="Refractive index"):
    """Visualize index models in a generic plot.

    Can be used as a quick visual check to see if the behavior of the index
    model is stable.

    Args:
        name (str): name of the library to be used in the plot title.
        dat (pd.DataFrame): Pandas DataFrame with the data to be plotted.
        xdat (str): key for the x-data. Default is 'Lambda'.
        ydat (str): key for the y-data. Default is 'Refractive index'.

    Returns:
        None
    """
    ylabels = [m for m in dat.columns if m != xdat]
    dat.plot(x=xdat, y=ylabels)
    plt.title(f"Neff and Ngrp for {name}")
    plt.xlabel(f"{xdat} (um)")
    plt.ylabel(f"{ydat}")
    plt.legend()
    plt.show()


def _input_check(order, Np):
    """Function to validate if the input parameters are valid.

    Ensures multi-parameter polynomial does not break down for the given order
    and number of input parameters "Np".
    If input values are invalid, returns proper information to adjust.

    Args:
        order (int): desired polynomial order.
        Np (int): Number of input parameters.

    Raises:
        Exception:
            If 'order' < 1 (no polynomial)
        Exception:
            If 'order' > 3 (possibility of over-fitting)
        Excetpion:
            If Np > 6 (possibility of over-fitting)

    Returns:
        bool: 'True' if order is within function limits.
    """
    min_order = 1
    max_order = 4
    max_Np = 6

    advice_order = f"Please use an order between {min_order} and {max_order}."
    advice_Np = f"Please use {max_Np} input parameters or less."
    if type(order) != int:
        raise Exception(f"Order should be an integer (current type: {type(order)}).")
    if order < 1:
        raise Exception(f"Order ({order}) is less than 1. " + advice_order)
        return False
    elif order > 4:
        raise Exception(
            f"Order ({order}) not available for multi-dimensional fit. " + advice_order
        )
        return False
    if Np > max_Np:
        raise Exception(
            f"Number of input parameters ({Np}) is too large for multi-dimensional fit. "
            + advice_Np
        )
    else:
        return True


def _generate_xmatrix(order, xvalues):
    """Function to generate the x-coordinates for the polynomial.

    This function generates an array of x-parameters from 'xvalues' to match the
    generated polynomial coefficients from the _generate_polycoef.

    Args:
        order (int): the order of the generated polynomial.
        xvalues (list): list of float input values for a single polynomial point.

    Raises:
        Exception:
            If 'order' is less than 1 (not possible to generate x-parameters).

    Returns:
        list: list of polynomial x-coordinates.
    """
    if _input_check(order=order, Np=len(xvalues)):
        out = [1,]  # Output list for x-values
        xvalues = list(xvalues)
        transmat = [[i] for i in xvalues]  # (Initial) transformation matrix.
        for m in transmat:  # Order = 1
            out += m

        for d in range(2, order + 1):
            newmat = []  # Build transformation matrix for next order.
            for n in range(len(xvalues)):
                mult = []  # Entry n of new transformation matrix.
                for k in range(n, len(transmat)):  # To prevent repetition of multiplication products
                    mult += [xvalues[n] * m for m in transmat[k]]
                out += mult  # Add multiplication products to output list
                newmat.append(mult)
            transmat = newmat  # Implement transformation matrix for next order.
    return out


def MultiDpoly(order, xi, coef, **kwargs):
    """Calculate output of a multi-dimensional polynomial.

    This function generates the output of a multi-dimensional polynomial function
    based on the given order, coefficients and input parameters.

    For the example having order = 2 and xi = [(x1a, x2a), (x1b, x2b), ...] 
    the polynomial coefficients should be in the form 
   
    coef = [a00, a01, a02, ... a05] 
    
    creating an output value following
    
    y = a00 + a01 * x1 + a02 * x2 + a03 * x1**2 + a04 * x1 * x2 + a05 x2 **2
    
    The coefficients can be calculated with e.g. sklearn.preprocessing PolynomialFeatures
    
    Args:
        order (int): polynomial degree.
        xi (list of tuples): input parameters for the functions.
            Defined as [(x1, y1, z1), (x2, y2, z2), ...].
        coef (numpy.array): polynomial coefficients of order 'order'.

    Raises:
        Exception:
            If the polynomial x-values and coefficients do not match in length.
            This indicates a mismatch in 'order' or xi values for the coefficients.
            Check the given coefficients or feed the raw data into this function.

    Returns:
        list: output values for the input points xi.
    """
    if type(xi[0]) in [tuple, list, np.ndarray]:
        try:
            xdummy = _generate_xmatrix(order=order, xvalues=list(xi[0]))
            xshape1 = len(xdummy)
            assert xshape1 == len(coef)
        except AssertionError:
            raise Exception(
                "Order or number of parameters do not match for coefficients and input values\n"
                + "Check the order or recalculate the polynomial coefficients."
            )
    elif _input_check(order, len(xi)):
        xshape1 = order

    xinput = np.zeros((len(xi), xshape1))
    for n in range(len(xi)):
        xinput[n, :] = np.array(_generate_xmatrix(order=order, xvalues=xi[n]))
    ymatrix = xinput * coef
    yvals = list(np.sum(ymatrix, axis=1))
    if len(yvals) == 1:
        yvals = yvals[0]
    return yvals
