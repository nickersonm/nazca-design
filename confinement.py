#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Draft dode by Stef van Uden.

Was part of the mutlilayer solver but need to be generalizzed first most likely,
or should overrule SolverBase function.

Ronald Broeke
"""


    # TODO: RB: Confine1D should be in a field processing class, not in the specific solver.
    #       Move to another module alltogether.
    def confinement1D(self, **kwargs):
        """
        TODO edit docstring
        Calculate the confinement of the ridge, based on (optional) custom parameters

        Parameters
        ----------
        **kwargs :
            wl (float): wavelength in um
            layers (array): array with stacks/materials.
            widths (array): widths of the stacks/materials.
            xmin (float): minimum x-value for the plot.
            xmax (float): maximum x-value for the plot.
            propOn (bool): enable the properties dictionary from the confinement calculation.
            plotOn (bool): enable the plot to visualize the confinement.
            edgeOn (bool): enable the confinement region in the plot.

        Returns
        -------
        IntEsq : Normalized confinement in the ridge layer.
                (Including properties if propOn is True)

        """
        datname = "1D Confinement"  # Type of data (for dictionary)

        self.wl = kwargs.pop("wl", self.wl)  # Wavelength in um
        layers = kwargs.pop("layers", self.Mat)
        self.Mat = layers
        widths = kwargs.pop(
            "widths", self.D
        )  # Obtain the width of the ridge from the stack
        self.D = widths
        Nself = self.Neff()
        width = widths[1]
        xmin = kwargs.pop("xmin", -widths[1] / 2 - widths[0])  # Minimum x in the plot
        xmax = kwargs.pop("xmax", widths[1] / 2 + widths[2])  # Maximum x in the plot
        xpoints = 101  # fixed by
        # xpoints = kwargs.pop('points',101)                  # Number of points to plot
        propOn = kwargs.pop("propOn", False)  # Value to enable properties-dict
        plotOn = kwargs.pop("plotOn", False)  # Value to enable/disable plots
        edgeOn = kwargs.pop("edgeOn", True)  # Value to enable/disable layer boundaries

        dx = (xmax - xmin) / (xpoints - 1)
        ind_low = int((-width / 2 - xmin) / dx)  # Lower index of the background layer
        ind_high = int((width / 2 - xmin) / dx)  # Higher index of the background layer

        # Generate a dictionary with the x-coordinates and E-field
        E_data = self.field1D
        x = E_data["x"]  # Obtain the x-coordinates from the dictionary
        E = E_data["E"]  # Obtain the E-field from the dictionary
        Esq = [m ** 2 for m in E]  # Calculate E-squared
        IntEsq = sum(
            [m * dx for m in Esq[ind_low:ind_high]]
        )  # Calculate the confinement in the Ridge layer
        IntEsq = IntEsq / sum([m * dx for m in Esq])  # Normalization to full range
        self.confinement = IntEsq

        if self.pol == 0:
            polstr = "TE"
        elif self.pol == 1:
            polstr = "TM"
        else:
            print(
                f"WARNING: polarization (pol={self.pol}) is not used for TE or TM. Continue at your own risk."
            )
            polstr = None

        if plotOn == True:  # If plots are enabled:
            fig, ax = plt.subplots()
            if edgeOn == True:  # Plot layer boundaries if enabled
                ax.vlines(
                    x=[self.D[0], self.D[0] + width], ymin=min(Esq), ymax=max(Esq)
                )
            ax.plot(x, Esq)  # Plot the E-squared against the x-coordinates
            ax.set_xlabel("x (um)")
            ax.set_ylabel("Normalized E^2")
            ax.set_title(
                f"Confinement in ridge: {round(IntEsq,5)}\n mode: {self.mode}, pol: {polstr}\n{self.name}"
            )
            plt.show()

        if propOn == True:
            IntEsq = {}
            IntEsq["Confinement_ridge"] = self.confinement
            IntEsq["Confinement_bg(sub+clad)"] = 1 - self.confinement
            IntEsq["properties"] = {
                "solver": self.name,
                "solver_accuracy(b)": self.accuracy,
                "data": datname,
                "keys": [m for m in IntEsq.keys()],
                "wl": self.wl,
                "pol": self.pol,
                "mode": self.mode,
                "materials": [m.name for m in self.Mat],
                "neff_(mat)": [m.Neff(wl=self.wl) for m in self.Mat],
                "width_(mat)": self.D,
                "settings": {
                    "xmin": min(E_data["x"]),
                    "xmax": max(E_data["x"]),
                    "points": self.accuracy,
                },
                "Neff": Nself,
                "Pmax": max(Esq),  # Maximum E-squared (E^2)
                "Ptot": sum(Esq)
                * (max(E_data["x"]) - min(E_data["x"]))
                / xpoints,  # Integral E^2 dx from xmin to xmax
            }
            for i in IntEsq["properties"].keys():
                print(f"{i}: {IntEsq['properties'][i]}")

        return IntEsq