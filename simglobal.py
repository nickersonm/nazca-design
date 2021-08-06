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
#sim need to be a singleton
#"better" to use not a class but put it in a module.


class Sim():
    def __init__(self, wl=1.55, pol=0, modes=(0, 1)):
        """"Initilialize golbl simulation settings.

        Args:
            wl (float); wavelength
            pol (int): polarization
            modes (tuple): list of modes to built a netlist for. Set to () for none.

        Returns:
            None
        """
        self.wl = wl
        self.pol = pol
        self.modes = modes
        #self.radius = radius

    def __str__(self):
        s = ("Global simulation settings:\n")
        s += (f"    wl    = {self.wl:3f}\n")
        s += (f"    pol   = {self.pol}\n")
        s += (f"    modes = {self.modes}\n")
        return s


sim = Sim()
#sim.wl = 1.55
#sim.pol = 0
#sim.modes = (0, 1)