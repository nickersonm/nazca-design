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

# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#
# =============================================================================
# (c) 2018 Katarzyna Lawniczuk, Ronald Broeke
# =============================================================================

"""
DIE template initialization
"""

from itertools import count
import nazca as nd
import nazca.bb_util as bbu
import nazca.geometries as geom


#==============================================================================
# Die template
#==============================================================================
class DesignArea():
    """Class containing a cell template for Die, Design Area and cleaving/dicing."""
    instnum = count()

    def __init__(self, name='DesignArea',
            length=4600, width=4000,
            cleave=100,
            coatingE='NO', coatingW='NO', coatingN=None, coatingS=None,
            die_org=(0, 0),
            pitch=12.5):
        """Construct a DesignArea object.

        The following terminology applies:
        -----------------
        "die": The physical chip after cleaving or dicing.

        "design area": Area in which to draw, which may be larger than the die area.
        The design area incorporates a cleaving/dice region.
        Waveguides may have to extend across the cleave/dice line.

        "cleave/dice line": Target location for cleaving or dicing.

        "cleaving/dicing region": Region around the cleave/dice-line
            indication the cleave or dice tolerance
        -----------------

        Each foundry tends to have its own definition for the gds origin in a
        a MPW unit cell, e.g. the center of the cell, the die bottom-left,
        or the bottom-left point outside the cleave-area.
        The 'die_org' keyword allows to position the die properly on the gds.

        Facet coatings can be indicated for each side of the die.
        A None value ignores any coating.


        Args:
            name (str): name of Design Area cell.
            length (float): (physical) die length in um in x-direction
            width (float): (physical) die width in um in y-direction
            die_org ((float, float)): (x, y) gds position of the die origon,
                i.e. the die's bottom-left corner.
                Note different foundries apply different origins.
            cleave (float): width of cleaving area
            coatingE (str): type of coating on east facet, options: AR | HR | NO | DC
            coatingW (str): type of coating on west facet, options: AR | HR | NO | DC
            coatingN (str): type of coating on north facet, options: AR | HR | NO | DC
            coatingS (str): type of coating on south facet, options: AR | HR | NO | DC
            pitch (float): pitch of optical IO at facet (default = 12.5)
                No IO are added for pitch = None

        Returns:
            None
        """
        self.name   = name
        self.length = length
        self.width  = width
        self.cleave = cleave
        self.coatingE = coatingE
        self.coatingW = coatingW
        self.coatingN = coatingN
        self.coatingS = coatingS
        self.die_org = die_org
        self.pitch = pitch
        self.layer_die = 2
        self.layer_die_user_area = 2
        self.layer_cleave_area = 2
        self.layer_cleave_line = 2


        self.coating_layers = {
            'AR': 'Coating_AR',
            'HR': 'Coating_HR',
            'NO': 'Coating_NO',
            'DC': 'Coating_DC',
        }

        self.coating_width = {
            'AR': 100,
            'HR': 100,
            'NO': 0.5*cleave,
            'DC': 100,
        }

        self.arrow = bbu.make_pincell(layer='package_pin')


    def size(self, length=4600, width=4000, cleave=100):
        """Set size of the (physical) die and the width of the cleaving area.

        Note that the 'die' is the physical size after cleaving.
        The 'design area' also includes the cleaving area and can be larger
        than the die.

        Returns:
            None
        """
        self.length = length
        self.width = width
        self.cleave = cleave


    def coating(self, coatingE='NO', coatingW='NO', coatingN=None, coatingS=None):
        """Define the facet coatings on the die.

        Available coating options: AR | HR | NO | DC
            AR: Anti-reflection
            HR: High-reflection (not available in standard MPW)
            NO: No coating (facets as-cleaved, default)
            DC: Don't care (specify if your circuit has no waveguides on that side

        Note that coating options must be aligned/verified with the technology used.

        Args:
            coatingE (str): East
            coatingW (str): West
            coatingN (str): North
            coatingS (str): South

        Returns None
        """
        self.coatingE = coatingE
        self.coatingW = coatingW
        self.coatingN = coatingN
        self.coatingS = coatingS


    def cell(self, name=None):
        """Create a foundry cell with cleave/dice line and area, and coating info.

        Returns:
            Cell: design area
        """
        if name is None:
            name = '{}_{}'.format(self.name, next(self.instnum))
        with nd.Cell(name=name) as C:
            hcleave = 0.5*self.cleave

            org_die = nd.Pin('org_die', show=False).put(self.die_org[0], self.die_org[1])
            nd.Pin('bbox', show=False).put(org_die.move(0, 0.5*self.width))
            bbu.put_boundingbox('bbox', self.length, self.width)

            #USER AREA:
            box = geom.rectangle(length=self.length-self.cleave,
                height=self.width-self.cleave, position=1)
            nd.Polygon(layer=self.layer_die_user_area, points=box).\
                put(org_die.move(hcleave, hcleave))

            #DIE AREA:
            box = geom.rectangle(length=self.length,
                height=self.width, position=1)
            nd.Polygon(layer=self.layer_die, points=box).put(org_die)

            #CLEAVE AREA:
            frame = geom.frame(sizew=self.cleave, sizel=self.length,
                sizeh=self.width)
            nd.Polygon(layer=self.layer_cleave_area, points=frame).\
                put(org_die.move(hcleave, hcleave))

            #CLEAVE_LINE
            cleaveline = [
                (0, 0),
                (self.length, 0),
                (self.length, self.width),
                (0, self.width),
                (0, 0)]
            nd.Polyline(
                layer=self.layer_cleave_line,
                width=0.1,
                pathtype=2,
                points=cleaveline).put(org_die)

            #Adding coatings to nd.Cell facets:
            pinE = nd.Pin().put(org_die.move(self.length, self.width, 180))
            pinW = nd.Pin().put(org_die.move(0))
            pinN = nd.Pin().put(org_die.move(0, self.width, -90))
            pinS = nd.Pin().put(org_die.move(self.length, 0, 90))

            coatings = {
                'E': (self.coatingE, pinE),
                'W': (self.coatingW, pinW),
                'N': (self.coatingN, pinN),
                'S': (self.coatingS, pinS)}
            for side, (coat, pin) in coatings.items():
                if coat in self.coating_layers.keys():
                    if coat is None:
                        continue
                    layer = self.coating_layers[coat]
                    if side in ['E', 'W']:
                        L = self.width
                    else:
                        L = self.length
                    box = [
                        (0, 0), (self.coating_width[coat], 0),
                        (self.coating_width[coat], L), (0, L)]
                    nd.Polygon(layer=layer, points=box).put(pin)
                elif coat is not None:
                    print("""Available coating options: AR | HR | NO | DC | NA
                        AR: Anti-reflection
                        HR: High-reflection (not available in standard MPW)
                        NO: No coating (facets as-cleaved, default)
                        DC: Don't care (specify if your circuit has no waveguides on that side)
                        NA: deprecated: use 'DC' instead""")
        return C


    def get_width(self):
        """Return width of the design area"""
        return self.width


    def get_length(self):
        """Return length of the design area"""
        return self.length


    def get_cleave(self):
        """Return cleave of the design area"""
        return self.cleave


   # def get_NIO(self):
   #     """Return number of optical IOs"""
   #     return round((self.width-self.cleave-self.pitch)/self.pitch)
