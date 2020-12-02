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
#
# Reference to low-level routines, should you want them.
#

"""Low-level Nazca routines."""

__all__ = ['GDSII_stream', 'GDSII_cell', 'GDSII_element', 'GDSII_record',
        'get_cell_annotation', 'get_cell_polyline', 'get_cell_polygon',
        'parameters_to_string', 'string_to_parameters']

from .gds_import import *
from .util import *
