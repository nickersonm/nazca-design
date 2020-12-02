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
# Author: Ronald Broeke 2017(c)
#-----------------------------------------------------------------------
"""
Convert BrightAWG coordinates.txt to nazca.txt

Reformat the AWG coordinates from the BrightAWG .txt outfile into Nazca format 
for all .txt. files in a specified directory.

    Rename port 'in' -> 'a'
    
    Rename port 'out' -> 'b'
    
    rotate ports 'in' 180 degree
    
    export with suffix .nazca.txt
"""

import os
from .netlist import Pin
from .bb_util import put_stub


def get_txt_files(dir):
    """Get all .txt. files in directory path 'dir'.

    Args:
        dir (str): root directory to scan

    Returns:
        list of str: list of filenames
    """
    root, dirs, files = next(os.walk(dir))
    txt = [os.path.join(root, f) for f in files if 's.txt' in f]
    return txt


def convert_file(file, xs, win, wout, warm=None):
    """Reformat the AWG coordinates from the OD BrightAWG .txt file into Nazca format.

    Output format: port, x, y, a, xs, w

    A header row is included in the output
    OD port names are renamed to nazca pin names, e.g. in0 -> a0, out4 -> b4.
    Also, input port will undergo a 180 degree rotation to adhere to the Nazca
    phylosophy of "outward pointing" pins.

    Args:
        file (str): filename
        xs (float): xsection of connection
        win (float): waveguide width input side
        wout (float): waveguide width output side

    Returns:
        list: list containing table with pin data
    """
def convert_file(file, xs, win, wout, warm=None):
    """Reformat the AWG coordinates from the OD BrightAWG .txt file into Nazca format.

    Output format: port, x, y, a, xs, w

    A header row is included in the output.
    OD port names are renamed to nazca pin names, e.g. in0 -> a0, out4 -> b4.
    Also, input port will undergo a 180 degree rotation to adhere to the Nazca
    philosophy of "outward pointing" pins.
    
    If the input file has the xsection and width included than those values
    are used. Otherwise the arguments  xs, win and wout are used.
        
    Args:
        file (str): filename
        xs (float): xsection of connection (if not in file)
        win (float): waveguide width input side (if not in file)
        wout (float): waveguide width output side (if not in file)

    Returns:
        list: list containing table with pin data
    """
    skiprows = 3
    replacing = [('[', ''), (']', ''), (' ', ''), (':', ','), ('\n', '')]
    output = ['port, x, y, a, xs, w']
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if i >= skiprows:
                for mapping in replacing:
                    line = line.replace(mapping[0], mapping[1])
                if 'in' in line:
                    line = line.replace('in', 'a')
                    parts = line.split(',')
                    parts[3] = str(float(parts[3])+180)
                    line = ','.join(parts)
                    width = win
                if 'out' in line:
                    line = line.replace('out', 'b')
                    width = wout
                if 'armi' in line:
                    line = line.replace('armi', 'c')
                    parts = line.split(',')
                    parts[3] = str(float(parts[3])+180)
                    line = ','.join(parts)
                    width = warm
                if 'armo' in line:
                    line = line.replace('armo', 'd')
                    width = warm
                if 'fpri' in line or 'fpro' in line:
                    parts = line.split(',')
                    parts[3] = str(float(parts[3])+180)
                    line = ','.join(parts)
                    width = None
                    xs = None
                parts = line.split(',')
                if len(parts) > 4:
                    xs = parts[4]
                if len(parts) > 5:
                    width = float(parts[5])
                line = '{},{},{}'.format(",".join(parts[:4]), xs, width)
                output.append(line)
    return output


def put_pins(file, xs, win, wout, warm=None, stubs=False):
    """Load the OD BrightAWG .txt file and place the pins in the active cell.

    This function skips the converstion to a file with Nazca coordinates.
    Can be usefull as a shortcut.

    To generate a new converted file use method 'convert' with the same
    parameters, except put_pins only accepts one file at a time.

    Args:
        file (str: filename
        xs (float): xsection of connection (if not in file)
        win (float): waveguide width input side (if not in file)
        wout (float): waveguide width output side (if not in file)

    Returns:
        None
    """
    pinlist = convert_file(file=file, xs=xs, win=win, wout=wout, warm=warm)
    #print(pinlist)
    for line in pinlist[1:]: # skip header
        name, x, y, a, xs, w = line.split(',')
        if name == 'org':
            continue
        if xs.lower() == 'none':
            xs = None
        if w.lower() == 'none':
            w = None
        else:
            w = float(w)
        Pin(name=name, width=w, xs=xs).put(x, y, a)
        if stubs:
            put_stub(name)


def convert(files, xs, win, wout):
    """Loop .txt input files and write converted nazca compatible pin output files.

    Args:
        files (str | list of str): filename(s)
        xs (float): xsection of connection (if not in file)
        win (float): waveguide width input side (if not in file)
        wout (float): waveguide width output side (if not in file)

    Returns:
        list of str: list of filenames of the converted files
    """

    if isinstance(files, str):
        files = [files]

    filenames = []
    for name in files:
        if 'coordinates.txt' in name:
            newname = name[:-4]+'_nazca.txt'
            print(newname)
            filenames.append(newname)
            updated = convert_file(name, xs, win, wout)
            with open(newname, 'w') as nf:
                for line in updated:
                    nf.write(line + '\n')
    if len(filenames) == 1:
        return filenames[0]
    else:
        return filenames
    

if __name__ == '__main__':
    path = '.'
    files = get_txt_files(path)
    xs = 'Deep'
    win = '1.5'
    wout = '2.0'
    convert(files, xs, win, wout)
