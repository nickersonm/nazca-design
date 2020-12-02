#!/usr/bin/env python3
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
# -*- coding: utf-8 -*-
# 2017 (c)  Ronald Broeke
# 2016-2018 (c) Xaveer Leijtens / csv2lyp

"""
Module to parse Klayout lyp files into csv via DataFrames: lyp2csv

Generate a matplotlib colormap from the csv table: csv2colormap
"""

"""
Structure in .lyp XML file:


layer-properties-tabs -> present in case of multiple tabs only
    layer-properties -> represent a tab view in Klayout
        properties
            GROUP
                LAYER2
                LAYER2
            GROUP
                GROUP
                    LAYER2
                    LAYER2
                GROUP
                    LAYER2
                    LAYER2
            LAYER1
            lAYER1

LAYER1:
    <properties>
    <...lyp_attr>
    <name>   -> layer name
    <source> -> layer/datatype

LAYER2:
    <group-members>
        <...lyp_attr>
        <name>   -> layer name
        <source> -> layer/datatype

GROUP:
    <group-members>
        <...lyp_attr>
        <name>   -> group-name
        <source> -> */*
"""

from collections import defaultdict
import pandas as pd
import xml.etree.ElementTree as Tree
import nazca.cfg as cfg
from nazca.logging import logger

lyp_attr   = ['layer', 'datatype', 'source', 'fill-color', 'frame-color',
    'frame-brightness', 'fill-brightness', 'dither-pattern', 'valid',
    'visible', 'transparent', 'width', 'marked', 'animation']
nazca_attr = ['layer', 'datatype', 'name', 'fill_color', 'frame_color',
    'frame_brightness', 'fill_brightness', 'dither_pattern', 'valid',
    'visible', 'transparent', 'width', 'marked', 'animation']


path = ''
doPrint = False
#==============================================================================
# lyp2csv
#==============================================================================
tabdict = defaultdict(list)
depth=0
def _parse_properties(lev1, infolevel=0):
    """Parse lyp tags <properties> and <group-member> levels."""
    global tabdict
    global depth

    depth += 1
    tabdict['depth'].append(depth)
    for lev2 in lev1:
        tag = lev2.tag
        value = lev1.find(tag).text
        if infolevel > 2:
            if tag == 'group-members': # remove linefeed
                value = ''
            print("{}{}: {}".format('  '*depth, tag, value))
        if tag == 'group-members':
            _parse_properties(lev2, infolevel=infolevel)
        else:
            if tag == 'source':
                value = value[:-2].split('/')
                tabdict['layer'].append(value[0])
                tabdict['datatype'].append(value[1])
            else:
                tabdict[tag].append(value)
    depth -= 1
    return None


def _parse_tab(lev1, infolevel=0):
    """Parse a lyp-file tab section.

    Rename Klayout attribute names to Nazca names
    and save the tab content to a csv file.

    Returns:
        str, Pandas DataFrame: tabname, Color data
    """
    global tabdict, path

    tabdict = defaultdict(list) #clear global dict for new tab
    tabname = 'nazca'
    for lev2 in lev1: # 'layer-properties'
        tag = lev2.tag
        if tag == 'properties':
            if infolevel > 1:
                print(tag)
            _parse_properties(lev2, infolevel=infolevel)
        elif tag == 'name': # tab name
            if lev1.find(tag).text is not None:
                tabname = lev1.find(tag).text
                if infolevel > 1:
                    print('tabname: {}'.format(tabname))

    df = pd.DataFrame.from_dict(tabdict)
    df.rename(columns=dict(zip(lyp_attr, nazca_attr)), inplace=True)

    #csv_filename = "{}_{}.csv".format(filenamebase, tabname)
    return tabname, df[['depth']+nazca_attr]
    #if infolevel > 0:
    #    print("Exported tab '{}' to file '{}'.".format(tabname, csv_filename))
    #return tabname, csv_filename


def lyp2csv(lypfile, infolevel=0):
    """Convert a Klayout .lyp file into a .csv Nazca color table.

    This function makes it possible to use a Klayout color definition (in xml)
    inside Nazca, e.g. for Matplotlib output.
    A .lyp file may contain multiple tabs
    and each tab is exported to a separate .csv file.

    Note that display method in Klayout and Matplotlib are not the same.
    Hence, a layout will have a different "style" to it.
    for things like:

        * tranparency and opaqueness
        * visualization of edges
        * stipple

    Args:
        lypfile (str): name of Klayout .lyp file

    Returns:
        dict: dictionary {<tabname>: <csv_filename>}
    """

    #tree = Tree.parse(_clean_lyp(filename))

    tree = Tree.parse(lypfile)
    root = tree.getroot()

    tabs = {}
    # Parse tags <layer-properties-tabs> and <layer-properties>
    # if >1 tabs the <layer-properties-tabs> tag is present with a tab name.

    if infolevel > 0:
        print('file: {}'.format(lypfile))

    if root.tag == 'layer-properties-tabs':
        if infolevel > 0:
            print('Multiple tabs found')
        for lev1 in root: # <layer-properties>
            tab, df = _parse_tab(lev1, infolevel=infolevel)
            csv_filename = "{}_{}.csv".format(lypfile[:-4], tab)
            df.to_csv(csv_filename, index=False, na_rep='')
            tabs[tab] = csv_filename
            if infolevel > 0:
                print('root:', lev1.tag)
                print("Exported tab '{}' to file '{}'.".format(tabname, csv_filename))
    else: # <layer-properties>
        if infolevel > 0:
            print('Single tab found')
            tab, df = _parse_tab(lev1, infolevel=infolevel)
            csv_filename = "{}_{}.csv".format(lypfile[:-4], tab)
            df.to_csv(csv_filename, index=False, na_rep='')
            tabs[tab] = csv_filename

    return tabs


def lyp2pd(lypfile, infolevel=0):
    """Convert a Klayout .lyp file into a Nazca color table.

    This function makes it possible to use a Klayout color definition (in xml)
    inside Nazca, e.g. for Matplotlib output.
    A .lyp file may contain multiple tabs
    and each tab is exported to a separate pandas dataframe.

    Note that display method in Klayout and Matplotlib are not the same.
    Hence, a layout will have a different "style" to it.
    for things like:

        * tranparency and opaqueness
        * visualization of edges
        * stipple

    Args:
        lypfile (str): name of Klayout .lyp file

    Returns:
        dict: dictionary {<tabname>: <tab_dataframe>}
    """

    #tree = Tree.parse(_clean_lyp(filename))

    tree = Tree.parse(lypfile)
    root = tree.getroot()

    tabs = {}
    # Parse tags <layer-properties-tabs> and <layer-properties>
    # if >1 tabs the <layer-properties-tabs> tag is present with a tab name.

    if infolevel > 0:
        print('file: {}'.format(lypfile))

    if root.tag == 'layer-properties-tabs':
        if infolevel > 0:
            print('Multiple tabs found')
        for lev1 in root: # <layer-properties>
            tab, df = _parse_tab(lev1, infolevel=infolevel)
            tabs[tab] = df
            if infolevel > 0:
                print('root:', lev1.tag)
    else: # <layer-properties>
        if infolevel > 0:
            print('Single tab found')
            tab, df = _parse_tab(lev1, infolevel=infolevel)
            tabs[tab] = df

    boolmap = {'false': False, 'true': True}
    bools = ['valid', 'visible', 'transparent', 'marked']
    for tab in tabs:
        for boolcolm in bools:
            tabs[tab][boolcolm] = tabs[tab][boolcolm].map(boolmap)

    return tabs



#==============================================================================
# csv2lyp
#==============================================================================
class Lypdata:
    """Class to export layer infomation to a Klayout lyp file.

    A lyp file is a layer properties file.
    """
    def __init__(self):
        """Initialize a Lypdate object.

        Returns:
            None
        """
        self.level = 0
        self.stack = []
        self.xml = []
        self.xml_tag('?xml version="1.0" encoding="utf-8"?')


    def xml_tag(self, tag, close=False):
        """Append xml tag.

        Args:
            tag (str):
            close (bool):

        Returns:
            None
        """
        self.xml.append('{}<{}{}>'.format(self.level*" ", close*"/", tag))


    def xml_value(self, tag, val=None):
        """Append xml value.

        Args:
            tag (str):
            val (): value

        Returns:
            None
        """
        if not val:
            self.xml.append('{}<{}/>'.format(self.level*" ", tag))
        else:
            self.xml.append('{}<{}>{}</{}>'.format(self.level*" ",tag,val,tag))


    def push(self, tag):
        self.xml_tag(tag)
        self.stack.append(tag)
        self.level += 1


    def pop(self):
        self.level -= 1
        self.xml_tag(self.stack.pop(), close=True)


    def layer2xml(self, row):
        """Append layer entry.

        Returns:
            None
        """
        for field in nazca_attr[3:]:
            self.xml_value(field.replace('_','-'), row[field])
        if row['layer'] == '*':
            src = '*/*@*'
        else:
            src = "{}/{}@1".format(row['layer'], row['datatype'])
        self.xml_value('name', row['name'])
        self.xml_value('source', src)


    def output(self, filename):
        """Save lyp file.

        Args:
            filename (str): output filename

        Returns:
            None
        """
        for i in range(len(self.stack)):
            self.pop() # Close remaining tags.
        with open(filename, 'w') as f:
            f.write('\n'.join(self.xml))
        print('wrote {}'.format(filename))


def csv2lyp(csvfiles, lypfile):
    """Write klayout .lyp layer properties file.

    Args:
        csvfiles (dict): dictionary {tab-name: filename} which contains all
            csv files to be written to a single lyp file. Each tab-name is
            the name of the tab in the layers view in KLayout and each
            filename is the filename of the csv file to be read.
        lypfile (str): filename of the .lyp file. This file will be created
            if it does not yet exist. If it already exists, it will be
            overwritten.

    Returns:
        None
    """
    lyp = Lypdata()

    if len(csvfiles) > 1:
        lyp.push('layer-properties-tabs')

    for tabname, filename in csvfiles.items():
        depth = 0
        df = pd.read_csv(filename, dtype=str, keep_default_na=False)
        df.loc[df.width == 'na', 'width'] = None # remove na from csv
        lyp.push('layer-properties') # Start of (new) tab.

        for index, row in df.iterrows():
            curdepth = int(row['depth'])
            for i in range(depth - curdepth + 1):
                lyp.pop() # Close previous
            depth = curdepth
            if depth == 1: # Oddity of lyp file.
                grouptag = 'properties'
            else:
                grouptag = 'group-members'
            lyp.push(grouptag) # Start of group
            lyp.layer2xml(row)
        while curdepth>1:
            lyp.pop()
            curdepth-=1            
        lyp.pop()
        lyp.xml_value('name', tabname)
        lyp.pop()
    lyp.output(lypfile)

def pd2lyp(tabs, lypfile):
    """Write klayout .lyp layer properties file.

    Args:
        csvfiles (dict): dictionary {tab-name: filename} which contains all
            csv files to be written to a single lyp file. Each tab-name is
            the name of the tab in the layers view in KLayout and each
            filename is the filename of the csv file to be read.
        lypfile (str): filename of the .lyp file. This file will be created
            if it does not yet exist. If it already exists, it will be
            overwritten.

    Returns:
        None
    """
    lyp = Lypdata()

    boolmap = {False: 'false', True: 'true'}
    bools = ['valid', 'visible', 'transparent', 'marked']
    for tab in tabs:
        for boolcolm in bools:
            tabs[tab][boolcolm] = tabs[tab][boolcolm].map(boolmap)

    if len(tabs) > 1:
        lyp.push('layer-properties-tabs')

    for tabname, df in tabs.items():
        depth = 0
        df.loc[df.width == 'na', 'width'] = None # remove na from csv
        lyp.push('layer-properties') # Start of (new) tab.

        for index, row in df.iterrows():
            curdepth = int(row['depth'])
            for i in range(depth - curdepth + 1):
                lyp.pop() # Close previous
            depth = curdepth
            if depth == 1: # Oddity of lyp file.
                grouptag = 'properties'
            else:
                grouptag = 'group-members'
            lyp.push(grouptag) # Start of group
            lyp.layer2xml(row)
        while curdepth>1:
            lyp.pop()
            curdepth-=1            
        lyp.pop()
        lyp.xml_value('name', tabname)
        lyp.pop()
    lyp.output(lypfile)




def nazca2csv(colors, filename):
    """Save nazca colormap to csv file.

    Args:
        colors (DataFrame): color information
        filename (str): name of csv output file

    Returns:
        None
    """
    #colors['depth'] = 1
    columns = [
        'depth',
        'layer',
        'datatype',
        'name',
        'fill_color',
        'frame_color',
        'frame_brightness',
        'fill_brightness',
        'dither_pattern',
        'valid',
        'visible',
        'transparent',
        'width',
        'marked',
        'animation']
#        'alpha',
    boolmap = {False: 'false', True: 'true'}
    bools = ['valid', 'visible', 'transparent', 'marked']
    for boolcolm in bools:
        colors.loc[:,boolcolm] = colors.loc[:,boolcolm].map(boolmap)
    colors[columns].to_csv(filename, index=False)
    return None


def nazca2lyp(filename):
    """Export the Nazca color DataFrame to .lyp file.

    Generate a csv file from the Nazca color DataFrame and a lyp file
    from that csv file

    Args:
        filename (str): basename of the output csv and lyp files

    Returns:
        None
    """
    csv_filename = '{}.csv'.format(filename)
    lyp_filename = '{}.lyp'.format(filename)

    nazca2csv(cfg.colors, csv_filename)
    csv2lyp({'nazca': csv_filename}, lyp_filename)
    return None


if __name__ == '__main__':
    from pprint import pprint

    fab = 'hhi'
    if fab == 'demo':
        path = 'demofab/'
        file = path + 'demofab_klayout_colors.lyp'
    elif fab =='smart':
        path = '../../pdk/smart/smart/'
        file = path + 'smart-nazca.lyp'
    elif fab =='hhi':
        path = '../../pdk/hhi/hhi/'
        file = path + 'hhi-3.lyp'
        tolyp = '../../pdk/hhi/hhi/table_colors_HHI-full.csv'
    tabs = lyp2csv(file, infolevel=3)
    pprint(tabs)

    f1 = 'demofab/table_colors_Light.csv'
    f2 = 'demofab/table_colors_Light.csv'

    csv2lyp({'Light1':f1}, 'demo-light1.lyp')
    csv2lyp({'Light1':f1, 'Light2':f2}, 'demo-light2.lyp')
    csv2lyp({'Light1':tolyp}, '../../pdk/hhi/hhi/demo-tolyp.lyp')

