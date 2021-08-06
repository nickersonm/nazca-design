#!/usr/bin/env python3
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
#-----------------------------------------------------------------------
#
#@Authors: Katarzyna Lawniczuk, Ronald Broeke 2017(c)
#@email: nazca@brightphotonics.eu


"""
**DEMO PDK building block definition.**

This module is intended to define the **building blocks** of a technology.
Building block definitions makes use of building block template functions.

This module also defines the default **Interconnect objects**,
which greatly facilitate making waveguide or metal interconnects between
building blocks.

To see the available building block and their details.
"""


import os
import nazca as nd
import nazca.pdk_template as pdk
import nazca.icons as icons


dir_path = os.path.dirname(os.path.abspath(__file__))
nd.pin2pin_drc_on()  # switch pin2pin DRC on
nd.cfg.group_connect = False  # switch off group connect in this file to save time (switch on again at bottom).
BBfunctioncalls = []  # Store PDK BB calls for test and documentation purposes

# Define the version info that will be added as annotation to each PDK cell
version = {
    'cellname': 'auto',  # 'auto' will replace 'auto' with the actual cellname
    'version': '0.5',  # version string
    'owner': 'nazca-design.org',  # building block owner
}


#==============================================================================
# Optional: Assign xsection objects to variable names for ease of use in this file.
#==============================================================================
xsShallow = nd.get_xsection('Shallow')  #: 'Shallow' xsection
xsDeep    = nd.get_xsection('Deep')     #: 'Deep' xsection
xsMetalDC = nd.get_xsection('MetalDC')  #: 'MetalDC' xsection
xsMetalRF = nd.get_xsection('MetalRF')  #: 'MetalRF' xsection


#==============================================================================
# Optional: Assign xsection parameter settings for ease of use in this file.
#==============================================================================
wshal  = xsShallow.width
Rshal  = xsShallow.radius
wdeep  = xsDeep.width
Rdeep  = xsDeep.radius
wmetdc = xsMetalDC.width
Rmetdc = xsMetalRF.radius
wmetrf = xsMetalDC.width
Rmetrf = xsMetalRF.radius
wact   = 2.0


#==============================================================================
# Optional: define (multiple) pin styles and assign them to xsections
# See cfg.py pinstylelabels for all style options
#==============================================================================
nd.add_pinstyle(
    name='default', 
    styledict={ 
        'size': 1.5,  # size of the pin shape (by default an arrow)
        'stub_length': 4.0,  # length of the stub in um
    },
)
nd.add_pinstyle(
    name='metal', 
    styledict={
        'shape': 'inv_pointer',  # shapes are defined as polygons in cfg.py
        'size': 3.0,
        'stub_length': 5.0,
    },
)
# Assign pinstyles to xsections:
# Note the pinstyle named "default" is assigned of none is specified.
xsMetalDC.pinstyle = 'metal'
xsMetalRF.pinstyle = 'metal'


#==============================================================================
# Define Interconnects (settings default to values in the xsection attributes).
#==============================================================================
shallow = nd.interconnects.Interconnect(xs='Shallow')  #: 'Shallow' interconnect
deep    = nd.interconnects.Interconnect(xs='Deep')     #: 'Deep' interconnect
metaldc = nd.interconnects.Interconnect(xs='MetalDC')  #: 'MetalDC' interconnect
metalrf = nd.interconnects.Interconnect(xs='MetalRF')  #: 'MetalRF' interconnect


#==============================================================================
#
# BUILDING BLOCKS
#
#==============================================================================


#==============================================================================
# waveguide transitions
#==============================================================================
d2s = pdk.Tp_XStransition(
    name='d2s',
    length=75,
    width=55,
    pinwidth={'a0':wdeep, 'b0':wshal},
    xs={'a0':'Deep', 'b0':'Shallow'},
    icon=icons.Tp_xsection_transition(layer1='DeepIcon', layer2='ShallowIcon')
)
BBfunctioncalls.append('d2s()')


s2d = pdk.Tp_XStransition(
    name='s2d',
    length=75,
    width=55,
    pinwidth={'a0':wshal, 'b0':wdeep},
    xs={'a0':'Shallow', 'b0':'Deep'},
    icon=icons.Tp_xsection_transition(layer1='ShallowIcon', layer2='DeepIcon')
)
BBfunctioncalls.append('s2d()')


a2s = pdk.Tp_XStransition(
    name='a2s',
    length=75,
    width=55,
    pinwidth={'a0':wact, 'b0':wshal},
    xs={'a0':'Active', 'b0':'Shallow'},
    icon=icons.Tp_xsection_transition(layer1='ActiveIcon', layer2='ShallowIcon')
)
BBfunctioncalls.append('a2s()')


s2a = pdk.Tp_XStransition(
    name='s2a',
    length=75,
    width=55,
    pinwidth={'a0':wshal, 'b0':wact},
    xs={'a0':'Shallow', 'b0':'Active'},
    icon=icons.Tp_xsection_transition(layer1='ShallowIcon', layer2='ActiveIcon')
)
BBfunctioncalls.append('s2a()')


#==============================================================================
# splitters
#==============================================================================
mmi1x1_sh = pdk.Tp_Modefilter(
    version=version,
    name='modefilter_sh',
    width=55,
    length=80,
    pinwidth={'a0':wshal,'b0':wshal},
    xs={'a0':'Shallow', 'b0':'Shallow'},
    icon=icons.Tp_icon_mmi(Nin=1, Nout=1, layer='ShallowIcon')
)
BBfunctioncalls.append('mmi1x1_sh()')


mmi1x2_sh = pdk.Tp_MMI1x2(
    name='mmi1x2_sh',
    width=55,
    length=80,
    pinwidth={'a0':wshal,'b0':wshal, 'b1':wshal},
    xs={'a0':'Shallow', 'b0':'Shallow', 'b1':'Shallow'},
    offset=2,
    icon=icons.Tp_icon_mmi(Nin=1, Nout=2, layer='ShallowIcon')
)
BBfunctioncalls.append('mmi1x2_sh()')


mmi2x2_sh = pdk.Tp_MMI2x2(
    name = 'mmi2x2_sh',
    width=60,
    length=350,
    pinwidth={'a0':wshal, 'a1':wshal, 'b0':wshal, 'b1':wshal},
    xs={'a0':'Shallow', 'a1':'Shallow', 'b0':'Shallow', 'b1':'Shallow'},
    offset=4,
    icon=icons.Tp_icon_mmi(Nin=2, Nout=2, layer='ShallowIcon')
)
BBfunctioncalls.append('mmi2x2_sh()')


mmi1x2_dp = pdk.Tp_MMI1x2(
    name='mmi1x2_dp',
    width=55,
    length=80,
    pinwidth={'a0':wdeep, 'b0':wdeep, 'b1':wdeep},
    xs={'a0':'Deep', 'b0':'Deep', 'b1':'Deep'},
    offset=2,
    icon=icons.Tp_icon_mmi(Nin=1, Nout=2, layer='DeepIcon')
)
BBfunctioncalls.append('mmi1x2_dp()')


mmi2x2_dp = pdk.Tp_MMI2x2(
    name = 'mmi2x2_dp',
    width=60,
    length=250,
    pinwidth={'a0':wdeep, 'a1':wdeep, 'b0':wdeep, 'b1':wdeep},
    xs={'a0':'Deep', 'a1':'Deep', 'b0':'Deep', 'b1':'Deep'},
    offset=4,
    icon=icons.Tp_icon_mmi(Nin=2, Nout=2, layer='DeepIcon')
)
BBfunctioncalls.append('mmi2x2_dp()')


#==============================================================================
# modulators
#==============================================================================
@pdk.hashme('abb_eopm_dc', 'length')  # Decorator to check for and resuse existing blocks
def abb_eopm_dc(length=750, contacts=2):
    """Create an electro-optic phase modulator cell.

    Args:
        length (float): length of the modulator section in um

    Returns:
        Cell: eopm element
    """
    varsinfo = {
        'length': nd.varinfo(min=100, default=length, max=2000, type='float', unit='um', doc='length of eopm'),
        'contacts': nd.varinfo(0, contacts, 6, 'int', 'null', 'number of contacts'),
    }
    nd.rangecheck(varsinfo)

    width = 60
    with nd.Cell(hashme=True) as C:
        C.version = version
        p1 = nd.Pin(name='a0', xs='Deep', width=wdeep, remark='optical').put(0, 0, 180)
        p2 = nd.Pin(name='b0', xs='Deep', width=wdeep, remark='optical').put(length, 0, 0)
        nd.connect_path(p1, p2, length)

        # x-positions of pads:
        L = 0.5*(length - wmetdc)
        xpos = []
        if contacts > 1:
            for n in range(contacts):
                xpos.append(0.5*wmetdc + L*2*n/(contacts-1))
        elif contacts == 1:
            xpos = [0.5*length]
        elif contacts == 0:
            pass
        for n, x in enumerate(xpos):
            pinname = 'c'+str(n)
            p1 = nd.Pin(name=pinname, xs='MetalDC', width=wmetdc).put(x, 0, 90)
            nd.put_stub(pinname)

        #nd.Pin(name='c0', xs='MetalDC', width=wmetdc).put(0.5*length, 5.0, 90)
        pdk.put_stub(['a0', 'b0', 'c0'])
        pdk.put_boundingbox('org', length, width, align='lc')
        icons.icon_strt(length, width, bufx=0, bufy=10, layer='MetalDCIcon').put('cc', 0, 'cc')
        icons.icon_strt(length, width, bufx=10, bufy=-8, layer='DeepIcon').put('cc', 0, 'cc')
    return C
BBfunctioncalls.append('abb_eopm_dc()')


#==============================================================================
# photodiodes
#==============================================================================
pd = pdk.Tp_PhotoDetector(
    name='pd_dp',
    length=50,
    width=55,
    pinwidth={'a0':wdeep, 'b0':wdeep, 'c0':wmetdc},
    xs={'a0':'Deep', 'b0':'Deep', 'c0':'MetalDC'},
    icon=icons.Tp_icon_strt(bufx=0, bufy=12.5, layer='MetalDCIcon')
)
BBfunctioncalls.append('pd()')


#==============================================================================
# laser elements
#==============================================================================
@pdk.hashme('soa', 'length', 'pad')
def soa(length=100, pad=True):
    """Create a SOA cell.

    Args:
        length (float): length of gain section in um

    Returns:
        Cell: SOA element
    """
    if pad:
        width = 200
    else:
        width = 50

    varsinfo = {
        'length': nd.varinfo(10, length, 2000, 'float', 'um', 'length of soa'),
    }
    nd.rangecheck(varsinfo)

    h = -70
    with nd.Cell(hashme=True) as C:
        C.version = version
        p1 = nd.Pin(name='a0', xs='Active', width=wact).put(0, h, 180)
        p2 = nd.Pin(name='b0', xs='Active', width=wact).put(length, h, 0)
        nd.connect_path(p1, p2, length)
        nd.Pin(name='c0', xs='MetalDC', width=wmetdc).put(0.5*length, 5.0, 90)
        pdk.put_stub(['a0', 'b0', 'c0'])
        pdk.put_boundingbox('org', length, width, align='lc')
        icons.icon_strt(length, width, bufx=0, bufy=10, layer='MetalDCIcon').put('cc', 0, 'cc')
        icons.icon_strt(length, width, bufx=10, bufy=-8, layer='ActiveIcon').put('cc', 0, h, 'cc')
    return C
BBfunctioncalls.append('soa()')


@pdk.hashme('soa_sh', 'length')
def soa_sh(length=100):
    """Create a SOA cell with Shallow io.

    Args:
        length (float): length of gain section in um

    Returns:
        Cell: SOA element
    """
    varsinfo = {
        'length': nd.varinfo(10, length, 2000, 'float', 'um', 'length of soa'),
    }
    nd.rangecheck(varsinfo)

    with nd.Cell(hashme=True) as C:
        C.version = version
        s2a_ = s2a().put(0)
        soa_ = soa(length).put()
        a2s_ = a2s().put()
        nd.Pin(name='a0', pin=s2a_.pin['a0']).put()
        nd.Pin(name='b0', pin=a2s_.pin['b0']).put()
        nd.Pin(name='c0', pin=soa_.pin['c0']).put()

        #pdk.put_stub(['a0', 'b0', 'c0'])
    return C
BBfunctioncalls.append('soa_sh()')


@pdk.hashme('dbr', 'length', 'pitch', 'duty_cycle')
def dbr(length=100, pitch=0.250, duty_cycle=0.5):
    """Create a DBR cell.

    Args:
        length (float): length of DBR section in um
        pitch (float): pitch (full period) of dbr in um
        duty_cycle (float): duty cycle of grating (default = 0.5)

    Returns:
        Cell: dbr element
    """
    width = 200
    h = -70
    varsinfo = {
        'length': nd.varinfo(min=0.0, default=length, max=1000.0, type='float', unit='um', doc='length of dbr'),
        'pitch': nd.varinfo(0.1, pitch, 0.5, 'float', 'um', 'length of soa '),
        'duty_cycle': nd.varinfo(0.0, duty_cycle, 1.0, 'float', 'um', 'length of soa'),
    }
    nd.rangecheck(varsinfo)

    icon = icons.Tp_icon_grating('ActiveIcon')
    with nd.Cell(hashme=True) as C:
        C.version = version
        p1 = nd.Pin(name='a0', xs='Active', width=wact).put(0, h, 180)
        p2 = nd.Pin(name='b0', xs='Active', width=wact).put(length, h, 0)
        nd.Pin(name='c0', xs='MetalDC', width=wmetdc).put(0.5*length, 5.0, 90)
        nd.connect_path(p1, p2, length)

        pdk.put_stub(['a0', 'b0', 'c0'])
        pdk.put_boundingbox('org', length, width, align='lc')
        icons.icon_strt(length, width, bufx=0, bufy=10, layer='MetalDCIcon').put('cc', 0, 'cc')
        icons.icon_strt(length, width, bufx=10, bufy=-8, layer='ActiveIcon').put('cc', 0, h, 'cc')
        icon(length, width, bufx=5, bufy=-40).put('cc', 0, h, 'cc')
    return C
BBfunctioncalls.append('dbr()')


@pdk.hashme('phase_shifter', 'length')
def phase_shifter(length=100):
    """Create a current injection phase shifting cell.

    Args:
        length (float): length of phase section in um

    Returns:
        Cell: phase-shifter element
    """
    nd.rangecheck({'length': nd.varinfo(10, length, 2000, 'float', 'um', 'length of phase_shifter'),})

    width = 200
    h = -70
    with nd.Cell(hashme=True) as C:
        C.version = version
        p1 =nd.Pin(name='a0', xs='Active', width=wact).put(0, h, 180)
        p2 = nd.Pin(name='b0', xs='Active', width=wact).put(length, h, 0)
        nd.connect_path(p1, p2, length)
        nd.Pin(name='c0', xs='MetalDC', width=wmetdc).put(0.5*length, 5.0, 90)
        pdk.put_stub(['a0', 'b0', 'c0'])
        pdk.put_boundingbox('org', length, width, align='lc')
        icons.icon_strt(length, width, bufx=0, bufy=10, layer='MetalDCIcon').put('cc', 0, 'cc')
        icons.icon_strt(length, width, bufx=10, bufy=-8, layer='ActiveIcon').put('cc', 0, h, 'cc')
    return C
BBfunctioncalls.append('phase_shifter()')


def Tp_isolation(name='isolation', length=100, xs=None, width=None, layer=None):
    """Template funtion."""
    @pdk.hashme(name, 'length')
    def isolation(length=length):
        """Create a p-isolation cell.

        Args:
            length (float): length of the isolation section in um

        Returns:
            Cell: isolation element
        """
        nd.rangecheck({'length': nd.varinfo(0, length, 2000, 'float', 'um', 'length of isolation'),})

        bb_width = 20
        with nd.Cell(hashme=True) as C:
            C.version = version
            p1 = nd.Pin(name='a0', xs=xs, width=width).put(0, 0, 180)
            p2 = nd.Pin(name='b0', xs=xs, width=width).put(length, 0, 0)
            nd.connect_path(p1, p2, length)
            pdk.put_stub()
            pdk.put_boundingbox('org', length, bb_width)
            icons.icon_strt(length, bb_width, bufx=4, bufy=-8, layer=layer).put('cc', 0, 'cc')
        return C
    return isolation
isolation_sh = Tp_isolation(name='isolation_sh', xs='Shallow', width=wshal, layer='ShallowIcon')
BBfunctioncalls.append('isolation_sh()')
isolation_act = Tp_isolation(name='isolation_act', xs='Active', width=wact, layer='ActiveIcon')
BBfunctioncalls.append('isolation_act()')


#==============================================================================
# io_optical
#==============================================================================
#io = pdk.Tp_IOtrapezoid(
#    name = 'io_sh',
#    length=100,
#    width=4.5,
#    angle=-7,
#    pinwidth={'a0':wshal, 'b0':wshal},
#    xs={'a0':'Shallow', 'b0':'Shallow'})
#BBfunctioncalls.append('io()')


#==============================================================================
# io_electrical
#==============================================================================
pad_dc = pdk.Tp_DCpad_rectangle(
    name='pad_dc',
    length=150,
    pinwidth={'c0':wmetdc},
    xs={'c0':'MetalDC'},
    icon=icons.Tp_icon_rounded_pad(bufx=20, bufy=20, layer=10)
)
BBfunctioncalls.append('pad_dc()')


pad_rf = pdk.Tp_RFpad(
    name='pad_rf',
    length=80,
    width=80,
    pinwidth={'c0':wmetrf},
    xs={'c0':'MetalRF'},
    icon=icons.Tp_icon_strt(bufx=0, bufy=5, layer=11)
)
BBfunctioncalls.append('pad_rf()')


# =============================================================================
# Composite building blocks
# =============================================================================
@pdk.hashme('eopm_dc', 'length', 'pads', 'sep')
def eopm_dc(length=1000, pads=False, sep=10):
    """EOPM with pad option.

    Args:
        length (float): modulator length in um
        pads (bool): draw pads if True (default = False)
        sep (float): separation of pad bbox from eopm stub.

    Returns:
        Cell: eopm element
    """
    # Define range checking:
    varinfo = {
        'length': nd.varinfo(min=0.0, default=length, max=2000.0, type='float', unit='um', doc='length of eopm'),
        'pads': nd.varinfo(min=None, default=pads, max=None, type='bool', unit='null', doc='flag to add pads or not'),
        'sep': nd.varinfo(min=0.0, default=sep, max=500.0, type='float', unit='um', doc='separtions of pads from the soa'),
    }
    nd.rangecheck(varinfo)

    contacts = 2
    with nd.Cell() as C:
        C.version = version
        E1 = abb_eopm_dc(length=length, contacts=contacts).put(0)
        if pads:
            for i in range(contacts):
                pname = 'c'+str(i)
                E2 = pad_dc().put('rc', E1.pin[pname].move(sep))
                metaldc.strt_p2p(pin1=E1.pin[pname], pin2=E2).put()
                E2.raise_pins(['c0'], [pname])
            E1.raise_pins(['a0', 'b0'])
        else:
            E1.raise_pins()
    return C
BBfunctioncalls.append('eopm_dc()')


@pdk.hashme('PD_dc', 'length', 'pads')
def pd_dc(length=100, pads=False, dy=0):
    """"PhotoDiode with pad option.

    Args:
        length (float): PD length in um
        pads (bool): draw pads if True (default = False)
        sep (float_: separation of pad bbox from pin stub.

    Returns:
        Cell: PD element
    """
    # Define range checking:
    varinfo = {
        'length': nd.varinfo(min=0.0, default=length, max=2000.0, type='float', unit='um', doc='length of eopm'),
        'pads': nd.varinfo(min=None, default=pads, max=None, type='bool', unit='null', doc='flag to add pads or not '),
        'dy': nd.varinfo(min=0.0, default=dy, max=500.0, type='float', unit='um', doc='separtions of pads from the soa'),
    }
    nd.rangecheck(varinfo)

    with nd.Cell() as C:
        C.version = version
        E1 = pd(length=length).put(0)
        if pads:
            E2 = pad_dc().put('rc', E1.pin['c0'].move(dy))
            metaldc.strt_p2p(pin1=E1.pin['c0'], pin2=E2).put()
            E2.raise_pins(['c0'])
            E1.raise_pins(['a0'])
        else:
            E1.raise_pins()
    return C
BBfunctioncalls.append('pd_dc()')


@pdk.hashme('io', 'shape', 'width', 'bend', 'deep')
def io(shape=None, width=3.0, bend=False, deep=False):
    """Create an IO cells that for edge coupling waveguide interfaces.

    Args:
        shape (str): shape of IO. values: None | 'simple' | 'tapered' | 'angled'
        bend (bool): make the IO connection horizontal (True) for connecting to layout
        deep (bool): add a transition to a deeply-etched waveguides (default = False)

    Returns:
        Cell: io element
    """
    options = ['simple', 'tapered', 'angled']
    if shape is None:
        shape = 'simple'
    if shape not in options:
        print("Error: io shape '{}'not recognized. Options are: {}".\
            format(shape, options))

    pin_indent = 0
    length = 100-pin_indent

    if shape == 'tapered':
        with nd.Cell(hashme=True) as C:
            C.version = version
            s = shallow.strt(length, width, arrow=False).put(0)
            iowg = shallow.taper(length, width, 3.0, arrow=False).put()
            if deep:
                shallow.strt(10, arrow=False).put()
                iowg = s2d().put()
            nd.Pin('a0', pin=s.pin['a0']).put()
            nd.Pin('b0', pin=iowg.pin['b0']).put()

    elif shape == 'angled':
        with nd.Cell(hashme=True) as C:
            C.version = version
            if not bend:
                s = io(width=4.5).put(0)
                iowg = shallow.taper(length, 4.5, 3.0, arrow=False).put()
                if deep:
                    shallow.strt(10, arrow=False).put()
                    nd.make_pincell().put()
                    iowg = s2d().put()
                else:
                    nd.make_pincell().put()
                nd.Pin('a0', pin=s.pin['a0']).put()
                nd.Pin('b0', pin=iowg.pin['b0']).put()

            else:
                s = io(width=4.5).put(0)
                shallow.taper(length, 4.5, 3.0, arrow=False).put()
                iowg = shallow.bend(angle=-7, arrow=False).put()
                if deep:
                    shallow.strt(10, arrow=False).put()
                    iowg = s2d().put()
                    nd.make_pincell().put()
                else:
                    nd.make_pincell().put()
                nd.Pin('a0', pin=s.pin['a0']).put()
                nd.Pin('b0', pin=iowg.pin['b0']).put()

    else: #shape == None or shape == 'simple':
        with nd.Cell(hashme=True) as C:
            C.version = version
            if not deep:
                iowg = shallow.strt(length, width=width).put(0)
                nd.Pin('a0', pin=iowg.pin['a0']).put()
                out = nd.Pin('b0', pin=iowg.pin['b0']).put()
            else:
                shallow.strt(length+10, width=width).put(0)
                iowg = s2d().put()
                nd.Pin('a0', pin=iowg.pin['a0']).put()
                out = nd.Pin('b0', pin=iowg.pin['b0']).put()
            nd.make_pincell().put(out)

    return C


@pdk.hashme('dbr_laser', 'Ldbr1', 'Ldbr2', 'Lsoa', 'Lpm')
def dbr_laser(Ldbr1=50, Ldbr2=500, Lsoa=750, Lpm=70):
    """Create a parametrized dbr laser building block."""
    Liso = 20
    with nd.Cell(hashme=True) as laser:
        laser.version = version
        #create an isolation cell for reuse
        iso = isolation_act(length=Liso)

        #draw the laser
        s2a1 = s2a().put(0)
        iso.put()
        dbr1 = dbr(length=Ldbr1).put()
        iso.put()
        soa1 = soa(length=Lsoa).put()
        iso.put()
        phase1 = phase_shifter(length=Lpm).put()
        iso.put()
        dbr2 = dbr(length=Ldbr2).put()
        iso.put()
        a2s1 = a2s().put()

        # add pins to the laser building block
        nd.Pin('a0', pin=s2a1.pin['a0']).put()
        nd.Pin('b0', pin=a2s1.pin['b0']).put()
        nd.Pin('c0', pin=dbr1.pin['c0']).put()
        nd.Pin('c1', pin=soa1.pin['c0']).put()
        nd.Pin('c2', pin=phase1.pin['c0']).put()
        nd.Pin('c3', pin=dbr2.pin['c0']).put()
        pdk.put_stub()

        length, y, a = nd.diff(laser.pin['a0'], laser.pin['b0'])
        pdk.put_boundingbox('org', length=abs(length), width=200, align='lb',
            move=(0, -30, 0))
    return laser
BBfunctioncalls.append('dbr_laser()')


@pdk.hashme('mzi', 'length', 'sep')
def mzi(length=1000, sep=50):
    with nd.Cell(hashme=True) as mziBB:
        mziBB.autobbox = True
        mziBB.version = version
        eopm = eopm_dc(length=length, pads=True, sep=40)
        mmi_left  = mmi2x2_dp().put('lc')
        deep.sbend(offset=sep).put(mmi_left.pin['b0'])
        eopm_top = eopm.put()
        deep.sbend(offset=-sep, length1=0.01).put()
        mmi_right = mmi2x2_dp().put()

        deep.sbend(offset=-sep).put(mmi_left.pin['b1'])
        eopm_bot = eopm.put(flip=True)
        deep.sbend_p2p(pin2=mmi_right.pin['a1']).put()

        nd.Pin('a0', pin=mmi_left.pin['a0']).put()
        nd.Pin('a1', pin=mmi_left.pin['a1']).put()
        nd.Pin('b0', pin=mmi_right.pin['b0']).put()
        nd.Pin('b1', pin=mmi_right.pin['b1']).put()
        nd.Pin('c0', pin=eopm_top.pin['c0']).put()
        nd.Pin('c1', pin=eopm_top.pin['c1']).put()
        nd.Pin('d0', pin=eopm_bot.pin['c0']).put()
        nd.Pin('d1', pin=eopm_bot.pin['c1']).put()
        pdk.put_stub()

        bb_length, y, a = nd.diff(mziBB.pin['a0'], mziBB.pin['b0'])
        pdk.put_boundingbox('org', length=abs(bb_length), width=2*sep+428)
    return mziBB
BBfunctioncalls.append('mzi()')


#==============================================================================
# GDS based BBs
#==============================================================================
with nd.Cell('awg1x4') as awg1x4:
    awg1x4.version = version
    awg = nd.load_gds(
        filename= os.path.join(dir_path, 'gdsBB', 'AWG_1x4.gds'),
        cellname='awg',
        newcellname='demo_awg',
        layermap={1:3, 20:'BB'},
    )
    awg.put()

    nd.Pin('a0', xs='Deep', width=wdeep).put((9.184851e-016, -5, -90))
    nd.Pin('b0', xs='Deep', width=wdeep).put((273.63575, -4.5798919, -85.569401))
    nd.Pin('b1', xs='Deep', width=wdeep).put((269.89849, -4.9532407, -88.523134))
    nd.Pin('b2', xs='Deep', width=wdeep).put((266.1423, -4.9532407, -91.476866))
    nd.Pin('b3', xs='Deep', width=wdeep).put((262.40504, -4.5798919, -94.430599))

    # lengths are just example values:
    nd.connect_path(awg1x4.pin['a0'], awg1x4.pin['b0'], 300, [(500, 140)])
    nd.connect_path(awg1x4.pin['a0'], awg1x4.pin['b1'], 500, [(600, 140)])
    nd.connect_path(awg1x4.pin['a0'], awg1x4.pin['b2'], 700, [(700, 140)])
    nd.connect_path(awg1x4.pin['a0'], awg1x4.pin['b3'], 900, [(800, 140)])
    
    nd.put_stub(['a0', 'b0', 'b1', 'b2', 'b3'])
BBfunctioncalls.append('awg1x4')


def BBcells():
    return [eval(F) for F in BBfunctioncalls]


#==============================================================================
# Create groups to provide a structure to put building blocks (BB) under.
# This is optional, though useful in navigation BBs.
#==============================================================================
bb_connectors     = pdk.Functional_group() #: inter-connectors
bb_transitions    = pdk.Functional_group() #: xsection transistions
bb_splitters      = pdk.Functional_group() #: power splitter and combiners
bb_combiners      = pdk.Functional_group() #: power splitter and combiners
bb_amplifiers     = pdk.Functional_group() #: SOAs
bb_modulators     = pdk.Functional_group() #: modulators
bb_photodiodes    = pdk.Functional_group() #: photo diodes
bb_mode_adapters  = pdk.Functional_group() #: spot size convertors
bb_io_electrical  = pdk.Functional_group() #: electrical IO (pads)
bb_io_optical     = pdk.Functional_group() #: optical IO (edge or surface)
bb_mode_filters   = pdk.Functional_group() #: mode filters
bb_awgs           = pdk.Functional_group() #: AWG multiplexers
bb_tapers         = pdk.Functional_group() #: tapers between different widths


bb_connectors.shallow = shallow
bb_connectors.deep = deep
bb_connectors.metaldc = metaldc
bb_connectors.metalrf = metalrf

bb_amplifiers.soa = soa
bb_awgs.awg1x4 = awg1x4
bb_io_electrical.pad_dc = pad_dc
bb_io_electrical.pad_rf = pad_rf
bb_io_optical.io = io

bb_modulators.eopm_dc = eopm_dc
bb_photodiodes.pd = pd

bb_mode_filters.mmi1x1_sh = mmi1x1_sh
bb_splitters.mmi1x2_sh = mmi1x2_sh
bb_splitters.mmi2x2_sh = mmi2x2_sh
bb_splitters.mmi1x2_dp = mmi1x2_dp
bb_splitters.mmi2x2_dp = mmi2x2_dp
bb_transitions.s2d = s2d
bb_transitions.d2s = d2s

nd.cfg.group_connect = True
