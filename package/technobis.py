# -*- coding: utf-8 -*-
"""
This is the G5 package template from Technobis IPPS.
The G5 package can accomodate ASPICs of diffferent sizes.
Minimum size is 3.9 mm x 3.9 mm.
Maximum size is 8.0 mm x 10.0 mm.

Based on and more details in:
http://www.technobis.com/files/8614/1458/3869/G5-rev1_1-design_rules_final.pdf

@author: Katarzyna Lawniczuk
"""
import pandas as pd
import nazca.cfg as cfg
import nazca as nd
import nazca.geometries as geom
import smart

print ('IPPS package loaded.')
print ('Technology options: none | smart | hhi.')

technology = {'none': (1001,    # boundary layer
                       4000,    # length
                       4000,    # hegiht
                       100,     # cleave
                       1),      # metal layer
              'smart': (60,     # boundary layer
                        4600,   # length
                        4000,   # hegiht
                        100,    # cleave
                        36),    # metal layer
              'hhi' : (60,      # boundary layer
                       6000,    # length
                       6000,    # hegiht
                       100,     # cleave
                       1),
                       }

def package(length=None, height=None, cleave=None, pitch=None,
            set_technology='none',
            drawDCbot=True, drawDCtop=True,
            setDClayer=None):
    """
    This is the G5 package from Technobis IPPS.
    """
    layer = technology[set_technology][0]
    if length is None:
        length = technology[set_technology][1]
    if height is None:
        height = technology[set_technology][2]
    if cleave is None:
        cleave = technology[set_technology][3]
    if setDClayer is None:
        setDClayer = technology[set_technology][4]
    if pitch is None:
        DCpitch = 180
    else: DCpitch = pitch

    ## Size limits
    if (length<3900):
        print('Chip is too small to be packaged with IPPS!')
    if (height<3900):
        print('Chip is too small to be packaged with IPPS!')
    if (length>10000):
        print('Chip is too large to be packaged with IPPS!')
    if (height>8000):
        print('Chip is too large to be packaged with IPPS!')

        ## DC limits
    if (length<5000):
        maxDC=21+1
    else:
        maxDC=35+1

    noDC = int((length-cleave*2)/DCpitch)
    ## checking DC limits
    if (noDC>=maxDC):
        noDC=maxDC

    def cell():
        nonlocal layer
        with nd.Cell(name='IPPS') as C:
            nd.put_pin(name='a0', connect=(0, 0, 180))
            nd.put_pin(name='b0', connect=(0, 0, 180))

            nd.polygon(layer=layer, points=geom.frame(cleave, length, height)).put((0,0))
        ## adding DC pads
            if set_technology is 'smart':

                edge=100
                DCsize=100

                if drawDCbot is True:
                    for i in range(0, noDC):
                        dcbot = smart.metal.DCpad().put((i*DCpitch+cleave+edge/2, cleave+edge/2,-90))
                        nd.put_pin(name='DCbot'+str(i), xs='Metal', connect=dcbot.pin['c0'])
                else: print('Set drawDCbot=True to get the bottom DC array.')

                if drawDCtop is True:
                    for i in range(0, noDC):
                        dctop = smart.metal.DCpad().put((i*DCpitch+cleave+edge/2, height-cleave-edge/2-DCsize,90))
                        nd.put_pin(name='DCtop'+str(i), xs='Metal', connect=dctop.pin['c0'])
                else: print('Set drawDCtop=True to get the top DC array.')



            else:
                edge=100
                DCsize=100
                layer=setDClayer
                if drawDCbot is True:
                    for i in range(0, noDC):
                        dcbot = nd.polygon(layer=layer, points=geom.square(DCsize)).put((i*DCpitch+cleave, cleave+edge/2,-90))
                        nd.put_pin(name='DCbot'+str(i), connect=(0,0,180))
                else: print('Set drawDCbot=True to get the bottom DC array.')
                if drawDCtop is True:
                    for i in range(0, noDC):
                        dctop = nd.polygon(layer=layer, points=geom.square(DCsize)).put((i*DCpitch+cleave+DCsize, height-cleave-edge/2-DCsize,90))
                        nd.put_pin(name='DCtop'+str(i), connect=(0,0,180))
                else: print('Set drawDCtop=True to get the top DC array.')

        return C
    return cell()



