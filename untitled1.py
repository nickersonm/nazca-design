#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 15:21:26 2020

@author: rgb
"""

import nazca as nd
nd.cfg.plt_figsize = 5

nd.add_xsection('test')
nd.add_layer2xsection('lay1', layer=(2, 2))

ic = nd.interconnects.Interconnect(xs='test', PCB=True)

ic.sbend(offset=100, radius=100).put(0)

nd.export_png()