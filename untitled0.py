#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 22:16:57 2020

@author: rgb
"""


import nazca as nd
import hhi

e1 = hhi.e1700.strt().put()
nd.e1700.bend_strt_bend_p2p(pin1=e1.pin['a0'], pin2=e1.pin['b0']).put()

nd.pathfinder(start=e1.pin['a0'])

