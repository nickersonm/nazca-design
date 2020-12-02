#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 21:39:19 2019

@author: rgb
"""

    # new lyp internal description
    # separate layer info in colID:
        {'layer': [layID: colID],
         'layer': [layID: colID],
         'layer': [layID: colID],
         'group':
           [groupname, colID,
            {'layer': [layID, colID],
             'layer': [layID, colID],
             'group':
               [groupname, colID,
                 {'layer': [layID, colID],
                  'layer': [layID, colID]}
               ]
             'layer': [layID, colID]}
           ]
        }