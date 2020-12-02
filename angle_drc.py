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
# Test replacement of gds cells from one file with those from another.
# This will be used for replacing black-box building blocks with their real
# implementation.
#
# (c) 2019 Ronald Broeke
#==============================================================================
"""
DRC on instance attributes angle and/or flip state.
"""

import os
from collections import OrderedDict, Counter
import math as m
import datetime
from pprint import pprint
import nazca as nd

drc_filename = 'drc.log'
blocklist_filename = 'bblist.log'

items = [] # DRC items
blocks = [] # versioned building blocks


def drc_angle(cellname, rules, flip_state, xya, flip):
    """Perform instantiation angle drc on a cell.

    This routine will add DRC error to a log file and to a global list to
    be able to export DRC positives to a gds file later.

    Args:
        cellname (str):
        rules (dict):
        flip_state (bool):
        xya (float, float, float):
        flip (bool):

    Returns:
        None
    """
    global cnt, items, drclog
    x, y, a = xya
    msg = ''
    if flip_state is None:
        flip_txt       = ' '
        flip_state_txt = ''
    else:
        flip_txt       = ' AND flip={} '.format(flip)
        flip_state_txt = 'AND flip={}'.format(flip_state)

    valuesOK = True
    domainsOK = True

    values = rules.get('values', None)
    if values is not None:
        values0 = values[:] # copy
        valuesOK = False
        if 0 in values0:
            values0.append(360)
        for value in values0:
            if m.isclose(a, value, rel_tol=1e-9, abs_tol=1e-6):
                valuesOK = True
                break
        if not valuesOK:
            msg += "      angle {:4f}{}not in allowed values: {} {}\n".\
                format(a, flip_txt, values, flip_state_txt)

    domains = rules.get('domains', None)
    if domains is not None:
        domainsOK = False
        for domain in domains:
            if (round(domain[0], 5) <= a) and (a <= round(domain[1], 5)):
                domainsOK = True
                msg = ''
                break
        if not domainsOK:
            msg += "      angle {:4f}{}not in allowed domains: {} {}\n".\
                format(a, flip_txt, domains, flip_state_txt)

    if (not domainsOK and values is None)\
            or (not valuesOK and domains is None)\
            or (not domainsOK and not valuesOK):
        cnt += 1
        items.append((cnt, cellname, xya, flip, msg))
        msg = "#{} cell '{}' @ ({:.3f}, {:.3f}, {:.3f}), flip={}\n".\
            format(str(cnt).zfill(4), cellname, x, y, a, flip) + msg
        #drclog.write(msg)
        for line in msg.split('\n'):
            if len(line) > 0:
                nd.logger.error(line)

def angle_drc(cell, rules=None, basename=None, version_layer=None):
    """Apply intantiation DRC rules to cell hierarchy and scan for black boxes.

    The format of angle DRC rules:

    rule_dict = {
      'angle':
        <cell_basename>:
          'noflip': # optional level if flip state has to be False
            'values': <list of allowed angles>
            'domain': <list allowed domains>
          'flip': # optional level if flip state has to be True
            'values': <list of allowed angles>
            'domain': <list allowed domains>111
      }

    Example::

    somerules = {
      'laser':
        'angle':
          'values': [0, 180]
      'wild_bb':
        'noflip':
          'domains: [[25, 50], [140, 260]]
        'flip':
          'values': [0]
          'domains': [[150, 210]]
      }
    instance_drc(rules=somerules, cell=<your Cell>)

    Args:
        rules (dict): drc rules for instantiation
        cell (Cell): cell (and subcells) to rundrc on
        version_layer (layer): layer containing the cell's version annotation

    Returns:
        Cell: cell with angle DRC result to overlay with the gds design.
    """
    global cnt, drclog, items

    if basename is None:
        basename = "{}".format(cell.cell_name)

    if rules is None:
        rules = nd.cfg.drc_instance_angle

    if version_layer is None:
        nd.logger.warning('version layer missing. Will not create a black box overview.')
    else:
        version_layer = nd.get_layer(version_layer)
    cnt = 0

    #drclogname = os.path.join(basename+'.'+drc_filename)
    #drclog = open(drclogname, 'w')

    #drclog.write("DRC on cell: {}\n".format(cell.cell_name))
    now = datetime.datetime.now()
    #drclog.write('datetime: {}\n'.format(now.strftime("%Y-%m-%d %H:%M")))

    # Get instance info from topcell down.
    positions = []
    for params in nd.cell_iter(cell, flat=True, infolevel=0):
        if params.cell_start:

            # 1. detect black box (besides DRC check):
            version = ''
            for anno, pos in params.iters['annotation']:
                if anno.layer == version_layer:
                    version = anno.text # TODO: what if multiple annotations exist?

            # 2. store position for DRC and possible version:
            positions.append((params.cell.cell_basename, params.cell.cell_name, params.transflip_glob, version))

    # cellname order needs to be reverse-ordered match the longest possible base name first:
    cell_rules_angle =  OrderedDict(sorted(rules['angle'].items(), reverse=True))
    for cellbasename, name, (trans, flip), version in positions:
        x, y, a = trans.xya()
        if version != '':
            blocks.append((cellbasename, version))
        for bb, rules in cell_rules_angle.items():
            if name.startswith(bb):
                angle_rules_flip = rules.get('flip', None)
                angle_rules_noflip = rules.get('noflip', None)
                if angle_rules_flip is not None and flip:
                    drc_angle(name, rules=angle_rules_flip, flip_state=True, xya=[x, y, a], flip=flip)

                elif angle_rules_noflip is not None and not flip:
                    drc_angle(name, rules=angle_rules_noflip, flip_state=False, xya=[x, y, a], flip=flip)

                elif flip and angle_rules_flip is None and angle_rules_noflip is not None:
                    cnt += 1
                    msg = "  flip=True not allowed\n"
                    items.append((cnt, name, [x, y, a], flip, msg))
                    msg = "#{} ERROR cell '{}' @ ({:.3f}, {:.3f}, {:.3f}), flip={}\n{}".\
                        format(str(cnt).zfill(4), name, x, y, a, flip, msg)
                    #drclog.write(msg)
                    nd.logger.error(msg)

                elif not flip and angle_rules_noflip is None and angle_rules_flip is not None:
                    cnt += 1
                    msg = "  flip=False not allowed\n"
                    items.append((cnt, name, [x, y, a], flip, msg))
                    msg = "#{} ERROR cell '{}' @ ({:.3f}, {:.3f}, {:.3f}), flip={}\n{}".\
                        format(str(cnt).zfill(4), name, x, y, a, flip, msg)
                    #drclog.write(msg)
                    nd.logger.error(msg)

                elif rules != {}:
                    drc_angle(name, rules=rules, flip_state=None, xya=[x, y, a], flip=flip)
                break # only the first hit (longest match of a cell in reversed ordered bb is the match
    nd.logger.info("angle violation items found: {}".format(cnt))
    #drclog.write("items found: {}".format(cnt))
    #drclog.close()
    #nd.logger.info("exported {}".format(drclogname))

    # Save DRC results as cells to gds to load on top of the design
    with nd.Cell('{}_drc'.format(cell.cell_name)) as DRC:
        for i in items:
            with nd.Cell(str(i[0]).zfill(4)) as C:
                bbox = nd.cfg.cellnames[i[1]].bbox
                points = [(bbox[0], bbox[1]), (bbox[0], bbox[3]), (bbox[2], bbox[3]), (bbox[2], bbox[1])]
                #print(i[1], i[2], flip, points)
                #print(i)
                nd.Polygon(points=points, layer=nd.cfg.drc_layer_instance_angle).put(0) # put(*i[2], flop=i[3])
            C.put(*i[2], flip=i[3])
    gdsname = os.path.join(basename+'.drc.gds')
    nd.export_gds(DRC, filename=gdsname, clear=False)

    # Save manifest
    if version_layer is not None:
        countBBs = Counter(blocks)
        bblistname = os.path.join(basename+ '.'+blocklist_filename)
        F = open(bblistname, 'w')
        F.write('#datetime: {}\n'.format(now.strftime("%Y-%m-%d %H:%M")))
        L, D, T = nd.cfg.layername2LDT[version_layer]
        F.write("#version_layer: ({}, {})\n".format(L, D))
        F.write("{:6} block_info\n".format('#units'))
        for bb in sorted(countBBs.keys()):
            F.write("{:6} {}\n".format(countBBs[bb], bb))
        F.close()
        nd.logger.info("exported {}".format(bblistname))

    return DRC

if __name__ == '__main__':
    # example on setting design rules on instantiation

    cell_rules = {
        'angle': {
            'mmi1x2_dp': { # cell-basename
                'values': [0, 90, 180, 270],
                'domains': None},
            'die2': {
                'values': None,
                'domains': [[0, 90]]},
            'soa': {
                'noflip': # add level if flip state is to be included
                    {'values': [0, 180],
                     'domains': [[90, 125]]},
                'flip': {'values': [0.1]}
                },
            'soa_sh': {'values': [0, 180]}
            }
        }


    # create/load a design to test:
    gds = '../tests/demofab-masks/demo_example6.gds'
    mask = nd.load_gds(gds)

    instance_drc(cell_rules, mask)
