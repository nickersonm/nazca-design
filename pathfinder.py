#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pathfinder functionality to trace a netlist through a layout.

@authors: Ronald Broeke, based on summer project by Stefan van Ieperen.

Copyright (c) 2019-2020 Bright Photonics B.V.
"""

import sys
import nazca as nd
from nazca import cfg

count = 0

def nb_filter_geo(start):
    """Yield all neighbouring nodes except for the origin, stubs and annotations.

    Iterates over the neighbours in the Cell.

    Args:
        start (Node): a Node from which you want the neighbours

    Yields:
        Node: iterator over start's filtered neighbours
    """
    for start_nb, _ in start.nb_geo:
        if not "org" in start_nb.name:
            yield start_nb


flip = 1
def find_optical_level(start, instance_stack=None, log=False, flip=1, trackertype='dis'):
    """Drill down from a node in the cell hierarchy until an edge is found or max depth is reached.

    This function finds if node <start> has an optical edge.
    If not, it is checked if node is connected to an instance to, if so,
    drill down into it.
    When drilling down this function creates a "C2I_map" attribute (cell2instance)
    in every Instance it passes to find the way back up from a Cell to the Instance along
    the celltree path it came from.

    Args:
        start (Node): a higher level Node from which you want to go down
        log (bool): if True print all nodes passed with an indent proportional to the cell level depth.

    Returns:
        Node, Node, list, float:
            cell_node (Node): start Node of an optical connection,
            nb_opt (Node): end node of optical connection,
            instance_stack (list): list of Instances traversed,
            length (float): length (value) of the connection
    """
    global count, F

    # keep track off / calculate absolute pin coordinates here
    #    with respect to the highest cell (of the start node)

    if instance_stack is None:
        instance_stack = []

    if start.cnode.instance is not None: # instance pin start
        if log:
            F.write(f"{count*'  '}pin '{start.name:10}' in inst '{start.cnode.instance.name:20}' I-{start.cnode.instance.id}\n")
        current_pin = start
        instance = start.cnode.instance
        instance_stack.append(instance)
        if instance.cnode.flip:
            flip *= -1
    else: # cell pin start
        current_pin = start

    #print(f"* start_pos: {start.fxya()}")

    # drill down until an optical edge is found or the bottom is reached:
    while True:
        if current_pin.cnode.instance is not None: # is an instance
            instance.C2I_map = {node.up:node for node in instance.pin.values()}
            cell_node = current_pin.up # go to cell node before drilling down
            if log:
                F.write(f"{count*'  '}pin '{cell_node.name:10}' in cell '{current_pin.cnode.cell.cell_name:20}' C-{current_pin.cnode.cell.id}\n")
        else:
            cell_node = current_pin

        # find all optical paths of cell_node, if any:
        next_nodes_opt = []
       # print(f"* cell_node:   {cell_node.xya()}")
        for nb, L, direction, sigtype, path in cell_node.path_nb_iter(sigtype=trackertype):
            # find geo position of optical nodes for visualisation
            if nb is None:
                x, y, a = 0, 0, 0 # termination
            else:
                x, y, a = nd.diff(cell_node, nb)
            posrel = nd.Pointer(x, y, a)
            if flip == -1:
                posrel.flip()
        #    print(f"  opt_position {posrel}")
            next_nodes_opt.append((nb, L, direction, posrel.xya(), path, sigtype))
        if len(next_nodes_opt) > 0:
            count += 1
            break

        # if no optical edge was found, drill down, if possible:
        for nb, trans in cell_node.nb_geo:
            if nb.cnode.instance is not None and not nb.cnode.cell.auxiliary:
                break # found instance level below
                # assume max *one* instance below can be connected
        if nb.cnode.instance is None:
            break # bottom reached
        else:
            count += 1
            current_pin = nb
            instance = nb.cnode.instance
            if instance.cnode.flip:
                flip *= -1

            # find transloc of
            instance_stack.append(instance)
            if log:
                F.write(f"{count*'  '}pin '{current_pin.name:10}' in inst '{current_pin.cnode.instance.name:20}' I-{current_pin.cnode.instance.id}\n")

    if len(next_nodes_opt) == 0:
       msg = "NO EDGE FOUND WHILE DRILLING DOWN: DEAD END. No path connection found in node {}.\n".format(cell_node)
       F.write(msg)
       #raise Exception(msg)
       return [None] * 4
    else:
       return cell_node, next_nodes_opt, instance_stack, flip


path = 0
def _pathfinder(start, end=None, log=False, logfilename='trace.log',trackertype='dis'):
    """Trace an optical connection by decending an ascending through cells with optical neighbors.

    If the start node is a ribbon pin (A0 or B0), it will trace all paths in
    a ribbon in an arrayed style. Otherwise it traces a singular path.
    opt_netlist holds the pin connection in a default dictionary. In principal,
    this dictionary can hold any data shared between the two pins.

    Args:
        start (Node): starting Node from where it will start the trace
        end (Node): if provided, the trace will stop et the ending node
        log (bool): if True, print all nodes it passes in an indented way
        trackertype (str): Type of connection to trace at the beginning. Default is 'dis' (distance)

    Returns:
        dict: a default dictionary with Pin connections and their lengths
    """
    def foo_print_instances(instance_stack):
        print(', '.join([i.name for i in instance_stack]))

    global path, count, F

    F = open(logfilename, 'w')
    # Ribbon: Trace all ribbon waveguides if start is an array pin
    if start.type is not None:  # obsolete
        # check if pin is an array pin
        if start.type in [3, 4]:
            rib_netlist = {}
            prefix = "a" if start.type == 3 else "b"
            # Iterate over all of the ribbon waveguides to trace them
            for i in range(N):
                pin = start.cnode.instance.pin[f"{prefix}{i}"]
                rib_netlist[pin.name, pin.cnode.instance.name] = pathfinder(pin, end=end, log=log)
            return rib_netlist

    opt_netlist = [] # TODO: describe
    paths = {}
    paths_endpin = {}
    count = 0
    loose_ends = []  # store settings for all neighbours of PIN1 in a tuple before travelling onto the next pin.
    visited = set()
    flip = 1
    start_mem = start

    # single path:
    instance_stack = None
    point_stack = [start.copy()]
    #print('\n')

    end = False
    backtrack = False
    while not end:
        if not instance_stack and backtrack:
            if log:
                F.write(f"{(count-1)*'  '}pin '{start.name:10}' in cell '{start.cnode.cell.cell_name:20}' C-{start.cnode.cell.id}\n")

        # Drill down from 'start' to a PIN1 with an optical edge between PIN1 and PIN2 at cell level.
        if not backtrack:
            PIN1, PIN2list, instance_stack, flip = find_optical_level(start, instance_stack, log, flip, trackertype)

            if PIN1 is None:
                #if log:
                    #F.write(f"{(count+1)*'  '}=== END PATH {path}: LOOP CLOSURE IN C-{PIN1.cnode.cell.id} on pin '{PIN1.name}'\n\n")
                if start_mem == start:
                    nd.logger.warning(f"Desired start pin {start.cnode.cell.cell_name} has no desired trackertype ({trackertype}).")
                else:
                    paths[path] = opt_netlist
                    paths_endpin[path] = pinid
                    path += 1
                if len(loose_ends) > 0:
                    backtrack = True
                    if log:
                        F.write(f"{(count)*'  '}=== BACK TRACK ===\n")
                    continue
                else:
                    end=True
                    break

            pinid = (tuple(instance_stack), PIN1, trackertype)
            if pinid in visited:
                if log:
                    F.write(f"{(count+1)*'  '}=== END PATH {path}: LOOP CLOSURE IN C-{PIN1.cnode.cell.id} on pin '{PIN1.name}'\n\n")
                paths[path] = opt_netlist
                paths_endpin[path] = pinid
                path += 1
                if len(loose_ends) > 0:
                   backtrack = True
                   if log:
                       F.write(f"{(count)*'  '}=== BACK TRACK ===\n")
                else:
                   end = True
                continue

            visited.add(pinid)
            PIN2tup = PIN2list.pop()
            for tup in PIN2list:
                loose_ends.append((
                    PIN1,
                    tup,
                    instance_stack.copy(),
                    count,
                    opt_netlist.copy(),
                    point_stack.copy(),
                    visited.copy(),
                    flip)
                )
        else:
            (   PIN1,
                PIN2tup,
                instance_stack,
                count,
                opt_netlist,
                point_stack,
                visited,
                flip
            ) = loose_ends.pop()
            if log:
                F.write(f"{(count-1)*'  '}pin '{PIN1.name:10}' in cell '{PIN1.cnode.cell.cell_name:20}' C-{PIN1.cnode.cell.id}\n")

        #PIN2, L, dir, point, line = PIN2tup
        old_trackertype = trackertype
        PIN2, L, dir, point, line, trackertype = PIN2tup
        backtrack = False
        x, y, a = point
        point_old = point_stack[-1]
        point_new = point_old.move(x, y, a+180)
        point_stack.append(point_new)

        if L is None:
            L = 0
            nd.logger.warning(f"Optical connect not None in cell {PIN2.cnode.cell.cell_name}. Setting to 0.")
            # TODO: raise a netlist error?
        opt_netlist.append((PIN1, PIN2, L, (point_old, point_new), line, (old_trackertype, trackertype)))

        if log:
            F.write(f"{(count)*'  '}=== JUMP OPTICAL lINK, length={L:0.3f} ===\n")
            #F.write(f"{(count-1)*'  '}pin '{PIN2.name:10}' in cell '{PIN2.cnode.cell.cell_name:20}' C-{PIN2.cnode.cell.id}\n")


        if start.cnode.instance is None and not instance_stack:
            if log:
                F.write(f"{count*'  '}=== BREAK OUT TOPCELL ===\n")
            paths[path] = opt_netlist
            paths_endpin[path] = (tuple(instance_stack), PIN2)
            break


        # Drill upward from cell PIN2 (with optical connection) until a connection to a next instance is found:
        while instance_stack:
            if log:
                F.write(f"{(count-1)*'  '}pin '{PIN2.name:10}' in cell '{PIN2.cnode.cell.cell_name:20}' C-{PIN2.cnode.cell.id}\n")

            nbs_side = []
            nbs_up = []
            instance = instance_stack[-1]
            try:
                pin2 = instance.C2I_map[PIN2] # instance pin2 for cell PIN2 up tree.
            except:
                if log:
                    F.write(f"{count*'  '}END PATH {path}: No way back to instance from cell '{PIN2.cnode.cell.cell_name}' for pin '{PIN2.name}'\n")
                paths[path] = opt_netlist
                paths_endpin[path] = (tuple(instance_stack), PIN2, trackertype)
                path += 1
                if len(loose_ends) > 0:
                    backtrack = True
                    if log:
                        F.write(f"{(count)*'  '}=== BACK TRACK ===\n")
                    count -= 1
                else:
                    end = True
                break
            count -= 1
            if log:
                F.write(f"{count*'  '}pin '{pin2.name:10}' in inst '{pin2.cnode.instance.name:20}' I-{pin2.cnode.instance.id}\n")
            # TODO: look for optical connections?
            for node in nb_filter_geo(pin2):
                if (
                    node.cnode.instance is not pin2.cnode.instance
                    and node.cnode.parent_cnode is pin2.cnode.parent_cnode
                    and node.width != 0
                    # check distance == 0
                    # check for begin in a specific xs?
                ):
                    nbs_side.append(node)
                elif pin2.cnode.parent_cnode is node.cnode:
                    nbs_up.append(node)
            if nbs_side:
                start = nbs_side[0] # TODO: iterate over the full set of side-nodes.
                if len(nbs_side) > 1:
                    raise Exception('Splitting Node')
                I = instance_stack.pop()
                if I.cnode.flip:
                    flip *= -1
                if log:
                    F.write(f"{(count)*'  '}=== JUMP TO INSTANCE ===\n")
                break
            elif nbs_up: # no neighbours found, back up: look for nb in parent
                I = instance_stack.pop()
                if I.cnode.flip:
                    flip *= -1

                PIN2 = nbs_up[0]
                if not instance_stack: # top cell
                    if log:
                        F.write(f"{(count-1)*'  '}END PATH {path} IN CELL '{nbs_up[0].cnode.cell.cell_name}' pin '{nbs_up[0].name}'\n\n")
                    paths[path] = opt_netlist
                    paths_endpin[path] = (tuple(instance_stack), PIN2, trackertype)
                    path += 1
                    start = PIN2
                    if len(loose_ends) > 0:
                        backtrack = True
                        if log:
                            F.write(f"{(count-1)*'  '}=== BACK TRACK ===\n")
                    else:
                        end = True

            else: # end of path, go back to loose ends
                paths[path] = opt_netlist
                paths_endpin[path] = (instance_stack, PIN2, trackertype)
                path += 1
                if log:
                    F.write(f"{count*'  '}END PATH {path}: in instance '{pin2.cnode.cell.cell_name}' for pin '{pin2.name}'\n\n")
                if len(loose_ends) > 0:
                   backtrack = True
                   if log:
                       F.write(f"{(count)*'  '}=== BACK TRACK ===\n")
                else:
                   end = True
                break

    F.close()
    #Ltot = 0
    #for p1, p2, L, point, line in opt_netlist:
    #    Ltot += L
    #print(f"Ltot = {Ltot}")
    return paths, paths_endpin

trace_opt = _pathfinder # for backward compatibility of tests


def print_paths(paths, width, layer, stdout=True, endpins=None, pathfilename=None):
    """Print a list of all paths and their lengths.

    Returns:
        None
    """
    handles = []
    if pathfilename is not None:
        handles.append(open(pathfilename, 'w'))
    if stdout:
        handles.append(sys.stdout)
    for F in handles:
        F.write("number of paths: {}:\n".format(len(paths)))
        F.write(f"{'num':3}: {'    length':10}   cellpath/pin\n")
        F.write(f"------------------------------\n")
    for pathnum, stuff in paths.items():
        points = []
        Ltot = 0
        for i, (p1, p2, L, point, line, track_tup) in enumerate(stuff):
            Ltot += L
            x0, y0 = point[0].xya()[:2]
            if i> 0:
                if (x0, y0) != points[-1]: # avoid overlapping point that will have lelngth 0 between them.
                    points.append((x0, y0))
            else:
                points.append((x0, y0))

            if line is not None:
                for (x, y) in line:
                    #print("line", line)
                    points.append((x0+x, y0+y))
        try:
            point
        except NameError:
            nd.logger.warning("path '{}' has no points".format(path))
        else:
            points.append(point[1].xya()[:2])
            if stdout:
                location = ''
                if endpins is not None:
                    tuptree, pin, tracktype = endpins[pathnum]
                    tree = [ins.cell.cell_name for ins in tuptree]
                    location = "{}/{}".format('/'.join(tree), pin.name)
                for F in handles:
                    F.write(f"{pathnum:3}: {Ltot:10.3f}   {location}\n")
            #print(points, '\n')
            nd.Polyline(points, layer=pathnum+layer, width=width, pathtype=1).put(0)
            nd.show_pin(nd.Pin().put(points[0]), radius=2*width, width=1*width, layer=pathnum+layer)
            nd.show_pin(nd.Pin().put(points[-1]), radius=1.5*width, width=1.5*width, layer=pathnum+layer)


def print_optical_paths(paths, width, layer, stdout=True, endpins=None, pathfilename=None):
    """Print a list of all paths and their lengths.

    Returns:
        None
    """
    handles = []
    if pathfilename is not None:
        handles.append(open(pathfilename, 'w'))
    if stdout:
        handles.append(sys.stdout)
    for F in handles:
        F.write("number of paths: {}:\n".format(len(paths)))
        F.write(f"{'num':3}: {'    length':10}   cellpath/pin\n")
        F.write(f"------------------------------\n")
    for pathnum, stuff in paths.items():
        polylines=[]
        points = []
        Ltot = 0
        for i, (p1, p2, L, point, line, track_tup) in enumerate(stuff):
            Ltot += L
            x0, y0 = point[0].xya()[:2]
            if i> 0:
                if (x0, y0) != points[-1]: # avoid overlapping point that will have lelngth 0 between them.
                    points.append((x0, y0))
            else:
                points.append((x0, y0))

            if line is not None:
                for (x, y) in line:
                    #print("line", line)
                    points.append((x0+x, y0+y))
            if track_tup[0]!=track_tup[1]:
                polylines.append((points,int(track_tup[0][-1])))
                points=[points[-1]]
        try:
            point
        except NameError:
            nd.logger.warning("path '{}' has no points".format(path))
        else:
            points.append(point[1].xya()[:2])
            polylines.append((points,int(track_tup[1][-1])))
            if stdout:
                location = ''
                if endpins is not None:
                    tuptree, pin, tracktype = endpins[pathnum]
                    tree = [ins.cell.cell_name for ins in tuptree]
                    location = "{}/{}".format('/'.join(tree), pin.name)
                for F in handles:
                    F.write(f"{pathnum:3}: {Ltot:10.3f}   {location}\n")
            #print(points, '\n')
            for poly,datatype in polylines:
                nd.Polyline(poly, layer=(pathnum+layer,datatype), width=width, pathtype=1).put(0)
                nd.show_pin(nd.Pin().put(poly[0]), radius=2*width, width=1*width, layer=(pathnum+layer,datatype))
                nd.show_pin(nd.Pin().put(poly[-1]), radius=1.5*width, width=1.5*width, layer=(pathnum+layer,datatype))


def pathfinder(start,
        end=None,
        stdout=True,
        pathfilename=None,
        log=False,
        logfilename='trace.log',
        width=2.0,
        layer=2000,
        trackertype='dis'
        ):
    """Find and visualize circuit paths.

    The paths contain the geometrical length between two pins.

    Args:
        start (Node): start Node to trace from
        end (Node): not in use
        stdout (bool): write path information to stdout (default = True)
        pathfilename (str): optional filename for saving pathfinder results
            (default=None)
        log (bool): generate log file of complete trace
        logfilename (str): filename of logfile
        width (float); width of the polyline visualizing the path in gds
        layer (int): start gds layer number to visualize the paths in gds.
            Path polylines will be placed in sequential layer numbers
        trackertype (str): Type of connection to trace at the beginning.
            Default is 'dis' (distance)

    Returns:
        None
    """
    print("start pin: {}.{}".format(start.cnode.cell.cell_name, start.name))
    paths, paths_endpin = _pathfinder(start=start, end=end, log=log,
        logfilename=logfilename,trackertype=trackertype)
    #print_endpins(paths_endpin)
    print_paths(paths, width=width, layer=layer, endpins=paths_endpin,
        pathfilename=pathfilename, stdout=stdout)


def path2lyp(start,
        end=None,
        stdout=True,
        pathfilename=None,
        log=False,
        logfilename='trace.log',
        width=2.0,
        layer=2000,
        trackertype='dis',
        separate_pol=False,
        group_colors=False,
        ):

        if isinstance(start, list):
            pass
        else:
            start=[start]

        color_dic={}

        if group_colors and separate_pol:
            nd.logger.warning(f"From path2Lyp: group_colors=True not available with separate_pol=True")

        paths={}
        paths_endpin={}
        for j,start_pin in enumerate(start):
            new_paths, new_paths_endpin = _pathfinder(start_pin, end=end, log=log, logfilename=f'{logfilename}.{j}',trackertype=trackertype)
            paths.update(new_paths)
            paths_endpin.update(new_paths_endpin)

            for i,end in new_paths_endpin.items():
                try:
                    L=[stuff[2] for stuff in new_paths[i]]
                    Ltot=sum(L)
                    path_id=str((start_pin.remark,trackertype,end[0][-1].pin[end[1].name].remark,end[2],f'{Ltot:.3f}'))
                    if stdout:
                        print(f'{i} : {str(start_pin.remark):20s} ({trackertype:4s}) --> {str(end[0][-1].pin[end[1].name].remark):20s} ({end[2]:4s}) : {Ltot:10.2f}')
                    if separate_pol:
                        max_pol=max([int(x[-1]) for stuff in new_paths[i] for x in stuff[5] ])
                        for j in range(max_pol+1):
                            nd.add_layer(name=f'Path{i:03} pol{j} : {str(start_pin.remark):20s} ({trackertype:4s}) --> {str(end[0][-1].pin[end[1].name].remark):20s} ({end[2]:4s}) : {Ltot:10.2f}',layer=(layer+i,j),dither_pattern=f'I{j}')
                    else:
                        lay_name=f'Path{i:03} : {str(start_pin.remark):20s} ({trackertype:4s}) --> {str(end[0][-1].pin[end[1].name].remark):20s} ({end[2]:4s}) : {Ltot:10.2f}'
                        if path_id in color_dic:
                            nd.add_layer(name=lay_name,layer=(layer+i,0),fill_color=color_dic[path_id],frame_color=color_dic[path_id])
                        else:
                            nd.add_layer(name=lay_name,layer=(layer+i,0))
                            if group_colors:
                                color_dic[path_id]=cfg.colors.loc[cfg.colors["name"].str.contains(f'Path{i:03}')]['fill_color'].values[0]
                except KeyError:
                    print(i,end)
                    nd.add_layer(name='Path%03i' % (i),layer=(3000+i,0))


        path_colors=cfg.colors[cfg.colors["name"].str.contains("Path")]
        #path_colors['width']=3
        path_colors.loc[:,'width']=3
        nd.nazca2csv(path_colors, 'path_colors.csv')
        nd.csv2lyp({'Paths': 'path_colors.csv'}, 'path_colors.lyp')
        if separate_pol:
            print_optical_paths(paths, width, layer, stdout=False, endpins=paths_endpin, pathfilename=pathfilename)
        else:
            print_paths(paths, width, layer, stdout=False, endpins=paths_endpin, pathfilename=pathfilename)

        return paths,paths_endpin



if __name__ == '__main__':
    from nazca.interconnects import Interconnect
    ic = Interconnect(radius=10)

    def splitter_1x2():
        """Mimmic a 1x2 power splitter for testing."""
        with nd.Cell('splitter', autobbox=True) as C:
            p1 = nd.Pin('a0').put(0, 0, 180)
            p2 = nd.Pin('b0').put(30, 10, 0)
            p3 = nd.Pin('b1').put(15, -4, 0)
            sb1 = ic.sbend_p2p(p1.rot(180), p2.rot(180)).put()
            sb2 = ic.sbend_p2p(p1.rot(180), p3.rot(180)).put()
            nd.connect_path(p1, p2, sb1.length_geo)
            nd.connect_path(p1, p3, sb2.length_geo)
            nd.put_stub()
        return C

    def splitter_2x2():
        """Mimmic a 2x2 power splitter for testing."""
        with nd.Cell('splitter2', autobbox=True) as C:
            p1 = nd.Pin('a0').put(0, 10, 180)
            p2 = nd.Pin('a1').put(0, -6, 180)
            p3 = nd.Pin('b0').put(30, 10, 0)
            p4 = nd.Pin('b1').put(30, -4, 0)
            sb1 = ic.sbend_p2p(p1.rot(180), p3.rot(180)).put()
            sb2 = ic.sbend_p2p(p1.rot(180), p4.rot(180)).put()
            sb3 = ic.sbend_p2p(p2.rot(180), p3.rot(180)).put()
            sb4 = ic.sbend_p2p(p2.rot(180), p4.rot(180)).put()
            nd.connect_path(p1, p3, sb1.length_geo)
            nd.connect_path(p1, p4, sb2.length_geo)
            nd.connect_path(p2, p3, sb3.length_geo)
            nd.connect_path(p2, p4, sb4.length_geo)
            nd.connect_path(p1, p2, 10)
            nd.connect_path(p3, p4, 8)
            nd.put_stub()
        return C


    with nd.Cell('TEST') as C:
        b1 = nd.bend().put(flip=True)
        sp0 = splitter_1x2().put()
        sp1 = splitter_1x2().put('b1', flip=True)
        sp2 = splitter_1x2().put(sp1.pin['a0'])
        sp3 = splitter_2x2().put()
        splitter_1x2().put()
        sp5 = splitter_2x2().put()
        ic.bend_strt_bend_p2p(sp5.pin['a1'], sp5.pin['b1'], ictype='ll').put()

        nd.Pin('a0', pin=b1.pin['a0']).put()
        print("\nPath from cell:")
        pathfinder(start=C.pin['a0'], log=True, pathfilename='paths.log',
            logfilename='trace_cell.log')


    c = C.put(0, 0, 90)
    path = 0
    print("\nPath from instance:")
    pathfinder(start=c.pin['a0'], log=True, logfilename='trace_inst.log', layer=3000)

    nd.export_gds(clear=False)
