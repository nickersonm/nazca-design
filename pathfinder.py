#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pathfinder functionality to trace a netlist through a layout.

@authors: Jerom Baas, Ronald Broeke, based on summer project by Stefan van Ieperen.

Copyright (c) 2019-2022 Bright Photonics B.V.
"""
from typing import Any, List, Optional, Dict
from copy import copy
from dataclasses import dataclass

# TODO look into making find_optical connections more readable

import nazca as nd
from nazca.util import ProtectedPartial, filter_eval


@dataclass
class Segment:
    """Class for storing the info a single path segment"""

    pin_in: nd.Node
    pin_out: nd.Node
    connection: Any
    position_in: Any
    position_out: Any
    upper_pin_in: Optional[nd.Node] = None
    upper_pin_out: Optional[nd.Node] = None
    tracker_kwargs_in: Optional[Dict] = None
    tracker_kwargs_out: Optional[Dict] = None


class Logger:
    """
    Class for logging outputs of the pathfinder.

    Args:
        filename (str): Name of the logfile to write to.
        logging_level (int): Level of logging. At 0 no logfile is created. At 1 a logfile is created with different
                            pins traversed. At 2, the act of moving up and down the hierarchy is also logged.
        depth (int): Attribute used to manage indentations in the logging file that correspond to hierarchical depth.
    Methods:
        log(): Writes a different log message to the logfile depending on the status string provided.
    """

    def __init__(
        self,
        filename="trace.log",
        logging_level=0,
        depth=0,
        trackertype='dis',
    ):
        self.filename = filename
        self.logging_level = logging_level
        self.depth = depth
        self.trackertype = trackertype

        if self.logging_level > 0:
            self.logfile = open(filename, "w")
        elif self.logging_level == -1:
            self.logfile = open(filename, "w")

    def log(
        self,
        pin,
        status,
        L=0,
        instance_stack=None,
        endpoint=None,
        connected_pins=None,
        paths_dict=None,
    ):
        """
        Outputs message to logging file if logging_level > 0. Outputs more detailed message if logging level > 1.
        The message depends on the status passed to the log function.
        Possible statuses:
            logging level 1:
                'start'
                'optical_connection'
                'pin_connection'
                'loop'
                'end_of_path'
                'end_of_path_splitting_node'
                'endpoint'
                'backtracking'
                'backtracked'
                'end_of_tracking'
            logging level 2:
                'drill_down'
                'drill_up'
                'no_pin_connection'
                'no_optical_connection'
                'reverse'
        Args:
            pin (node): current pin at the time of logging.
            status (str): string that determines what kind of message is logged.
            L (float): length of optical connection found. Used only if status is 'optical_connection'.
            instance_stack (list): Instance stack of instances traversed. Used for hierarchy in backtracking.
            endpoint (nd.Pin, nd.Cell or nd.Instance): Used for status 'endpoint', to specify which specific
                                                        user-provided endpoint was reached.
            connected_pins (list): List of connected pins. This is used to output the list of connections when status is
                                   'end_of_path_splitting_node' so that debugging is easier.
            paths_dict (dict): dictionary of final paths
        """
        if self.logging_level > 0:
            if status == "start":
                self.depth = 0
                self.logfile.write(f"Start of tracking paths...\n")
                self.depth = 1
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth * '  '}Starting path tracking at pin '{pin.name}' in instance of "
                        f"cell '{pin.cnode.cell.cell_name}'.\n"
                    )
                except:
                    self.logfile.write(
                        f"{self.depth*'  '}Starting path tracking at pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n"
                    )
            elif status == "optical_connection":
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth * '  '}Found optical connection to pin '{pin.name}' in instance "
                        f"of cell '{pin.cnode.cell.cell_name}' with length {L}.\n"
                    )
                    if self.logging_level > 1:
                        self.logfile.write(
                            f"{self.depth * '  '}------------------------------------------------------"
                            f"-----------------\n"
                        )
                        self.logfile.write(
                            f"{self.depth * '  '}Looking for pin connections to pin '{pin.name}' "
                            f"in instance of cell '{pin.cnode.cell.cell_name}'...\n"
                        )
                except:
                    self.logfile.write(
                        f"{self.depth*'  '}Found optical connection to pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}' with length {L}.\n"
                    )
                    if self.logging_level > 1:
                        self.logfile.write(
                            f"{self.depth * '  '}------------------------------------------------------"
                            f"-----------------\n"
                        )
                        self.logfile.write(
                            f"{self.depth * '  '}Looking for pin connections to pin '{pin.name}' in "
                            f"cell '{pin.cnode.cell.cell_name}'...\n"
                        )
            elif status == "pin_connection":
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth*'  '}Found pin-pin connection to pin '{pin.name}' in instance of "
                        f"cell '{pin.cnode.cell.cell_name}'.\n"
                    )
                    if self.logging_level > 1:
                        self.logfile.write(
                            f"{self.depth * '  '}------------------------------------------------------"
                            f"-----------------\n"
                        )
                        self.logfile.write(
                            f"{self.depth * '  '}Looking for optical connections to pin '{pin.name}' "
                            f"in instance of cell '{pin.cnode.cell.cell_name}'...\n"
                        )
                except:
                    self.logfile.write(
                        f"{self.depth * '  '}Found pin-pin connection to pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n"
                    )
                    if self.logging_level > 1:
                        self.logfile.write(
                            f"{self.depth * '  '}------------------------------------------------------"
                            f"-----------------\n"
                        )
                        self.logfile.write(
                            f"{self.depth * '  '}Looking for optical connections to pin '{pin.name}' "
                            f"in cell '{pin.cnode.cell.cell_name}'...\n"
                        )
            elif status == "loop":
                self.depth = 1
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth * '  '}Detected loop at pin '{pin.name}' in instance of cell "
                        f"'{pin.cnode.cell.cell_name}'.\n\n"
                    )
                except:
                    self.logfile.write(
                        f"{self.depth*'  '}Detected loop at pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n\n"
                    )
            elif status == "end_of_path":
                self.depth = 1
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth * '  '}Reached end of path at pin '{pin.name}' in instance of "
                        f"cell '{pin.cnode.cell.cell_name}'.\n\n"
                    )
                except:
                    self.logfile.write(
                        f"{self.depth*'  '}Reached end of path at pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n\n"
                    )
            elif status == "end_of_path_splitting_node":
                self.depth = 1
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth * '  '}Splitting nodes: connections to:\n"
                    )
                    for connection in connected_pins:
                        self.logfile.write(
                            f"{self.depth * '  '}\t- {connection}\n"
                        )
                    self.logfile.write(
                        f"{self.depth * '  '}Terminated path at pin '{pin.name}' in instance of "
                        f"cell '{pin.cnode.cell.cell_name}'.\n\n"
                    )
                except:
                    self.logfile.write(
                        f"{self.depth * '  '}Splitting nodes: connections to:\n"
                    )
                    for connection in connected_pins:
                        self.logfile.write(
                            f"{self.depth * '  '}\t- {connection}\n"
                        )
                    self.logfile.write(
                        f"{self.depth*'  '}Terminated path at pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n\n"
                    )
            elif status == "endpoint":
                self.depth = 1
                # figure out if endpoint is pin, cell or instance
                if isinstance(endpoint, nd.Node):
                    self.logfile.write(
                        f"{self.depth * '  '}Reached endpoint: pin is endpoint pin.\n"
                    )
                    try:
                        # check if pin has an instance
                        endpoint.cnode.instance.name
                        # if pin has an instance, slightly change output message
                        self.logfile.write(
                            f"{self.depth * '  '}Reached end of path at pin '{endpoint.name}' "
                            f"in instance of cell '{endpoint.cnode.cell.cell_name}'.\n\n"
                        )
                    except:
                        self.logfile.write(
                            f"{self.depth * '  '}Reached end of path at pin '{endpoint.name}' in cell "
                            f"'{endpoint.cnode.cell.cell_name}'.\n\n"
                        )
                elif isinstance(endpoint, nd.Cell):
                    self.logfile.write(
                        f"{self.depth * '  '}Reached endpoint: pin "
                        f"in cell '{endpoint.cell_name}'.\n"
                    )
                    try:
                        # check if pin has an instance
                        pin.cnode.instance.name
                        # if pin has an instance, slightly change output message
                        self.logfile.write(
                            f"{self.depth * '  '}Reached end of path at pin '{pin.name}' in instance of "
                            f"cell '{pin.cnode.cell.cell_name}'.\n\n"
                        )
                    except:
                        self.logfile.write(
                            f"{self.depth * '  '}Reached end of path at pin '{pin.name}' in cell "
                            f"'{pin.cnode.cell.cell_name}'.\n\n"
                        )
                elif isinstance(endpoint, nd.Instance):
                    self.logfile.write(
                        f"{self.depth * '  '}Reached endpoint: pin in {endpoint}.\n"
                    )
                    try:
                        # check if pin has an instance
                        pin.cnode.instance.name
                        # if pin has an instance, slightly change output message
                        self.logfile.write(
                            f"{self.depth * '  '}Reached end of path at pin '{pin.name}' in instance of "
                            f"cell '{pin.cnode.cell.cell_name}'.\n\n"
                        )
                    except:
                        self.logfile.write(
                            f"{self.depth * '  '}Reached end of path at pin '{pin.name}' in cell "
                            f"'{pin.cnode.cell.cell_name}'.\n\n"
                        )
            elif status == "backtracking":
                self.depth = 0
                self.logfile.write(f"{self.depth*'  '}Backtracking...\n")
            elif status == "backtracked":
                if self.logging_level > 1:
                    self.depth = len(instance_stack) + 1
                else:
                    self.depth = 1
                self.logfile.write(
                    f"{self.depth*'  '}Backtracked to pin '{pin.name}' in cell "
                    f"'{pin.cnode.cell.cell_name}'.\n"
                )
            elif status == "end_of_tracking":
                self.depth = 0
                self.logfile.write(
                    f"\n\nTracking done...\n\nFound {len(pin)} paths with the following endpins:\n"
                )
                for i, p in enumerate(pin):
                    if len(paths_dict[i]) > 0:
                        try:
                            p.cnode.instance.name
                            self.logfile.write(
                                f"{i+1}. Pin {p.name} in instance of {p.cnode.cell.cell_name} at"
                                f" {paths_dict[i][-1].position_out}.\n"
                            )
                        except:
                            self.logfile.write(
                                f"{i + 1}. Pin {p.name} in cell {p.cnode.cell.cell_name} at"
                                f" {paths_dict[i][-1].position_out}.\n"
                            )

        if self.logging_level > 1:
            if status == "drill_up":
                self.depth -= 1
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth*'  '}Move up to pin '{pin.name}' in instance of "
                        f"cell '{pin.cnode.cell.cell_name}'.\n"
                    )
                except:
                    self.logfile.write(
                        f"{self.depth * '  '}Move up to pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n"
                    )
            elif status == "drill_down":
                self.depth += 1
                try:
                    # check if pin has an instance
                    pin.cnode.instance.name
                    # if pin has an instance, slightly change output message
                    self.logfile.write(
                        f"{self.depth * '  '}Move down to pin '{pin.name}' in instance of "
                        f"cell '{pin.cnode.cell.cell_name}'.\n"
                    )
                except:
                    self.logfile.write(
                        f"{self.depth * '  '}Move down to pin '{pin.name}' in cell "
                        f"'{pin.cnode.cell.cell_name}'.\n"
                    )
            elif status == "no_optical_connection":
                self.logfile.write(
                    f"{self.depth * '  '}No optical connection of type '{self.trackertype}' found...\n"
                )
            elif status == "no_pin_connection":
                self.logfile.write(f"{self.depth * '  '}No connected pins found...\n")
            elif status == "reverse":
                self.logfile.write(f"{self.depth * '  '}Reversed tracking. Drilling down first...\n")


class Paths:
    """
    The Paths class creates an object that contains information on the different paths from a starting pin.

    Attributes:
        endpins (list): A list containing the endpins of the different paths.
        segments (dict): A dictionary containing the raw output data from the pathfinder.
        endpins_dict (dict): A dictionary containing the endpins and their positional information.
        number_of_paths (float): The number of different paths found by the pathfinder.
        startpin (node): The starting pin from which the paths were traced.
        lengths (list): A list containing the lengths of the different paths.

    Methods:
        show_paths(): A function that visualizes the traced paths when the layout is exported to a gds.
    """

    def __init__(
        self,
        startpin,
        endpoints=None,
        logfilename="trace.log",
        logging_level=0,
        reverse=False,
        allow_splitting_nodes=False,
        trackertype="dis",
        instance_stack=None,
        tracker_kwargs=None,
        show=False,
        show_layer_number=2000,
        show_width=7,
    ):
        """
        Initializes the paths by finding the endpins and paths from a startpin using the _findpaths() function,
        calculating their lengths and the number of paths.

        Args:
            startpin (Node): Startpoint for the pathfinding algorithm.
            logfilename (str): Name of the logfile.
            endpoints (list): Optional List of pins, cells and/or instances where the pathfinder will terminate the paths
            logging_level (int): Level of logging. At 0 no logfile is created. At 1 a logfile is created with
                    different pins traversed. At 2, the act of moving up and down the hierarchy is also logged.
            reverse (bool): If true, pathfinder will start tracking in the reverse direction. The default direction is
                    to look for optical connections at the start. If reverse is true, the pathfinder will look for
                    pin2pin connections first.
            allow_splitting_nodes (bool): If true, the pathfinder will handle a splitting path in pin2pin connections.
                    By default, this is false: If this happens, the pathfinder will terminate the path and raise an error
            trackertype (str): Type of optical connections to look for.
            instance_stack (list): Optional instance stack that a user can provide. This allows the pathfinder to break
                    out of the starting level.
            show (bool): Show path in the cell. default=False.

        Returns:
            None
        """

        if startpin is None:
            msg = "You must specify a start pin to trace paths from."
            raise RuntimeError(msg)

        if tracker_kwargs is None:
            tracker_kwargs = {}

        # find the paths
        endpins, paths_dict = self._findpaths(
            startpin=startpin,
            endpoints=endpoints,
            logfilename=logfilename,
            logging_level=logging_level,
            reverse=reverse,
            allow_splitting_nodes=allow_splitting_nodes,
            trackertype=trackertype,
            instance_stack=instance_stack,
            tracker_kwargs=tracker_kwargs,
        )

        self.endpins = endpins
        self.segments = paths_dict
        self.endpins_dict = self._generate_endpins_dictionary(paths_dict, startpin, instance_stack)
        #self.lengths = self._get_lengths(self.segments)
        self.number_of_paths = len(self.segments)
        self.startpin = startpin
        self.pathlayers = set()  # keep track of layers added by pathfinder to allow pathfinder to reuse them.
        if show:
            self.show(layer_number=show_layer_number, width=show_width)


    @staticmethod
    def get_pin_coordinate_up_instance_stack(pin, instance_stack):
        """Get the xy coordinates of the pin goin up the instance stack"""
        flipstate = False
        pointer = nd.Pointer(0, 0, 0)
        for instance in instance_stack:
            move = copy(instance.cnode.pointer)
            if flipstate:
                move.flip()
            pointer.move_ptr(move)
            flipstate = not flipstate if instance.cnode.flip else flipstate

        move = copy(pin.pointer)
        if flipstate:
            move.flip()
        pointer.move_ptr(move)

        return pointer.xy()

    def show_paths(
        self,
        layer_number=2000,
        width=5.0,
        paths2cell=None,
    ):
        """
        Visualize paths in the layout via a polyline.

        Args:
            layer_number (int): layer number for path 1.
            width (float): width of the lines drawn between the nodes of a path.
            paths2cell (nd.Cell): Cell to which to add the paths. Default is the active cell.

        Returns:
            None
        """

        if paths2cell is not None:
            nd.cfg.cells.append(paths2cell)
            nd.cfg.patchcell = True

        points_dict = {}
        for i, segments in self.segments.items():
            if len(segments) == 0:
                points_dict[i] = []
                continue
            # for each path, create a new layer for visualization
            layername = f"path{layer_number+i}"
            if layername not in nd.cfg.reuse_pathlayers:
                layer = nd.add_layer(layername, (layer_number + i, 0))
                nd.cfg.reuse_pathlayers.add(layername)
            else:
                layer = nd.get_layer(layername)
            points = []
            for segment in segments:
                points.append(segment.position_in)
            points.append(segments[-1].position_out)

            # place polyline that visualizes the paths
            nd.Polyline(points, layer=layer, width=width, pathtype=1).put(0, 0, 0)

            # place pins that visualize the startpin and endpin of each path.
            nd.show_pin(pin=points[0], radius=2 * width, width=1 * width, layer=layer)
            nd.show_pin(pin=points[-1], radius=1.5 * width, width=1.5 * width, layer=layer)

            # update dictionary of the positions. Can be used for debugging purposes.
            points_dict[i] = points

        if paths2cell is not None:
            nd.cfg.cells.pop()
            nd.cfg.patchcell = False

        return points_dict
    show = show_paths

    def _findpaths(
        self,
        startpin,
        endpoints=None,
        logfilename="trace.log",
        logging_level=0,
        reverse=False,
        allow_splitting_nodes=False,
        trackertype="dis",
        instance_stack=None,
        tracker_kwargs=None,
    ):
        """
        Function to find the paths originating from a certain starting pin. The tracing follows the following algorithm:
        Step 1: Look for optical connections in a cell and move to a connection. For example, if the start pin is pin a0
         of a straight section step 1 will yield the b0 pin of that same straight section and the length from a0 to b0
        Step 2: Look for pin to pin connections from a pin. This step yields pins connected to the pin found in step 1.

        If in step 1 multiple optical connections are found, the state for those connections is stored as loose ends.
        When the end of the current path is reached, the algorithm will backtrack to one of the loose ends and continue
        tracking from there. Tracking continues until all paths have reached an endpoint.

        Args:
            startpin (Node): Startpoint for the pathfinding algorithm.
            logfilename (str): Name of the logfile.
            endpoints (list): Optional List of pins, cells and/or instances where the pathfinder will terminate the paths
            logging_level (int): Level of logging. At 0 no logfile is created. At 1 a logfile is created with
                    different pins traversed. At 2, the act of moving up and down the hierarchy is also logged.
            reverse (bool): If true, pathfinder will start tracking in the reverse direction. The default direction is
                    to look for optical connections at the start. If reverse is true, the pathfinder will look for
                    pin2pin connections first.
            allow_splitting_nodes (bool): If true, the pathfinder will handle a splitting path in pin2pin connections.
                    By default, this is false: If this happens, the pathfinder will terminate the path and raise an error
            trackertype (str): Type of optical connections to look for.
            instance_stack (list): Optional instance stack that a user can provide. This allows the pathfinder to break
                    out of the starting level.
        Returns:
            endpins (list): List of the endpins (Nodes) of all paths that were found.
            paths_dict (dict): Dictionary of all paths that were found. Each path entry contains a list of [pin, length]
                        where pin is the next pin and length is the optical length between that pin and the previous pin
        """

        # initialize logger
        logger = Logger(filename=logfilename, logging_level=logging_level, trackertype=trackertype)

        # define parameters for tracking of paths
        # instance_stack stores hierarchical information, which is used to look for pin connections up in the hierarchy
        if instance_stack is None:
            instance_stack = []
            current_pin = startpin
        else:
            instance_stack = instance_stack.copy()
            try:
                current_pin = startpin.up
            except AttributeError:
                current_pin = startpin

        # visited pins stores all pins previously visited, which is used for detecting looped paths
        visited_pins = []
        # loose end stores the tracking state of all but one connection when multiple connections are found.
        loose_ends = []
        # boolean for decision-making:
        backtracking = False

        if not reverse:
            perform_step_1 = True
            perform_step_2 = False
        else:
            perform_step_1 = False
            perform_step_2 = True

        tracking = True

        # endpins stores the endpins of all paths traced, this is an output
        endpins = []

        # variables for output
        path_id = 0
        pinlist = []
        paths_dict = {}

        # before tracking, check if startpin is an endpin
        endpoint_reached, endpoint = self._is_endpoint(current_pin, endpoints)
        if endpoint_reached:
            # End of path. Stop performing path tracking of this path
            tracking = False

            # end of path, but we need to go up the instance stack in some way before output
            current_pin, instance_stack = self._drill_up(current_pin, instance_stack)
            # update logger and output lists
            logger.log(current_pin, endpoint=endpoint, status="endpoint")

            # if the endpoint is a pin, add that pin to endpins and paths_dict instead of the current pin
            if isinstance(endpoint, nd.Node):
                endpins.append(endpoint)
            else:
                endpins.append(current_pin)
            paths_dict[path_id] = pinlist

        for instance in instance_stack:
            instance.C2I_nodemap = {
                node.up: node for node in instance.pin.values()
            }  # Cell2Instance mapping for the instance nodes.

        logger.log(current_pin, status="start")

        # if reversing, we need to drill down to ensure that any pin2pin connection is found
        if reverse:
            logger.log(current_pin, status="reverse")
            current_pin, instance_stack = self._drill_down(
                current_pin,
                instance_stack,
                logger=logger,
            )

        while tracking:

            # if end of path is reached and there are still loose ends, track paths from one of the stored loose ends
            if backtracking:
                # obtain state of one of the loose ends
                (
                    current_pin,
                    instance_stack,
                    visited_pins,
                    pinlist,
                    perform_step_1,
                    perform_step_2,
                    tracker_kwargs,
                ) = loose_ends.pop()
                logger.log(
                    current_pin, status="backtracked", instance_stack=instance_stack
                )
                # stop backtracking and start tracking from the current pin
                backtracking = False

                # check if the pin backtracked to is actually an endpin:
                endpoint_reached, endpoint = self._is_endpoint(
                    current_pin, endpoints, instance_stack
                )
                if endpoint_reached:
                    # End of path. Stop performing path tracking of this path
                    perform_step_1 = False
                    perform_step_2 = False
                    # end of path, but we need to go up the instance stack in before output
                    current_pin, instance_stack = self._drill_up(
                        current_pin, instance_stack
                    )
                    # update logger and output lists
                    logger.log(current_pin, endpoint=endpoint, status="endpoint")

                    # if the endpoint is a pin, add that pin to endpins and paths_dict instead of the current pin
                    if isinstance(endpoint, nd.Node):
                        endpins.append(endpoint)
                    else:
                        endpins.append(current_pin)

                    paths_dict[path_id] = pinlist
                    path_id += 1

                    # check if there are any loose ends to track
                    if not loose_ends:
                        # no loose ends, stop tracking
                        break
                    else:
                        # loose ends, start backtracking
                        backtracking = True
                        logger.log(current_pin, status="backtracking")
                    # set current pin to None as we are currently not tracking from a pin
                    current_pin = None

            if perform_step_1:
                # STEP 1: Find optical connections
                tracker_kwargs_in = tracker_kwargs
                input_pin = current_pin
                (
                    segment_start,
                    optical_connections,
                    instance_stack,
                ) = self._find_optical_connections(
                    current_pin,
                    instance_stack,
                    logger=logger,
                    trackertype=trackertype,
                    tracker_kwargs=tracker_kwargs_in,
                )
                # deal with possibility of no optical connections
                if optical_connections is None:
                    # End of path. Stop performing path tracking of this path
                    perform_step_1 = False
                    perform_step_2 = False

                    # update logger and output lists
                    logger.log(current_pin, status="no_optical_connection")
                    logger.log(current_pin, status="end_of_path")
                    endpins.append(current_pin)
                    paths_dict[path_id] = pinlist
                    path_id += 1
                    # check if there are any loose ends to track
                    if not loose_ends:
                        # no loose ends, stop tracking
                        break
                    else:
                        # loose ends, start backtracking
                        backtracking = True
                        logger.log(current_pin, status="backtracking")
                    # set current pin to None as we are currently not tracking from a pin
                    current_pin = None

                # found optical connection(s), start tracking from (one of) the optical connections
                else:
                    # flow control: need to perform step 1 next
                    perform_step_1 = False
                    perform_step_2 = True

                    # move to new pin, update loose ends
                    (
                        current_pin,
                        connection_length,
                        direction,
                        posrel,
                        path,
                        sigtype,
                        tracker_kwargs_out,
                    ) = optical_connections.pop()
                    logger.log(
                        current_pin, status="optical_connection", L=connection_length
                    )
                    connection_length = (
                        ProtectedPartial(
                            filter_eval(connection_length), **tracker_kwargs
                        )
                        if callable(connection_length)
                        else connection_length
                    )
                    tracker_kwargs = tracker_kwargs_out

                    # update the pinlist for this path
                    pinlist_copy = pinlist.copy()
                    position_in = self.get_pin_coordinate_up_instance_stack(segment_start, instance_stack)
                    position_out = self.get_pin_coordinate_up_instance_stack(current_pin, instance_stack)
                    pinlist.append(
                        Segment(
                            segment_start,
                            current_pin,
                            connection_length,
                            position_in,
                            position_out,
                            upper_pin_in=input_pin,
                            tracker_kwargs_in=tracker_kwargs_in,
                            tracker_kwargs_out=tracker_kwargs_out,
                        )
                    )

                    # store remaining connections as loose ends that will be visited after path end is reached
                    for _ in optical_connections:
                        (
                            optical_connection,
                            connection_length,
                            direction,
                            posrel,
                            path,
                            sigtype,
                            tracker_kwargs_out,
                        ) = _

                        _pinlist_copy = pinlist_copy.copy()
                        position_in = self.get_pin_coordinate_up_instance_stack(segment_start, instance_stack)
                        position_out = self.get_pin_coordinate_up_instance_stack(optical_connection, instance_stack)
                        _pinlist_copy.append(
                            Segment(
                                segment_start,
                                optical_connection,
                                connection_length,
                                position_in,
                                position_out,
                                upper_pin_in=input_pin,
                                tracker_kwargs_in=tracker_kwargs_in,
                                tracker_kwargs_out=tracker_kwargs_out,
                            )
                        )
                        loose_ends.append(
                            [
                                optical_connection,
                                instance_stack.copy(),
                                visited_pins.copy(),
                                _pinlist_copy,
                                perform_step_1,
                                perform_step_2,
                                tracker_kwargs_out,
                            ]
                        )

                # If a current pin was found, perform checks:
                if current_pin is not None:
                    # Check if current pin was already visisted.
                    # generate a unique pin_id to differentiate between pins in different instances of the same cell
                    pin_id = (tuple(instance_stack), current_pin)
                    if pin_id in visited_pins:
                        # End of path. Stop performing path tracking of this path
                        perform_step_1 = False
                        perform_step_2 = False
                        # path is a loop, drill up from current pin
                        current_pin, instance_stack = self._drill_up(
                            current_pin, instance_stack
                        )

                        # update logger and output lists
                        logger.log(current_pin, status="loop")
                        endpins.append(current_pin)
                        # pinlist.append([current_pin, connection_length])
                        paths_dict[path_id] = pinlist
                        path_id += 1
                        # check if there are any loose ends to track
                        if not loose_ends:
                            # no loose ends, stop tracking
                            break
                        else:
                            # loose ends, start backtracking
                            backtracking = True
                            logger.log(current_pin, status="backtracking")
                        # set current pin to None as we are currently not tracking from a pin
                        current_pin = None
                    # if no loop, store the current pin for the sake of future loop detection
                    else:
                        visited_pins.append((tuple(instance_stack), current_pin))

                    # check if current pin is at an endpoint
                    if current_pin is not None:
                        endpoint_reached, endpoint = self._is_endpoint(
                            current_pin, endpoints, instance_stack
                        )
                        if endpoint_reached:
                            # End of path. Stop performing path tracking of this path
                            perform_step_1 = False
                            perform_step_2 = False
                            # end of path, but we need to go up the instance stack in some way before output
                            current_pin, instance_stack = self._drill_up(
                                current_pin, instance_stack
                            )
                            # update logger and output lists
                            logger.log(current_pin, endpoint=endpoint, status="endpoint")

                            # if the endpoint is a pin, add that pin to endpins and paths_dict
                            if isinstance(endpoint, nd.Node):
                                endpins.append(endpoint)
                            paths_dict[path_id] = pinlist
                            path_id += 1

                            # check if there are any loose ends to track
                            if not loose_ends:
                                # no loose ends, stop tracking
                                break
                            else:
                                # loose ends, start backtracking
                                backtracking = True
                                logger.log(current_pin, status="backtracking")
                            # set current pin to None as we are currently not tracking from a pin
                            current_pin = None

            if perform_step_2:
                # STEP 2: Find pin connections from current_pin
                (
                    connected_pins,
                    upper_output_pin,
                    instance_stack,
                ) = self._find_connected_pins(
                    current_pin, instance_stack, logger=logger
                )
                if reverse:
                    reverse = False
                else:
                    pinlist[-1].upper_pin_out = upper_output_pin
                # Found no connected pins
                if connected_pins is None:
                    # End of path. Stop performing path tracking of this path
                    perform_step_1 = False
                    perform_step_2 = False
                    # end of path, but we need to go up the instance stack in some way before output
                    current_pin, instance_stack = self._drill_up(
                        current_pin, instance_stack
                    )
                    # update logger and output lists
                    logger.log(current_pin, status="end_of_path")
                    endpins.append(current_pin)
                    # pinlist.append([current_pin, connection_length])
                    paths_dict[path_id] = pinlist
                    path_id += 1

                    # check if there are any loose ends to track
                    if not loose_ends:
                        # no loose ends, stop tracking
                        break
                    else:
                        # loose ends, start backtracking
                        backtracking = True
                        logger.log(current_pin, status="backtracking")
                    # set current pin to None as we are currently not tracking from a pin
                    current_pin = None

                # Found connections
                else:
                    # flow control: need to perform step 1 next
                    perform_step_1 = True
                    perform_step_2 = False

                    if len(connected_pins) == 1:
                        # move to new pin. For pin to pin connection the connection length is zero
                        connected_pin_1 = connected_pins.pop()

                        if not instance_stack:
                            # store the previous pin. We want to return this pin in the case of splitting nodes
                            # update current pin with the connected pin that we found.
                            current_pin = connected_pin_1
                            logger.log(current_pin, status="pin_connection")
                        # if we are not at the top level, try to go to the top level first
                        else:
                            try:
                                current_pin, instance_stack = self._drill_up(
                                    connected_pin_1, instance_stack
                                )
                                logger.log(current_pin, status="pin_connection")
                            except:
                                current_pin = connected_pin_1
                                logger.log(current_pin, status="pin_connection")
                    else:
                        # store other connected pins as loose ends, if splitting nodes are allowed
                        if allow_splitting_nodes:
                            nd.main_logger(f'Splitting node at pin {current_pin}, allow_splitting_nodes is set '
                                           f'to True.', 'warning')
                            connected_pin_1 = connected_pins.pop()

                            if not instance_stack:
                                # store the previous pin. We want to return this pin in the case of splitting nodes
                                # update current pin with the connected pin that we found.
                                current_pin = connected_pin_1
                                logger.log(current_pin, status="pin_connection")
                            # if we are not at the top level, try to go to the top level first
                            else:
                                try:
                                    current_pin, instance_stack = self._drill_up(
                                        connected_pin_1, instance_stack
                                    )
                                    logger.log(current_pin, status="pin_connection")
                                except:
                                    current_pin = connected_pin_1
                                    logger.log(current_pin, status="pin_connection")

                            for connected_pin in connected_pins:
                                # add this pin and connection length (0 for pin to pin connection) to copy of pinlist.
                                pinlist_copy = pinlist.copy()
                                # store the loose ends
                                loose_ends.append(
                                    [
                                        connected_pin,
                                        instance_stack.copy(),
                                        visited_pins.copy(),
                                        pinlist_copy,
                                        perform_step_1,
                                        perform_step_2,
                                        tracker_kwargs,
                                    ]
                                )
                        # if splitting nodes are not allowed (which is the default setting), we need to terminate this
                        # path and raise an error
                        else:
                            nd.main_logger(f'Splitting node at pin {current_pin}. Terminated path.','error')
                            perform_step_1 = False
                            perform_step_2 = False
                            # terminate path, but we need to go up the instance stack in some way before output
                            current_pin, instance_stack = self._drill_up(
                                current_pin, instance_stack
                            )
                            # update logger and output lists
                            logger.log(current_pin, status="end_of_path_splitting_node", connected_pins=connected_pins)
                            endpins.append(current_pin)
                            # pinlist.append([current_pin, connection_length])
                            paths_dict[path_id] = pinlist
                            path_id += 1

                            # check if there are any loose ends to track
                            if not loose_ends:
                                # no loose ends, stop tracking
                                break
                            else:
                                # loose ends, start backtracking
                                backtracking = True
                                logger.log(current_pin, status="backtracking")
                            # set current pin to None as we are currently not tracking from a pin
                            current_pin = None

                    # check if new current pin is at an endpoint
                    if current_pin is not None:
                        endpoint_reached, endpoint = self._is_endpoint(
                            current_pin, endpoints
                        )
                        if endpoint_reached:
                            # End of path. Stop performing path tracking of this path
                            perform_step_1 = False
                            perform_step_2 = False

                            # end of path, but we need to go up the instance stack in some way before output
                            current_pin, instance_stack = self._drill_up(
                                current_pin, instance_stack
                            )
                            # update logger and output lists
                            logger.log(current_pin, endpoint=endpoint, status="endpoint")

                            # if the endpoint is a pin, add that pin to endpins and paths_dict instead of the current pin
                            if isinstance(endpoint, nd.Node):
                                endpins.append(endpoint)
                                # pinlist.append([endpoint, connection_length])
                            else:
                                endpins.append(current_pin)
                                # pinlist.append([current_pin, connection_length])
                            paths_dict[path_id] = pinlist
                            path_id += 1

                            # check if there are any loose ends to track
                            if not loose_ends:
                                # no loose ends, stop tracking
                                break
                            else:
                                # loose ends, start backtracking
                                backtracking = True
                                logger.log(current_pin, status="backtracking")
                            # set current pin to None as we are currently not tracking from a pin
                            current_pin = None

        logger.log(pin=endpins, status="end_of_tracking", paths_dict=paths_dict)

        # return the endpins that were found during the tracking
        return endpins, paths_dict


    def _generate_endpins_dictionary(
            self,
            paths_dict,
            startpin,
            instance_stack=None,
    ):
        """
        Function to generate an endpins dictionary that contains the endpins and their position at the highest
        hierarchical level.

        args:
            paths_dict (dict): Dictionary of the paths, obtained from self._findpaths().
            startpin (nd.Pin): Startpin from which tracking is performed.
            instance_stack (list): List of instances above the startpin.

        returns:
            endpins_dict (dict): dictionary containing a dictionary of the endpin and position for each path.
        """
        endpins_dict = {}

        if instance_stack is None:
            instance_stack = []

        for i, segments in paths_dict.items():
            # if paths were found, use the segments to obtain each endpin and position
            if len(segments) > 0:
                endpins_dict[i] = {
                    'endpin': segments[-1].pin_out,
                    'position': segments[-1].position_out
                }
            # if no paths were found (termination on startpin), extract the position of the startpin.
            else:
               endpins_dict[i] = {
                    'endpin': startpin,
                    'position': self.get_pin_coordinate_up_instance_stack(startpin, instance_stack)
                }

        return endpins_dict

    def _get_lengths(
        self,
        segments=None,
    ):
        """
        Get the lengths of all the paths stored in a paths dictionary

        Args:
            segments (dict): A dictionary containing the paths found using _findpaths()
        Returns:
            lengths (list): A list of all the pathlengths.
        """
        if segments is None:
            segments = self.segments

        lengths = []
        for i in range(len(segments)):
            lengths.append(pathlength(segments[i]))

        return lengths

    def _drill_down(
            self,
            start,
            instance_stack=None,
            logger=None,
    ):
        """Drill down from a node in the cell hierarchy until the max depth is reached

        """
        if instance_stack is None:
            instance_stack = []

        if start.cnode.instance is not None:  # instance pin start
            current_pin = start
            instance = start.cnode.instance
            instance_stack.append(instance)
        else:  # cell pin start
            current_pin = start

        while True:
            if current_pin.cnode.instance is not None:  # pin resides in an instance
                instance.C2I_nodemap = {
                    node.up: node for node in instance.pin.values()
                }  # Cell2Instance mapping for the instance nodes.
                cell_node = (
                    current_pin.up
                )  # go to corresponding cell node of current_pin before drilling down
            else:
                cell_node = current_pin

            #   Do not drill down in auxiliary cells that have (by definition) no netlist function.
            for nb, trans in cell_node.nb_geo:
                if nb.cnode.instance is not None and not nb.cnode.cell.auxiliary:
                    break  # found instance level below represented via node nb.
                    # assume max *one* instance below can be connected
            if nb.cnode.instance is None:
                break  # bottom reached, leave while loop
            else:
                current_pin = nb
                instance = nb.cnode.instance
                instance_stack.append(instance)
                if logger is not None:
                    logger.log(current_pin, status="drill_down")

        start_segment = cell_node
        return start_segment, instance_stack

    def _find_optical_connections(
        self,
        start,
        instance_stack=None,
        logger=None,
        flip=1,
        trackertype="dis",
        end_cell_name="",
        tracker_kwargs=None,
    ):
        """Drill down from a node in the cell hierarchy until a valid edge is found or max depth is reached.

        This function finds if node <start> has an optical edge.
        If not, it is checked if node is connected to an instance to, if so,
        drill down into it.
        When drilling down this function creates a "C2I_nodemap" attribute (cell2instance)
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
        # TODO improve this docstring, change names of variables to be more descriptive, add comments explaining how
        #  this works
        # TODO refactor?

        if instance_stack is None:
            instance_stack = []

        if start.cnode.instance is not None:  # instance pin start
            current_pin = start
            instance = start.cnode.instance
            instance_stack.append(instance)
            if instance.cnode.flip:
                flip *= -1
        else:  # cell pin start
            current_pin = start

        # drill down until a signal edge is found or the bottom is reached:
        while True:
            next_nodes_opt = []  # store signal connections.

            if current_pin.cnode.instance is not None:  # pin resides in an instance
                instance.C2I_nodemap = {
                    node.up: node for node in instance.pin.values()
                }  # Cell2Instance mapping for the instance nodes.
                cell_node = (
                    current_pin.up
                )  # go to corresponding cell node of current_pin before drilling down
                if cell_node.cnode.cell.cell_name == end_cell_name:
                    break  # end of path, leave while loop.
            else:
                cell_node = current_pin

            # Find (and store) all signal paths of cell_node, if any:
            for nb, L, direction, sigtype, path, extra_out in cell_node.path_nb_iter(
                sigtype=trackertype,
                extra=tracker_kwargs,
            ):
                # find xya position of the optical nodes for visualisation
                if nb is None:
                    x, y, a = 0, 0, 0  # termination
                else:
                    x, y, a = nd.diff(cell_node, nb)
                posrel = nd.Pointer(x, y, a)
                if flip == -1:
                    posrel.flip()
                next_nodes_opt.append(
                    (nb, L, direction, posrel.xya(), path, sigtype, extra_out)
                )
            if len(next_nodes_opt) > 0:
                if logger is not None:
                    logger.log(cell_node, status="drill_down")
                break  # found optical link, leave while loop

            # If no optical edge was found, drill down, if possible.
            #   Do not drill down in auxiliary cells that have (by definition) no netlist function.
            for nb, trans in cell_node.nb_geo:
                if nb.cnode.instance is not None and not nb.cnode.cell.auxiliary:
                    break  # found instance level below represented via node nb.
                    # assume max *one* instance below can be connected
            if nb.cnode.instance is None:
                break  # bottom reached, leave while loop
            else:
                current_pin = nb
                instance = nb.cnode.instance
                if instance.cnode.flip:
                    flip *= -1
                instance_stack.append(instance)
                if logger is not None:
                    logger.log(current_pin, status="drill_down")

        start_segment = cell_node
        if len(next_nodes_opt) == 0:
            return [None] * 3
        else:
            return start_segment, next_nodes_opt, instance_stack

    def _find_connected_pins(
        self,
        pin,
        instance_stack,
        logger,
    ):
        """
        Function that looks for pin-to-pin connections to the input pin.

        Args:
            pin (Node): pin to find pin-to-pin connections for.
            instance_stack (list): List of parent instances/cells of the input pin
            logger (Logger): Logger object used for logging purposes

        Returns:
            adjacent_pins (list): list of all pins connected to the input pin. If no connected pins are found,
                                    this returns None.
            instance_stack (list):List of parent instances/cells of the output pins.  If no connected pins are found,
                                    returns the initial instance stack.
        """

        # save the instance stack. If no connect pins are found, we need to output the old instance stack
        old_instance_stack = instance_stack.copy()

        # while there are instances in the instance stack, go up in hierarchy and look for pin connections
        while instance_stack:
            adjacent_pins = []
            parent_pins = []
            # move up one step
            instance = instance_stack[-1]
            pin = instance.C2I_nodemap[pin]

            # look at all pins connected to this pin
            for node in self._nb_filter_geo(pin):
                # if node is in a different instance with the same parent, it is a connected node
                if (
                    node.cnode.instance is not pin.cnode.instance
                    and node.cnode.parent_cnode is pin.cnode.parent_cnode
                    and node.width != 0
                ):
                    adjacent_pins.append(node)
                # otherwise, it is possible that the node is a parent node.
                elif (
                    pin.cnode.parent_cnode is node.cnode
                    and pin.x == node.x
                    and pin.y == node.y
                ):
                    parent_pins.append(node)

            # if we find a connected pin, return that connected pin and an updated instance stack
            if adjacent_pins:
                instance_stack.pop()
                logger.log(pin, status="drill_up")
                return adjacent_pins, pin, instance_stack

            # if we don't find a connected pin but do find a parent node,
            # move to parent node and update the instance stack
            elif parent_pins:
                pin = parent_pins[0]
                instance = instance_stack.pop()

                logger.log(pin, status="drill_up")

                # if instance stack is empty, that means there is no pin to pin connection.
                if instance_stack == []:
                    logger.log(pin, status="no_pin_connection")
                    # return None for adjacent pins and the old instance stack.
                    adjacent_pins = None
                    instance_stack = old_instance_stack
                    return adjacent_pins, pin, instance_stack

            # if we don't find any adjacent nodes or parent nodes, that means there are no pin to pin connections
            else:
                logger.log(pin, status="no_pin_connection")
                # return None for adjacent pins and the old instance stack.
                adjacent_pins = None
                instance_stack = old_instance_stack
                return adjacent_pins, pin, instance_stack

        logger.log(pin, status="no_pin_connection")
        # return None for adjacent pins and the old instance stack.
        adjacent_pins = None
        instance_stack = old_instance_stack
        return adjacent_pins, pin, instance_stack

    def _drill_up(
        self,
        pin,
        instance_stack,
    ):
        """
        Function to go to the highest hierarchical level saved in the instance stack corresponding to pin.

        Args:
            pin (Node): pin from which we move up.
            instance_stack (list): List of parent instances/cells of the input pin.
        Returns:
            pin (Node): top level pin.
            instance_stack (list): Updated list of parent instances/cells of the input pin.
        """

        # continue trying to go up in hierarchy until the instance_stack is empty
        while instance_stack != []:
            # consider the next instance:
            instance = instance_stack[-1]
            try:
                pin = instance.C2I_nodemap[pin]
            # if the top level is the instance of a cell, it's not possible to find a pin one level up for the
            # last entry of the instance stack
            except:
                break

            # update the instance stack
            instance_stack.pop()

            # look at all pins connected to this pin
            for node in self._nb_filter_geo(pin):
                # if the node is a parent of the pin, consider that node
                if pin.cnode.parent_cnode is node.cnode:
                    # if the instance stack is empty, we've reached the top level
                    if instance_stack == []:
                        return node, instance_stack
                    # if the instance stack is not empty, move to the parent node
                    else:
                        pin = node

        # return the updated pin and instance stack
        return pin, instance_stack

    def _is_endpoint(self, pin, endpoints, input_instance_stack=None):
        """
        Function to check if a pin is an endpoint. A pin is an endpoint if it or its parent pins is the same as one of
        the endpoints, or is in a cell that is the same as one of the endpoints, or if it is in an instance that is the
        same as one of the endpoints.

        First, this functions drills all the way down the hierarchy, while storing all instances traversed. Then, this
        drills all the way up through all the instances traversed (and possible instances provided as
        input_instance_stack), checking at each corresponding pin whether it is an endpoint. If an endpoint is
        encountered, this function returns True. If no endpoint is encountered in any of the instances, this function
        returns False.

        Args:
            pin (Node): pin to check for
            endpoints (list): list of endpoints, which can be pins, cells, or instances.
            input_instance_stack (list): list containing a history of traversed instances

        Returns:
            True, endpoint: if the pin is an endpoint, return True and the endpoint
            False, None: if the pin is not an endpoint, return False and None
        """
        # TODO consider refactoring

        if endpoints is None:
            return False, None

        # if endpoints were specified, look if the pin is one of the endpoints
        if endpoints is not None:

            # convert endpoints to a list if only a single endpoint was provided
            try:
                len(endpoints)
            except:
                endpoints = [endpoints]

            # for each endpoint, check if the pin is an endpoint or the pin's cell is an endpoint
            for endpoint in endpoints:
                # generate instance_stack if necessary
                if input_instance_stack is None:
                    instance_stack = []
                else:
                    instance_stack = input_instance_stack.copy()

                # setup for drilling down
                if pin.cnode.instance is not None:  # instance pin start
                    current_pin = pin
                    instance = pin.cnode.instance
                    instance_stack.append(instance)
                else:  # cell pin start
                    current_pin = pin

                # drill down
                while True:
                    if (
                        current_pin.cnode.instance is not None
                    ):  # pin resides in an instance
                        instance.C2I_nodemap = {
                            node.up: node for node in instance.pin.values()
                        }  # Cell2Instance mapping for the instance nodes.
                        cell_node = (
                            current_pin.up
                        )  # go to corresponding cell node of current_pin before drilling down
                    else:
                        cell_node = current_pin

                    # Do not drill down in auxiliary cells that have (by definition) no netlist function.
                    for nb, trans in cell_node.nb_geo:
                        if (
                            nb.cnode.instance is not None
                            and not nb.cnode.cell.auxiliary
                        ):
                            break  # found instance level below represented via node nb.
                            # assume max *one* instance below can be connected
                    if nb.cnode.instance is None:
                        if (
                            current_pin.cnode.instance is not None
                        ):  # pin resides in an instance
                            instance.C2I_nodemap = {
                                node.up: node for node in instance.pin.values()
                            }  # Cell2Instance mapping for the instance nodes.
                            pin = (
                                current_pin.up
                            )  # go to corresponding cell node of current_pin before drilling down
                        else:
                            pin = current_pin

                        break  # bottom reached, leave while loop
                    else:
                        current_pin = nb
                        instance = nb.cnode.instance
                        instance_stack.append(instance)

                # next, drill up and check if any pin there is an endpoint
                while instance_stack != []:
                    if isinstance(endpoint, nd.Node):
                        if pin == endpoint:
                            return True, endpoint
                    elif isinstance(endpoint, nd.Cell):
                        if pin.cnode.cell.cell_name == endpoint.cell_name:
                            instance = instance_stack[-1]
                            return True, instance.C2I_nodemap[pin]
                    elif isinstance(endpoint, str):
                        if pin.cnode.cell.cell_name == endpoint:
                            instance = instance_stack[-1]
                            return True, instance.C2I_nodemap[pin]
                    elif isinstance(endpoint, nd.Instance):
                        if pin.cnode.instance == endpoint:
                            return True, pin
                    else:
                        msg = f"Provided endpoint {endpoint} is not a pin, cell or instance."
                        raise Exception(msg)

                    # consider the next instance:
                    instance = instance_stack[-1]
                    try:
                        pin = instance.C2I_nodemap[pin]
                    except:
                        break

                    if isinstance(endpoint, nd.Node):
                        if pin == endpoint:
                            return True, endpoint
                    elif isinstance(endpoint, nd.Cell):
                        if pin.cnode.cell.cell_name == endpoint.cell_name:
                            instance = instance_stack[-1]
                            return True, instance.C2I_nodemap[pin]
                    elif isinstance(endpoint, str):
                        if pin.cnode.cell.cell_name == endpoint:
                            instance = instance_stack[-1]
                            return True, instance.C2I_nodemap[pin]
                    elif isinstance(endpoint, nd.Instance):
                        if pin.cnode.instance == endpoint:
                            return True, pin

                    # look at all pins connected to this pin
                    parent_pins = []
                    for node in self._nb_filter_geo(pin):
                        # if the node is a parent of the pin, consider that node
                        if (pin.cnode.parent_cnode is node.cnode
                            and pin.x == node.x
                            and pin.y == node.y
                        ):
                            # store the parent node
                            parent_pins.append(node)

                            # check if the node itself is an endpoint
                            if node == endpoint:
                                return True, endpoint
                            # check if the node resides in an instance. That cell could be the endpoint
                            elif node.cnode.instance is not None:
                                node = node.up
                                try:
                                    if node.cnode.cell.cell_name == endpoint.cell_name:
                                        return True, endpoint
                                except:
                                    pass
                                # check if pin's instance is an endpoint
                                try:
                                    if node.cnode.instance == endpoint:
                                        return True, endpoint
                                except:
                                    pass

                    # if any parent pins were found, move to the first parent pin that was found
                    if parent_pins != []:
                        instance_stack.pop()
                        pin = parent_pins[0]

            # if none of the pins where endpoints, return false
            return False, None

    @staticmethod
    def _nb_filter_geo(
        start,
    ):
        """Yield all neighbouring nodes except for the origin, stubs and annotations.

        Iterates over the neighbours in the Cell.

        Args:
            start (Node): a Node from which you want the neighbours

        Yields:
            start_nb (Node): iterator over start's filtered neighbours
        """
        for start_nb, _ in start.nb_geo:
            if "org" not in start_nb.name:
                yield start_nb


def pathlength(path: List[Segment], **kwargs):
    """Returns the length of the paths

    The kwargs provided are fed connections if they are functions
    """
    lengths = [
        seg.connection(**kwargs) if callable(seg.connection) else seg.connection
        for seg in path
    ]
    return sum(lengths)


def findpath(
    start,
    end=None,
    end_cell_name=None,
    stdout=False,
    pathfilename=None,
    append=False,
    log=False,
    logfilename="trace.log",
    width=2.0,
    layer=2000,
    tracker="dis",
    show=True,
    paths2cell=None,
):
    """
    Old function that should not be used anymore. This version will raise a warning to inform the user that the
    new wrapper should be used instead.

    The behaviour is simular to the old one, in the sense that it returns a list of segments for each path

    Raises an error when called.
    """
    nd.logger.error(
        "Findpath function is deprecated and will be removed. Use the class Paths() instead."
    )

    paths = Paths(
        start,
        endpoints=end or end_cell_name,
        logfilename=logfilename,
        logging_level=0,
        reverse=False,
        trackertype=tracker,
        instance_stack=None,
        tracker_kwargs=None,
    )

    paths.show_paths(
        layer_number=layer,
        width=width,
        paths2cell=paths2cell,
    )

    return paths.segments


if __name__ == "__main__":
    pass
