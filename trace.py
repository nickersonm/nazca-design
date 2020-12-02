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

# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
#==============================================================================
# (c) 2017 Ronald Broeke
#==============================================================================

"""Module for tracing path lengths.

The example below shows how to start and stop a trace and how to query it
for the length of all the components (cells) that are put in the layout between the
start and stop calls. Note that cells must have a length_geo attribute, which is
the case for at least all components in mask_elements and all interconnects.

Example::

    import nazca as nd

    nd.trace.trace_start()
    nd.strt(length=100).put()
    nd.bend(angle=90).put()
    nd.trace.trace_stop()
    print(nd.trace.trace_length())

"""


from collections import defaultdict

trace_closed = None
trace_cnt = -1
trace_id_list = []
traces = defaultdict(list)
infolevel = 0

def trace_start(name=None):
    """Start recording a trace.

    If no trace <name> is provided an ordinal number will be selected.

    Args:
        name (str | int): optional id of the trace to stop

    Returns:
        None
    """
    global trace_id, trace_cnt, trace_ids
    trace_cnt += 1

    if name is None:
        while trace_cnt in trace_id_list:
            trace_cnt += 1
        name = trace_cnt

        if infolevel > 0:
            if len(trace_id_list) == 0:
                print('\n')
            print("{}trace open: {}".format('  '*len(trace_id_list), trace_cnt))

    elif name in trace_id_list:
        print("Warning: starting already active trace_id = {} in trace_start.".\
           format(name))

        return
    trace_id_list.append(name)
    #print('tracelist:', trace_id_list)


def trace_stop(name=None):
    """Stop recording on trace.

    If no trace_id is given the last started trace is stopped.

    Args:
        trace_id (str | int): id of the trace to stop

    Returns:
        None
    """
    global trace_closed
    if name is None:
        name = trace_id_list[-1]
    #print('stop trace:', name)
    trace_closed = name
    if name == trace_id_list[-1]:
        del trace_id_list[-1]
        #try:
        #    trace_id_list[-1]
        #except:
        #    pass
        if infolevel > 0:
            print('{}trace close: {}, size: {}'.\
                format('  '*len(trace_id_list), name, len(traces[name])))
    else:
        print("Warning: Trying to stop already stopped or not existing trace id={}.".\
            format(trace_id))
    return None


def trace_append(elm):
    """Add an element to the active trace."""
    traces[trace_id_list[-1]].append(elm)


def trace_length(name=None):
    """Calculate the trace length.

    If no trace name is provided in <name> the last closed trace will be used.

    Returns:
        float: trace length
    """
    if name is None:
        name = trace_closed
    length = 0
    if infolevel > 0:
        print('{}trace length: {}'.format('  '*len(trace_id_list), name))
    for elm in traces[name]:
        try:
            L = elm.cnode.cell.length_geo
            length += L
            if infolevel > 0:
                print("{0}{1}: {2}, {3}".format('  '*len(trace_id_list), name, elm, L))
        except AttributeError:
            if infolevel > 0:
                print("{}{}: No trace elements found".format('  '*len(trace_id_list), name))

    if infolevel > 0:
        print("  length = {}".format(length))
    return length
