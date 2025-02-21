# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: nodes.py
#  Copyright (C) 2019-2025 Shane Parker
#
#  This file is part of Molecular Blender.
#
#  Molecular Blender is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  Molecular Blender is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Molecular Blender; see COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#

"""Utility functions to manage material nodes"""

def float_lerp(a, b, t):
    return (1.0 - t) * a + t * b

def sum_heights(nodes_array):
    result = 0
    for node in nodes_array:
        result = result + node.height
    return result

def sum_widths(depth_nodes):
    result = 0
    for depth in depth_nodes:
        max_width = 0
        for node in depth_nodes[depth]:
            if max_width < node.width:
                max_width = node.width
        result = result + max_width
    return result

def calc_priority_by_type(node):
    if node.type == 'NEW_GEOMETRY' or node.type == 'TEX_COORD' or node.type == 'GROUP_INPUT':
        return -6
    if node.type == 'VALUE' or node.type == 'ATTRIBUTE':
        return -5
    if node.type == 'SEPXYZ':
        return -4
    if node.type == 'SEPHSV' or node.type == 'SEPRGB' or node.type == 'BLACKBODY':
        return -3
    if node.type == 'MATH' or node.type == 'VECT_MATH':
        return -2
    if node.type == 'COMBXYZ':
        return -1
    if node.type == 'COMBHSV' or node.type == 'COMBRGB':
        return 1
    if node.type == 'MIX_RGB' or node.type == 'HUE_SAT':
        return 2
    if node.type == 'TEX_IMAGE' or node.type == 'TEX_MUSGRAVE' or node.type == 'TEX_BRICK' or node.type == 'TEX_NOISE' or node.type == 'TEX_VORONOI':
        return 3
    if node.type == 'BSDF_DIFFUSE' or node.type == 'BSDF_PRINCIPLED' or node.type == 'EMISSION':
        return 4
    if node.type == 'HOLDOUT' or node.type == 'VOLUME_SCATTER' or node.type == 'VOLUME_ABSORPTION':
        return 5
    if node.type == 'MIX_SHADER':
        return 6
    if node.type == 'OUTPUT_MATERIAL' or node.type == 'OUTPUT_LAMP' or node.type == 'GROUP_OUTPUT':
        return 7

    return 0

def calc_priority_by_socket(node):
    if len(node.inputs) == 0:
        return -9999
    if len(node.outputs) == 0:
        return 9999

    result = 0
    for in_socket in node.inputs:
        if in_socket.is_linked:
            for link in in_socket.links:
                if link.is_valid:
                    if len(link.from_node.inputs) == 0:
                        result -= 1
                    else:
                        result += 2

    for out_socket in node.outputs:
        if out_socket.is_linked:
            for link in out_socket.links:
                if link.is_valid:
                    if len(link.to_node.outputs) == 0:
                        result += 10
                    else:
                        result -= 1

    return result

priority_calcs = { "type" : calc_priority_by_type, "socket" : calc_priority_by_socket }

def arrange_nodes(node_array, calc_priority_by, horiz_padding=0.125, vert_padding=0.125):
    # Create a dictionary where the key is the
    # depth and the value is an array of nodes.
    calc_priority = priority_calcs[calc_priority_by]

    depth_nodes = {}
    for node in node_array:

        depth = calc_priority(node)
        if depth in depth_nodes:

            # Add the node to the node array at that depth.
            depth_nodes[depth].append(node)
        else:

            # Begin a new array.
            depth_nodes[depth] = [node]

    # Add padding to half the width.
    extents_w = (0.5 + horiz_padding) * sum_widths(depth_nodes)
    t_w_max = 0.5
    sz0 = len(depth_nodes)
    if sz0 > 1:
        t_w_max = 1.0 / (sz0 - 1)

    # List of dictionary KVPs.
    depths = sorted(depth_nodes.items())
    depths_range = range(0, sz0, 1)
    for i in depths_range:
        nodes_array = depths[i][1]
        t_w = i * t_w_max
        x = float_lerp(-extents_w, extents_w, t_w)

        extents_h = (0.5 + vert_padding) * sum_heights(nodes_array)
        t_h_max = 0.5
        sz1 = len(nodes_array)
        if sz1 > 1:
            t_h_max = 1.0 / (sz1 - 1)

        top = t_h_max * (sz1 - 1)
        nodes_range = range(0, sz1, 1)
        for j in nodes_range:
            node = nodes_array[j]
            t_h = top - j * t_h_max
            y = float_lerp(-extents_h, extents_h, t_h)
            half_w = 0.5 * node.width
            half_h = 0.5 * node.height
            node.location.xy = (x - half_w, y - half_h)
