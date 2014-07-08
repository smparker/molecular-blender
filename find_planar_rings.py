#  Molecular Blender
#  Filename: find_planar_rings.py
#  Copyright (C) 2014 Shane Parker, Joshua Szekely
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

import bpy

# Creates a graph structure from list of atoms and bonds
def createGraph(atomlist,connectivitylist):
    graph = {}
    for x in range(len(atomlist)+1):
        graph[x] = set()
    for x,y in connectivitylist:
        graph[x].update(set([y]))
        graph[y].update(set([x]))
    return graph

# Check if vertices of loop are coplanar
def isPlanar(atomlist,cycle):
    if len(cycle) is 3:
        return True
    else:
        testlist = list(cycle)
        for x in range(3,len(testlist)):
            C1 = atomlist[ testlist[x-3] ].position
            C2 = atomlist[ testlist[x-2] ].position
            C3 = atomlist[ testlist[x-1] ].position
            C4 = atomlist[ testlist[x-0] ].position
            threshold = (C3-C1).dot((C2-C1).cross(C4-C3))
            if abs(threshold) > 0.1:
                return False
        return True

# Basic depth first search algorithm, from pulled from the Internet
def DFS(graph,start,goal):
    stack = [(start,[start])]
    while stack:
        (vertex,path) = stack.pop()
        for next in graph[vertex] - set(path):
            if next == goal:
                yield path + [next]
            else:
                stack.append((next,path + [next]))

# Finds all unique paths of a vertex back to itself
def findSpecificCycles(graph,vertex):
    cycles = []
    fixedcycles = []
    for goal in graph[vertex]:
        paths = list(DFS(graph,vertex,goal))
        for path in paths:
            if set(path) not in cycles:
                cycles.append(set(path))
                fixedcycles.append(path)
    return cycles,fixedcycles

# Finds all closed walks of length > 2, less than 7
def findAllUniqueCycles(graph):
    cycles = []
    for x in list(graph):
        possibleCycles,extras = findSpecificCycles(graph,x)
        for cycle in possibleCycles:
            if set(cycle) not in cycles and 2 < len(cycle) < 7:
                cycles.append(set(cycle))
    return cycles

# Given a list of atoms and bonds, determine where aromatic rings are and plot them
def plotRings(context,atomlist,bondlist,options):
    graph = createGraph(atomlist,bondlist)
    cycles = findAllUniqueCycles(graph)
    planarCycles = [x for x in cycles if isPlanar(atomlist,x)]
    for num,pCycle in enumerate(planarCycles):
        objname  = "Ring" + str(num)
        meshname = objname + "mesh"
        ringMesh = bpy.data.meshes.new(meshname)
        verts = []
        bonds = []
        size = len(pCycle)
        # Get the properly ordered cycle
        extras,fixedcycles = findSpecificCycles(graph,list(pCycle)[0])
        fpCycle = []
        for fc in fixedcycles:
            if set(fc) == pCycle:
                fpCycle = fc
                break
        for atomIndex,meshIndex in zip(fpCycle,range(size)):
            verts.append(atomlist[atomIndex].position)
            bonds.append([meshIndex, (meshIndex+1)%size])
        ringMesh.from_pydata(verts,bonds,[range(size)])
        ringMesh.update()
        ringObj = bpy.data.objects.new(objname,ringMesh)
        ringObj.data = ringMesh
        context.scene.objects.link(ringObj)