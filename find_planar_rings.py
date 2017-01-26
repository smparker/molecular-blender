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

from .util import stopwatch, Timer

import bpy

# Creates a graph structure from list of atoms and bonds
def createGraph(atomlist,connectivitylist):
    graph = [ set() for x in range(len(atomlist)) ]

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

# Basic depth first search algorithm, pulled from the Internet
def DFS(graph, start, goal, pathrange=(2,7)):
    stack = [(start,[start])]
    while stack:
        (vertex,path) = stack.pop()
        # any successes after this point will add one to current length of path
        if len(path)+1 < pathrange[1]:
            for nextvert in graph[vertex] - set(path):
                if nextvert == goal:
                    if pathrange[0] < len(path)+1 < pathrange[1]:
                        yield path + [nextvert]
                else:
                    stack.append((nextvert,path + [nextvert]))

# Finds all unique paths of a vertex back to itself
def findSpecificCycles(graph, vertex, pathrange=(2, 7)):
    cycles = []
    fixedcycles = []
    for goal in graph[vertex]:
        for path in DFS(graph,vertex,goal,pathrange=pathrange):
            if set(path) not in cycles:
                cycles.append(set(path))
                fixedcycles.append(path)
    return cycles,fixedcycles

# Finds all closed walks of with lengths bounded by path range. Defaults to paths
# of length > 2 and < 7
def findAllUniqueCycles(graph, pathrange=(2,7)):
    cycles = []
    for i in range(len(graph)):
        possibleCycles,extras = findSpecificCycles(graph,i,pathrange)
        for cycle in possibleCycles:
            if set(cycle) not in cycles and pathrange[0] < len(cycle) < pathrange[1]:
                cycles.append(set(cycle))
    return cycles

# Given a list of atoms and bonds, determine where aromatic rings are and plot them
def plotRings(context, molecule, options):
    timer = Timer()
    atomlist = molecule.atoms
    bondlist = [(b.iatom.index, b.jatom.index) for b in molecule.bonds]
    timer.tick_print("generate lists")
    graph = createGraph(atomlist,bondlist)
    timer.tick_print("create graphs")
    cycles = findAllUniqueCycles(graph)
    timer.tick_print("find unique cycles")
    planarCycles = [x for x in cycles if isPlanar(atomlist,x)]
    timer.tick_print("check cycle planarity")

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
    timer.tick_print("plot all planar cycles")
