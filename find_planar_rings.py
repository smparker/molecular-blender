# Creates a graph structure from list of atoms and bonds
def createGraph(atomlist,connectivitylist):
	graph = {}
	for x in range(len(atomlist)):
		graph[x+1] = set()
	for x,y in connectivitylist:
		graph[x].update(set([y]))
		graph[y].update(set([x]))
	return graph

#Check if vertices of loop are coplanar
def isPlanar(atomlist,cycle):

# Find center of mass of loop
def findCenter(atomlist,cycle):

# Given a list of atoms and bonds, determine where aromatic rings are and plot them
def plotRings(atomlist,bondlist):


#test case: two fused 5-membered rings
graph = {1: set([2,5,8]),
		2: set([1,3]),
		3: set([2,4]),
		4: set([3,5]),
		5: set([1,4,6]),
		6: set([5,7]),
		7: set([6,8]),
		8: set([1,7])}
atoms = [1,2,3,4,5,6,7,8]
bonds = [(1,2),(2,3),(3,4),(4,5),(1,5),(1,8),(8,7),(7,6),(6,5)]

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
	for goal in graph[vertex]:
		paths = list(DFS(graph,vertex,goal))
		for path in paths:
			if set(path) not in cycles:
				cycles.append(set(path))
	return cycles

# Finds all closed walks of length > 2, less than 7
def findAllUniqueCycles(graph):
	cycles = []
	for x in list(graph):
		possibleCycles = findSpecificCycles(graph,x)
		for cycle in possibleCycles:
			if set(cycle) not in cycles and 2 < len(cycle) < 7:
				cycles.append(set(cycle))
	return cycles

findAllUniqueCycles(graph)
