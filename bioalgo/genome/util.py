from collections import Counter

def all_kmers_of(text, k):
	for i in range(len(text) - k + 1):
		yield text[i:i+k]

def eulerian_path(edges):
	inDegree = Counter()
	outDegree = Counter()
	for edge in edges:
		outDegree[edge[0]] += 1
		inDegree[edge[1]] += 1

	start = [i for i in outDegree.keys() if outDegree[i] - inDegree[i] == 1][0]
	end = [i for i in inDegree.keys() if inDegree[i] - outDegree[i] == 1][0]

	edges.append((end, start))

	eulerianCycle = eulerian_cycle(edges)
	return make_eulerian_path(eulerianCycle, start, end)

def eulerian_cycle(edges):
	start = edges[0][0]
	cycle = form_a_cycle(start, edges)
	while edges:
		newStart = select_new_starting_point_in_cycle(cycle, edges)
		cycle = rotate_cycle(cycle, newStart)
		cycle = cycle + form_a_cycle(newStart, edges)

	return cycle

def make_eulerian_path(eulerianCycle, start, end):
	for i, path in enumerate(eulerianCycle):
		if path == (end, start):
			return eulerianCycle[i+1:] + eulerianCycle[:i]

def form_a_cycle(start, edges):
	curr = start
	cycle = []
	while True:
		edge = randomly_pop_edge_starting_from(curr, edges)
		curr = edge[1]
		cycle.append(edge)
		if curr == start:
			break

	return cycle

def randomly_pop_edge_starting_from(start, edges):
	for edge in edges:
		if edge[0] == start:
			edges.remove(edge)
			return edge

def select_new_starting_point_in_cycle(cycle, edges):
	for edge in edges:
		if edge[0] in [edge[0] for edge in cycle]:
			return edge[0]

def rotate_cycle(cycle, newStart):
	for i, edge in enumerate(cycle):
		if edge[0] == newStart:
			cycle = cycle[i:] + cycle[:i]
			break

	return cycle
