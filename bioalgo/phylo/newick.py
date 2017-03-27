from bioalgo.phylo.newicktree import NewickTree, NewickNode

def parse(string):
	root, _ = read_internal_node(string, 0)
	return NewickTree(root)

def read_internal_node(string, curr):
	assert string[curr] == '(', 'Invalid format.'
	node = NewickNode()

	# parse through the string representing child nodes
	curr += 1
	while string[curr] not in [')', ';']:
		# another child exists
		if string[curr] == ',':
			curr += 1
		# another internal-node child exists
		elif string[curr] == '(':
			child, curr = read_internal_node(string, curr)
			node.add_child(child)
		# another leaf-node child exists
		else:
			child, curr = read_leaf_node(string, curr)
			node.add_child(child)

	# set name and distance of the node
	curr += 1
	name, distance, curr = parse_name_and_distance(string, curr)
	node.set_name(name)
	node.set_distance(distance)

	return node, curr

def read_leaf_node(string, curr):
	node = NewickNode()
	
	name, distance, curr = parse_name_and_distance(string, curr)
	node.set_name(name)
	node.set_distance(distance)

	return node, curr

def parse_name_and_distance(string, curr):
	readingName = True
	name, distance = None, None
	characters = []
	while string[curr] not in [',', ')', ';']:
		if string[curr] == ':':
			name = ''.join(characters) if characters else None
			readingName = False
			characters = []
			curr += 1
			continue

		characters.append(string[curr])
		curr += 1

	if readingName:
		name = ''.join(characters) if characters else None
	else:
		distance = float(''.join(characters))

	return name, distance, curr