class NewickTree:
	
	def __init__(self, root):
		self.root = root

class NewickNode:

	def __init__(self):
		self.name = None
		self.children = []
		self.distance = None

	def set_name(self, name):
		self.name = name

	def set_distance(self, distance):
		self.distance = distance

	def add_child(self, child):
		self.children.append(child)