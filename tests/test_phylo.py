import sys; sys.path.append('../')
import bioalgo.phylo.newick as newick
import unittest

class NewickTest(unittest.TestCase):

	def setUp(self):
		pass

	def test_read_internal_node(self):
		node, curr = newick.read_internal_node('(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;', 0)
		self.assertEqual(node.name, 'F')
		self.assertEqual([n.name for n in node.children], ['A', 'B', 'E'])
		self.assertEqual([n.name for n in node.children[2].children], ['C', 'D'])

		self.assertRaises(AssertionError, newick.read_internal_node, 'A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;', 0)

	def test_read_leaf_node(self):
		node, curr = newick.read_leaf_node('A:0.1,', 0)
		self.assertEqual(node.name, 'A')
		self.assertEqual(node.distance, 0.1)
		self.assertEqual(curr, 5)

	def test_parse(self):
		tree1 = newick.parse('(cat)dog;')
		tree2 = newick.parse('(dog:42,cat:33);')

		self.assertEqual(tree1.root.name, 'dog')
		self.assertEqual(tree2.root.name, None)
		self.assertEqual(tree2.root.children[0].name, 'dog')
		self.assertEqual(tree2.root.children[0].distance, 42)
