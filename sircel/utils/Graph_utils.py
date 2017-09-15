"""
Akshay Tambe
Pachter and Doudna groups
UC Berkeley

Dropseq_utils.py
A few useful objects for working with single-cell barcode
data as de Bruijn graphs
"""
import sys
import numpy as np

class Edge:
	"""
	An edge in a de Bruijn graph
	Attributes
		node (string)
		neighbor (string)
		weight (float)
	To be a valid edge, all the nucleotides in node and neighbor must 
		be the same, except for the first nt of the node / 
		last nt of the edge
	"""
	def __init__(self, _node, _neighbor, _weight):
		self.node = _node
		self.neighbor = _neighbor
		self.weight = _weight
		assert (self.is_valid_edge()), \
			print('Not a valid edge: %s\t%s' % (self.node, self.neighbor))
	
	def get_name(self):
		"""
		Returns a tuple (node, neighbor) for this edge
		"""
		return(self.node, self.neighbor)
	
	def get_weight(self):
		return self.weight
	
	def get_sequence(self):
		return(self.node[0] + self.neighbor)
		
	def is_valid_edge(self):
		return(
			(self.node[1:] == self.neighbor[0:-1]) and \
			(self.weight >= 0.0))
	
	def is_self_edge(self):
		""""
		Returns true if node and neigobor are the same
		"""
		return(self.node == self.neighbor)

class Path:
	"""
	A path through a de Bruijn graph
	Attributes
		edges (list): a list of Edge objects
	
	To be a valid path, edges[i].node must be equal to edges[i-1].neighbor
		for all entries in self.edges
	"""
	def __init__(self, _edges):
		self.edges = _edges
		assert (self.is_valid_path()), \
			print('Not a valid path\n%s' % self.get_nodes_ordered())
	
	def get_weight(self):
		"""
		Returns float
			the weight of the lowest-weight edge
		"""
		return(self.get_lowest_weight_edge().get_weight())
	
	def get_cycle_weight(self):
		"""
		Returns float
			0 if the path is not a cycle
			the path weight otherwise
		"""
		if(not self.is_cycle()):
			return 0.0
		else:
			return self.get_weight()
	
	def get_lowest_weight_edge(self):
		return(
			min(self.edges,
			key = lambda edge: edge.get_weight()))
	
	def get_sequence(self):
		edges = self.edges
		seq = ''
		seq = (edges[0]).node
		for edge in edges[1:]:
			seq += (edge.node)[-1]
		return seq
	
	def get_sequence_circular(self):
		"""
		Gets the sequence for a cycle by starting at a node beginning with '$'
		"""
		seq = self.get_sequence()
		if(not self.is_cycle()):
			return seq
		
		start_ind = 0
		if('$' in seq):
			start_ind = seq.index('$') + 1
		seq_len = self.get_length() - 1
		seq = seq[start_ind:] + seq[0:start_ind]
		return(seq[0:seq_len])
	
	def get_start_node(self):
		return self.edges[0].node
	
	def get_end_node(self):
		return self.edges[-1].neighbor
	
	def get_nodes_ordered(self):
		nodes = [edge.node for edge in self.edges]
		if(not self.is_cycle()):
			nodes.append(self.edges[-1].neighbor)
		return sorted(nodes)
			
	def get_length(self):
		return len(self.edges)
	
	def is_cycle(self):
		return (self.get_start_node() == self.get_end_node())
	
	def is_valid_path(self):
		if(len(self.edges) == 1):
			return True
		for i in range(0, len(self.edges) - 1):
			neighbor = self.edges[i].neighbor
			next_node = self.edges[i+1].node
			if(neighbor != next_node):
				print(neighbor, next_node)
				return False
		return True
		
	def is_possible_cycle(self, cycle_length):
		"""
		Returns whether this path might be part of a cycle of fixed length
			If path_length > (cycle_length - kmer_size) there is some substring
			that is shared between the last position(s) of the path's end node and
			the first position(s) of its start node
		"""		
		kmer_size = len(self.get_start_node())
		if(self.get_length() < (cycle_length - kmer_size)):
			return True
		if(self.is_cycle()):
			return True
		
		overlap = self.get_length() - (cycle_length - kmer_size)
		if((self.get_start_node()[0:overlap]) == 
				(self.get_end_node()[kmer_size - overlap:])):
			return True
		return False

class Graph:
	def __init__(self, _edges_lst):
		#create dict of edge objects
		self.edges = {}
		self.self_edges = {}
		for edge in _edges_lst:			
			if(edge.is_self_edge()):
				key = edge.node
				self.self_edges[key] = edge
			else:
				key = edge.get_name()
				self.edges[key] = edge
		
	def get_num_edges(self):
		return(
			len(self.edges.values()) + \
			len(self.self_edges.values()))
	
	def get_total_weight(self):
		total_weight = sum( \
				[edge.get_weight() for edge in self.edges.values()]) + \
			sum( \
				[edge.get_weight() for edge in self.self_edges.values()])
		return total_weight
		
	def get_edges_sorted(self):
		edges_list = sorted(
			self.edges.values(),
			key = lambda edge: edge.get_weight(),
			reverse=True)
		return edges_list
	
	def get_outgoing_edges(self, node):
		alphabet = ['A', 'C', 'G', 'T', '$']
		outgoing_edges = []
		for a in alphabet:
			neighbor = ''.join([node[1:] + a])
			key = (node, neighbor)
			if(key in self.edges):
				outgoing_edges.append(self.edges[key])
		return outgoing_edges
	
	def get_outgoing_edges_sorted(self, node):
		return sorted(
			self.get_outgoing_edges(node), 
			key = lambda edge: edge.get_weight(), 
			reverse = True)#sorted by edge weight (descending)]

	def find_cyclic_path(self, current_path, expected_path_length):
		"""
		Args
			current_path (Path):
				a Path (of variable length)
			expected_path_length (int):
				length of cycles to look for (# edges in cyclic path)
		Returns 
			Cycle (boolean): whether or not a cycle was found
			Path (Path): a path object representing the best path

		Recursively finds high-weight cyclic paths of fixed length
			Initialize by passing a path of length 1 (a single edge)
			Recursion terminates if:
				If current path length is larger than expected path length
				If current path cannot be a cycle of expected length
				If current path is a cycle of expected length
			Otherwise recursively:
				Return highest weight cyclic path from each outgoing node
		
		This recursion has a fixed upper bound for branching. For any starting edge, 
			The total number of possible paths is given by: 
				4**(expected_path_length - kmer_size)
		"""
		#check if path is too long (fail condition)		
		if(current_path.get_length() > expected_path_length + 1):
			return(False, current_path)	
		#check that this path is possibly part of a cycle (fail condition)
		if(not current_path.is_possible_cycle(expected_path_length)):
			return(False, current_path)
		#check if path is a cycle of expected length (success condition)
		if(current_path.is_cycle() and 
			np.fabs(current_path.get_length() - expected_path_length) <= 1):
			return(True, current_path)
		#check if path is a smaller cycle, and if so check for self edges
		if(current_path.is_cycle() and 
			current_path.get_length() < expected_path_length):
			
			(possible_self_edges,
				updated_path) = self.check_possible_self_edges(
				current_path, 
				expected_path_length)
			if(not possible_self_edges):
				return(False, current_path)
			else:
				return(True, updated_path)

		#else return best cycle of all outgoing edges (continue recursion)
		current_node = current_path.get_end_node()
		outgoing_edges = self.get_outgoing_edges_sorted(current_node)
		#outgoing_edges is sorted in descending order by weight
		if(len(outgoing_edges) == 0):
			return (False, current_path)
		best_path = current_path
		for (i, outgoing_edge) in enumerate(outgoing_edges):
			downstream_path = Path(current_path.edges + [outgoing_edge])
			is_cycle, path = self.find_cyclic_path(
				downstream_path, 
				expected_path_length)
			if(path.get_cycle_weight() > best_path.get_cycle_weight()):
				best_path = path
			if(i < len(outgoing_edges) - 1 and \
				best_path.is_cycle() and \
				best_path.get_cycle_weight() > outgoing_edges[i+1].get_weight()):
				break
		return(best_path.is_cycle(), best_path)
	
	def check_possible_self_edges(self, current_path, expected_path_length):
		possible_self_edges = []
		for edge in current_path.edges:
			node = edge.node
			if(node in self.self_edges.keys()):
				self_edge = self.self_edges[node]
				possible_self_edges.append(self_edge)
		if(len(possible_self_edges) == 0):
			return (False, current_path)
		#else
		best_self_edge = max(
			possible_self_edges, 
			key = lambda edge: edge.get_weight())
		num_self_edges_needed = expected_path_length - current_path.get_length()
		
		new_edges = []
		min_edge_weight = current_path.get_cycle_weight()
		for i in range(current_path.get_length()):
			this_edge = current_path.edges[i]
			if(this_edge.node == best_self_edge.node):
				for j in range(num_self_edges_needed):
					new_edge = Edge(this_edge.node, this_edge.node, min_edge_weight)
					new_edges.append(new_edge)
			new_edges.append(this_edge)
		return (True, Path(new_edges))
	
	def find_all_cyclic_paths(self, start_node, start_neighbor, expected_path_length):
		key = (start_node, start_neighbor)
		if(key not in self.edges):
			raise StopIteration
		
		start_edge = self.edges[key]
		initial_path = Path([start_edge])
		while(True):
			cycle, path = self.find_cyclic_path(initial_path, expected_path_length)			
			if(path.get_length() != expected_path_length):
				raise StopIteration
			elif(cycle):
				yield path
			#decrement edges in graph by cycle weight
			cycle_weight = path.get_cycle_weight()
			for edge in path.edges:
				key = edge.get_name()
				if(key in self.edges.keys()):
					self.edges[key].weight -= cycle_weight

			
		
		

