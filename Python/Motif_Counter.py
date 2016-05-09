import networkx as nx
# import matplotlib.pyplot as plt
import itertools
from collections import *
import networkxgmml as xgmml
import numpy
import math


def import_network(file_loc):
	return nx.read_edgelist(file_loc, create_using= nx.DiGraph(), nodetype = str, data = [('weight', float)])

def get_motif_dic():
	"""
	Return a dictionary of all 13 triads. Rtype: [string, subgraph]
	"""
	motifs = {
		"s0" : nx.DiGraph([(2,0), (2,1)]),
		"s1" : nx.DiGraph([(2,0), (1,2)]),
		"s2" : nx.DiGraph([(2,0), (1,2), (2,1)]),
		"s3" : nx.DiGraph([(1,0), (2,0)]),
		"s4" : nx.DiGraph([(1,0), (2,0), (2,1)]),
		"s5" : nx.DiGraph([(1,0), (2,0), (2,1), (1,2)]),
		"s6" : nx.DiGraph([(0,2), (1,2), (2,1)]),
		"s7" : nx.DiGraph([(0,2), (2,0), (1,2), (2,1)]),
		"s8" : nx.DiGraph([(0,2), (2,1), (1,0)]),
		"s9" : nx.DiGraph([(0,2), (2,0), (2,1), (1,0)]),
		"s10" : nx.DiGraph([(0,2), (2,0), (1,2), (1,0)]),
		"s11" : nx.DiGraph([(0,2), (2,0), (1,2), (2,1), (1,0)]),
		"s12" : nx.DiGraph([(0,2), (2,0), (1,2), (2,1), (1,0), (0,1)])
	}

	return motifs
	

def motif_profile_random(network, trials = 100, motifs = get_motif_dic()):
	"""
	Calculates a motif profile from n trials of randomized config network based on the degree seq of 
	'network'.

	Returns a tuple: ([motif_name, avg_frequency], [motif_name, standard_deviation])
	"""
	m_count = dict(zip(motifs.keys(), [[] for m in motifs.keys()]))
	m_avg = {}
	m_deviation = {}

	for i in range(trials):
		config = config_network(network)
		new_rand = count_motifs(config, motifs)
		for motif in new_rand.keys():
			m_count[motif].append(new_rand[motif])

	for motif in m_count.keys():
		m_avg[motif] = numpy.mean(m_count[motif])
		m_deviation[motif] = numpy.std(m_count[motif])

	return (m_avg, m_deviation)

def z_score(m_real, m_rand, m_dev):
	"""
	Calculates the Z_Score of a motif profile of a real network, 
	and a motif profile of a random network.
	"""
	z_scores = dict(zip(m_real.keys(), [0 for m in m_real.keys()]))
	for motif in m_real.keys():
		z = (m_real[motif] - m_rand[motif])/(float(m_dev[motif]))
		z_scores[motif] = z
	return z_scores

def count_motifs(network, motifs):
	"""
	Returns a dictionary of motif counts, keyed on names and valued on occurances. [string, int]
	"""
	m_count = dict(zip(motifs.keys(), [0 for x in motifs.keys()]))
	sub_graphs = generate_triads_2(network)

	for sub_graph in sub_graphs:
		for motif_name in motifs.keys():
			if nx.is_isomorphic(sub_graph, motifs[motif_name]):
				m_count[motif_name] += 1
				break

	return m_count

def generate_triads_1(network):
	"""
	Generate triads by making all node combinations, filtering < 2 edges. 
	Not suitable for large graphs. 

	"""
	nodes = network.nodes()
	node_combos = list(itertools.combinations(nodes, 3))
	triads = []

	for combo in node_combos:
		sub_graph = network.subgraph(combo)
		if nx.is_weakly_connected(sub_graph):
			# print combo
			# print sub_graph.edges()
			triads.append(sub_graph)
	# print len(triads)
	return triads

def generate_triads_2(network):
	# Treat the graph as undirected
	u_network = nx.Graph(network)
	visited = set()
	triads = []
	edges = u_network.edges()

	for node in u_network.nodes():
		neighbors = u_network.neighbors(node)
		if len(neighbors) < 2: 
			continue
		pairs = itertools.combinations(neighbors, 2)
		for pair in pairs:
			nodes = sorted([node, pair[0], pair[1]])
			if str(nodes) not in visited:
				visited.add(str(nodes))
				triads.append(network.subgraph(nodes))
	return triads

def config_network(network):
	in_seq = sorted(network.in_degree().values())
	out_seq = sorted(network.out_degree().values())
	config = nx.directed_configuration_model(in_seq, out_seq, nx.DiGraph())

	# Remove self loops and parallel edges
	config = nx.DiGraph(config)
	config.remove_edges_from(config.selfloop_edges())
	return config

def significance_profile(z_scores):
	profile = dict(zip(z_scores.keys(), [0 for x in z_scores.keys()]))
	a = math.sqrt(reduce(lambda x,y: x + y**2, z_scores.values()))
	for motif in z_scores.keys():
		profile[motif] = z_scores[motif]/float(a)
	return profile

def plot_profile(profile): 
	vals = [profile["s0"], profile["s1"], profile["s2"], profile["s3"], profile["s4"], profile["s5"], profile["s6"], profile["s7"], profile["s8"], profile["s9"], profile["s10"],  profile["s11"],  profile["s12"]]
	print vals
	plot = plt.plot(vals)
	plt.show()

if __name__ == "__main__":
	print "Beginning main"
	network1 = import_network("word_association_graph_DSF.txt")
	network2 = xgmml.XGMMLReader(open("florentineBiz.xgmml", "r"))
	# t1 = generate_triads_1(network2)
	# print "-------"
	# t2 = generate_triads_2(network2)

	m_real = count_motifs(network1, get_motif_dic())
	quit()
	m_rand = motif_profile_random(network2)
	z_scores = z_score(m_real, m_rand[0], m_rand[1])
	profile = significance_profile(z_scores)
	plot_profile(profile)
	print profile
