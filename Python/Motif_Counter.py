import networkx as nx
import matplotlib.pyplot as plt
import itertools
from collections import *
import networkxgmml as xgmml
import numpy
import pickle
import math


LOG = True

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
	if LOG:
		print "\tEntering randomization..."
	m_count = dict(zip(motifs.keys(), [[] for m in motifs.keys()]))
	m_avg = {}
	m_deviation = {}

	for i in range(trials):
		if LOG:
			print "\t\tCalculating random trial " + str(i) + "of" + str(trials)
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
	count = 0

	if LOG:
		print "\t\t Counting triads"
	for nodes in sub_graphs:
		sub_graph = network.subgraph(nodes)
		possibles = {}
		for motif_name in motifs.keys():
			count += 1
			if LOG and count % 5000 == 0:
				print "\t\t\t" + str(count/float(len(sub_graphs) * 13))
			if nx.faster_could_be_isomorphic(sub_graph, motifs[motif_name]):
				possibles[motif_name] = sub_graph

		# Only one possibility? Add it
		if len(possibles.keys()) == 1:
			m_count[possibles.keys()[0]] += 1
		else: 
			for possible_name in possibles:
				if nx.is_isomorphic(sub_graph, motifs[possible_name]):
					m_count[possible_name] += 1
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
			triads.append(sub_graph)
	return triads

def generate_triads_2(network):
	if LOG:
		print "\t\t\t Generating Triads"
	# Treat the graph as undirected
	u_network = nx.Graph(network)
	visited = set()
	triads = []
	edges = u_network.edges()
	count = 0
	for node in u_network.nodes():
		count += 1
		if LOG and count % 20 == 0:
			print "\t\t\t\t" + str(count/float(len(u_network.nodes())))
		neighbors = u_network.neighbors(node)
		if len(neighbors) < 2: 
			continue
		pairs = itertools.combinations(neighbors, 2)
		for pair in pairs:
			nodes = sorted([node, pair[0], pair[1]])
			if str(nodes) not in visited:
				visited.add(str(nodes))
				triads.append(nodes)
	return triads

def config_network(network):
	in_seq = network.in_degree().values()
	out_seq = network.out_degree().values()
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

def show_convergence(network, trials = 5):
	# def motif_profile_random(network, trials = 100, motifs = get_motif_dic()):
	# """
	# Calculates a motif profile from n trials of randomized config network based on the degree seq of 
	# 'network'.

	# Returns a tuple: ([motif_name, avg_frequency], [motif_name, standard_deviation])
	profiles = []
	for i in range(trials):
		if LOG:
			print "\tCalculating randomized convergence" + str(i) + "of" + str(trials)
		newProfile = motif_profile_random(network, 1)[0]
		if i != 0:
			for motif in newProfile.keys():
				newProfile[motif] += profiles[-1][motif]
		profiles.append(newProfile)

	# Generate average of results
	for i in range(trials):
		profile = profiles[i]
		for key in profile.keys():
			profile[key] *= 1/float(i+1)

	return profiles

def get_z_scores(network, trials = 100):
	m_real = count_motifs(network, get_motif_dic())
	m_rand = motif_profile_random(network, trials)
	z_scores = z_score(m_real, m_rand[0], m_rand[1])
	return z_scores
	

if __name__ == "__main__":
	print "Beginning main"
	eng_words = import_network("../datasets/eng_words.txt")
	jap_words = import_network("../datasets/japanese_words.txt")
	fr_words = import_network("../datasets/french_words.txt")
	sp_words = import_network("../datasets/spanish_words.txt")

	yeast = import_network("../datasets/yeast.txt")
	ecoli = import_network("../datasets/ecoli.txt")

	social = import_network("../datasets/prisoninter.txt")

	p1 = import_network("../datasets/protein_1.txt")
	p2 = import_network("../datasets/protein_2.txt")
	p3 = import_network("../datasets/protein_3.txt")
	
	e1 = import_network("../datasets/s208.txt")
	e2 = import_network("../datasets/s420.txt")
	e3 = import_network("../datasets/s838.txt")

	networks = {"eng_words" : eng_words, "jap_words" : jap_words, "fr_words" : fr_words, "sp_words" : sp_words, "yeast" : yeast, "ecoli" : ecoli, "social" : social,"p1" : p1, "p2" : p2, "p3" : p3, "e1" : e1, "e2" : e2, "e3" : e3}
	

	def plot_conv():
		s4_conv = pickle.load(open("e1_conv.txt"))
		plot = plt.plot(s4_conv)


	# Show how long the config takes to stabilize
	for name in {"eng_words": eng_words, "yeast": yeast, "social": social, "p1": p1, "e1": e1}.keys():
		if LOG:
			print "---------"
			print "Testing convergence for " + str(name)
		network = networks[name]
		if name == "eng_words":
			y_conv = show_convergence(network, 3)
		else:
			y_conv = show_convergence(network, 100)
		s4_conv = [c["s4"] for c in y_conv]
		pickle.dump(s4_conv, open(name + "_conv","wb"))	
	quit()

	# TODO: Plot the convergence profiles of each type of thing
		# Test the z_score calculation - decent?
		# Setup z_score save for each network




	z = get_z_scores(network, 10)
	ecoli_z = get_z_scores(ecoli,10)
	print ecoli_z
	ecoli_profile = significance_profile(ecoli_z)
	quit()




	y_conv = show_convergence(ecoli, 100)
	s4_conv = [c["s4"] for c in y_conv]
	pickle.dump(s4_conv, open("s4_conv_ecoli", "wb"))	

	y_conv = show_convergence(yeast, 100)
	s4_conv = [c["s4"] for c in y_conv]
	pickle.dump(s4_conv, open("s4_conv_yeast", "wb"))

	y_conv = show_convergence(words, 100)
	s4_conv = [c["s4"] for c in y_conv]
	pickle.dump(s4_conv, open("s4_conv_words", "wb"))

	















	# convergence = show_convergence(network2, 1000);
	# # s2_conv = [c["s2"] for c in convergence]
	# pickle.dump(s1_conv, open("conv", "wb"))

	# s1_conv = pickle.load(open("conv","rb"))
	# print s1_conv
	# # print s1_conv
	# plt.show()



	quit()	
	# plot_profile(profile)
	# print profile
