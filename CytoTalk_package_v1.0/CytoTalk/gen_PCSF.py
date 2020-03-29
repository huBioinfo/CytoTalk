# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 12:09:32 2018

@author: Administrator
"""
import numpy
import pcst_fast
import datetime

Edge = numpy.loadtxt("../../IntegratedNet_edge.txt", dtype = 'int64', delimiter = "\t")
Cost = numpy.loadtxt("IntegratedNet_edgeCost.txt", dtype = 'float')
Prize = numpy.loadtxt("../IntegratedNet_nodePrize.txt", dtype = 'float')
rootNode = numpy.loadtxt("../../RootNode.txt", dtype = 'int64')

# vertices, edges = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
# root: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
# pruning: The standard GW pruning method is 'gw', which is also the default.
# verbosity_level: an integer indicating how much debug output the function should produce. Note!!!should include this parameter.
PCSF_node, PCSF_edge =pcst_fast.pcst_fast(Edge, Prize, Cost, rootNode, 1, 'strong', 0)

numpy.savetxt("PCSF_Node.txt", PCSF_node, fmt='%d')  #-correct
numpy.savetxt("PCSF_Edge.txt", PCSF_edge, fmt='%d')


