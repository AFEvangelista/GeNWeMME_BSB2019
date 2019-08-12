import sys
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


'''
Input: the gene network file (reading using rb)
Output: a networkX graph
'''
def read_gene_network_nx(gene_network_file):
    mx = nx.read_edgelist(gene_network_file, delimiter='\t')
    return mx

'''
Input: a gene network (netowrkX), a threshold and a a list of genes considered in the analisys
Output: a Related Gene Network (RGN)
'''
def create_related_gene_network(gene_network_nx, threshold, genes):
    rgn = nx.Graph()
    
    for g1, g2 in gene_network_nx.edges:
        if g1 in genes and g2 in genes:
            rgn.add_edge(g1, g2, weight=1.0)
    
    similarity = nx.jaccard_coefficient(gene_network_nx)
    for g1, g2, coefficient in similarity:
        #print('({}, {}) -> {:.4f}'.format(g1, g2, coefficient))
        # adding in RGN only edges with weight greater than threshold and if genes are in the BMM
        if coefficient > threshold and g1 in genes and g2 in genes: 
            rgn.add_edge(g1, g2, weight=coefficient)

    
    ## routine to print the graph
    #pos = nx.circular_layout(rgn)  # positions for all nodes
    ## nodes
    #nx.draw_networkx_nodes(rgn, pos, node_size=300)
    ## edges
    #nx.draw_networkx_edges(rgn, pos, width=1)
    #labels = nx.get_edge_attributes(rgn,'weight')
    #nx.draw_networkx_edge_labels(rgn, pos, font_size=8,edge_labels=labels)
    ## labels
    #nx.draw_networkx_labels(rgn, pos, font_size=10)

    #plt.axis('off')
    #plt.show()
    
    return rgn

def create_related_gene_network_from_weighted(gene_network_nx, threshold, genes):
    rgn = nx.Graph()

    for g1, g2, edge in gene_network_nx.edges(data=True):
        if edge["weight"] > threshold and g1 in genes and g2 in genes:
            rgn.add_edge(g1, g2, weight=edge["weight"])
   
    return rgn
