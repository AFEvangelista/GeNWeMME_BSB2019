import sys
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it


'''
Input: a gene network (netowrkX) and a size k
Output: a list of connected k-size components (networkx type)
'''
def find_connected_components(rgn, k):
    connected_components = []

    combinations = it.combinations(rgn, k) # generating combinations from the nodes of the graph
    for comb in combinations:
        component = rgn.subgraph(comb) # each nodes combinations, extract a component (subgraph)
        if nx.is_connected(component): # check if the component is connected
            connected_components.append(component)
        
    return connected_components

'''
Input: Weight Mutation Matrix (WMM)
Output: the sum of the weights of the genes in WMM

It is the Relevance Function I(c) of the method
'''
def mutation_relevance_function(wmm, genes): # I(c)
    s = 0
    for gene in genes:
        s = s + wmm[gene].sum()
    return s

'''
Input: Related Gene Network (RGN) and a list of genes
Output: The sum of edges of the genes in the RGN

It is the Related Gene Function R(c) of the method
'''
def genes_relation_function(genes): # R(c)
    edge_weights = nx.get_edge_attributes(genes,'weight')
    s = 0
    for edge in edge_weights:
        s = s + edge_weights[edge]
    return s

'''
Input: A dictionary, with the alteration_matrix and a list of genes
Output: The dendrix weight function score, for the set of genes passed in the parameter

It is the Mutual Exclusivity Weight Function W(c) of the method
'''
def dendrix_weight_function(alteration_matrix, genes): # W(c)
    cont_samples = 0  # samples in which g is mutated
    sample_set = set()  # samples with at least one gene mutated
    for g in genes:
        for sample in alteration_matrix:
            if g in alteration_matrix[sample]:
                cont_samples = cont_samples + 1
                sample_set.add(sample)
    coverage_overlap = cont_samples - len(sample_set)
    weight = len(sample_set) - coverage_overlap
    if weight < 0:
        weight = 0
    return weight

'''
Input: a component and the alteration matrix
Output: number and frequency of samples that the component is all mutated
'''
def get_mutated_samples_component(component, am):
    count = 0
    
    for sample in am:
        if set(component) <= set(am[sample]): # if component is a subset of mutated genes of the sample
            count = count + 1
    frequency = round((count / len(am))*100, 2)
    return count, frequency


'''
Input: wmm, alteration matrix, components and thresholds
Output: the list of components sorted by the ranking function
'''
def rank_function(wmm, am, connected_components, alpha, beta, delta):
    result_df = pd.DataFrame(columns=["number_and_frequency_of_mutated_samples",
                                      "mutation_relevance_raw", "mutation_relevance",
                                      "genes_relation_raw", "genes_relation",
                                      "me_weight_raw", "me_weight",
                                      "score"])
    result_df.set_index = "component"
    for component in connected_components:
        component_in_a_list = list(component)
        
        mutation_relevance_raw = mutation_relevance_function(wmm, component_in_a_list) # considering the type of mutation
        genes_relation_raw = genes_relation_function(component) # considering the weights in the RGN
        me_weight_raw = dendrix_weight_function(am, component_in_a_list) # weight function from Dendrix method
        
        c = ', '.join(component_in_a_list)
        
        count, frequency = get_mutated_samples_component(component_in_a_list, am)
        formated_str = str(count) + " (" + str(frequency) + "%)"
        result_df.loc[c] = [formated_str, mutation_relevance_raw, 0, genes_relation_raw, 0, me_weight_raw, 0, 0]

        
    mutation_relevance_max = result_df["mutation_relevance_raw"].max()
    genes_relation_max = result_df["genes_relation_raw"].max()
    me_weight_max = result_df["me_weight_raw"].max()
    
    for component in result_df.index:
        result_df.at[component, "mutation_relevance"] = result_df.at[component, "mutation_relevance_raw"] / mutation_relevance_max
        result_df.at[component, "genes_relation"] = result_df.at[component, "genes_relation_raw"] / genes_relation_max
        result_df.at[component, "me_weight"] = result_df.at[component, "me_weight_raw"] / me_weight_max
        result_df.at[component, "score"] = alpha * result_df.at[component, "mutation_relevance"] + beta * result_df.at[component, "genes_relation"] + delta * result_df.at[component, "me_weight"]

    result_df = result_df.sort_values(by='score', ascending=False)
   
    return result_df
