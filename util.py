import os
import networkx as nx
import collections


'''
Input: Binary Mutation Matrix (BMM) or weight Mutation Matrix (WMM). Both are dataframes
Output: Mutation Matrix in a text file
'''
def create_txt_file_mutation_matrix(matrix, output_file_name):
    matrix.to_csv(output_file_name, sep='\t')

'''
Input: alteration matrix (dictionary)
Output: alteration matrix in a file
'''
def create_txt_file_alteration_matrix(matrix, output_file_name):
    output_file = open(output_file_name, 'w')
    for sample in matrix:
        output_file.write(sample + '\t' + '\t'.join(matrix[sample]) + '\n') # removing brackets and commas
    output_file.close()
   
   
'''
Input: Alteration matrix
Output: A simple report about the samples
'''
def get_sample_report(am, output_file_name):
    output_file = open(output_file_name, 'w')
    
    num_samples = len(am)
    output_file.write("Number of samples: {}\n".format(num_samples))
    for sample in am:
        output_file.write("{}: {} mutations\n".format(sample, len(am[sample])))
        
    output_file.close()
    return 0

'''
Input: Binary Mutation matrix
Output: Output: A simple report about the genes
'''
def get_gene_report(bmm, output_file_name):
    output_file = open(output_file_name, 'w')
    
    genes = list(bmm)
    output_file.write("Number of genes: {}\n".format(len(genes)))
    for gene in genes:
        q = bmm[gene].sum()
        f = round((q / len(list(bmm.index))) * 100, 2)
        output_file.write("{}: mutated in {} ({}%) samples\n".format(gene, q, f))
    
    output_file.close()
    return 0

'''
Input: Related Gene Network
Output: Information about the RGN
'''
def get_network_report(rgn, output_file_name):
    print("Number of nodes (genes) in the network: {}".format(len(list(rgn.nodes))))
    print("Number of edges (interactions) in the network: {}".format(len(list(rgn.edges))))

   
'''
Input: output folder and the name of each mutation matrix"
Output: true if all files exists
'''
def files_exists(folder, wmm, bmm, am):
    if os.path.isfile(folder + "/" + wmm) and os.path.isfile(folder + "/" + bmm) and os.path.isfile(folder + "/" + am):
        return True
    return False
    
