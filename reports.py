import networkx as nx
import pandas as pd




'''

'''
def sample_report(am):
    num_samples = len(am)
    print("Number of samples: {}".format(num_samples))
    for sample in am:
        print("{}: {} mutations".format(sample, len(am[sample])))
    return 0
        





