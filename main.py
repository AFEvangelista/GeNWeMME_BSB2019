import datetime
import os
import sys
import configparser
import pandas as pd
import networkx as nx
import step1
import step2
import step4
import util
import logging
import reports



def main():
    input_parameters_file = sys.argv[1]
    cp = configparser.ConfigParser() 
    cp.read(input_parameters_file)
    
    cp_input_type = cp["INPUT_TYPE"]
    input_type = int(cp_input_type["TYPE"])
    
    cp_input = cp["INPUT"]
    
    if input_type == 1:
        input_maf_file_name = cp_input["MAF_FILE_NAME"]
        input_gene_network_file_name = cp_input["GENE_NETWORK_FILE_NAME"]
    elif input_type == 2:
        input_bmm_file_name = cp_input["BMM_FILE_NAME"]
        input_wmm_file_name = cp_input["WMM_FILE_NAME"]
        input_rgn_file_name = cp_input["RGN_FILE_NAME"]
    elif input_type == 3:
        input_maf_file_name = cp_input["MAF_FILE_NAME"]
        input_wgn_file_name = cp_input["WEIGHTED_GENE_NETWORK_FILE_NAME"]
    else: # input_type == 4
        input_bmm_file_name = cp_input["BMM_FILE_NAME"]
        input_wmm_file_name = cp_input["WMM_FILE_NAME"]
        input_wgn_file_name = cp_input["WEIGHTED_GENE_NETWORK_FILE_NAME"]
    
    
    pre_processing_threshold_phi = float(cp_input["PRE_PROCESSING_THRESHOLD_PHI"])
    rgn_threshold_gamma = float(cp_input["RGN_THRESHOLD_GAMMA"])
    component_size_k = int(cp_input["COMPONENT_SIZE_K"])
    mutation_relevance_alpha = float(cp_input["MUTATION_RELEVANCE_ALPHA"])
    related_genes_beta = float(cp_input["RELATED_GENES_BETA"])
    mutual_exclusivity_delta = float(cp_input["MUTUAL_EXCLUSIVITY_DELTA"])
    
    cp_mutation_weights = cp["VARIANT_CLASSIFICATION_WEIGHTS"]
    mutation_weights = {}
    mutation_weights["Nonsense_Mutation"] = float(cp_mutation_weights["NONSENSE"])
    mutation_weights["Missense_Mutation"] = float(cp_mutation_weights["MISSENSE"])
    mutation_weights["Splice_Site"] = float(cp_mutation_weights["SPLICE"])
    mutation_weights["Frame_Shift_Del"] = float(cp_mutation_weights["FRAMESHIFT_DEL"])
    mutation_weights["Frame_Shift_Ins"] = float(cp_mutation_weights["FRAMESHIFT_INS"])
    mutation_weights["In_Frame_Del"] = float(cp_mutation_weights["INFRAME_DEL"])
    mutation_weights["In_Frame_Ins"] = float(cp_mutation_weights["INFRAME_INS"])
    mutation_weights["3'UTR"] = float(cp_mutation_weights["_3UTR"])
    mutation_weights["5'UTR"] = float(cp_mutation_weights["_5UTR"])
    mutation_weights["Nonstop_Mutation"] = float(cp_mutation_weights["NONSTOP_MUTATION"])
    mutation_weights["Translation_Start_Site"] = float(cp_mutation_weights["TRANSLATION_START_SITE"])
    
    cp_output = cp["OUTPUT"]
    output_folder = cp_output["OUTPUT_FOLDER"]
    output_wmm_file_name = cp_output["WMM_FILE"]
    output_bmm_file_name = cp_output["BMM_FILE"]
    output_am_file_name = cp_output["AM_FILE"]
    output_rgn_file_name = cp_output["RGN_FILE"]
    output_components_file_name = cp_output["COMPONENTS_FILE"]
    output_ranking_file_name = cp_output["RANKING_FILE"]
    output_sample_report_file_name = cp_output["SAMPLE_REPORT_FILE"]
    output_gene_report_file_name = cp_output["GENE_REPORT_FILE"]
    output_network_report_file_name = cp_output["NETWORK_REPORT_FILE"]
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', filename=output_folder + '/' + 'log.txt')
    logging.info('******** Starting ********')
    
    if input_type == 1 or input_type == 3:
        # Step 1
        maf = pd.read_csv(input_maf_file_name, sep="\t", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
        step1.filter_genes_min_freq(maf, pre_processing_threshold_phi)
        wmm = step1.get_weighted_mutation_matrix(maf, mutation_weights)
        bmm = step1.get_binary_mutation_matrix(maf)
        am = step1.get_alteration_matrix(bmm)
        logging.info('Input Type 1 or 3 - End of Step 1')
    else: # input_type == 2 or input_type == 4:
        # Step 1 and 2 - Using files generated before
        wmm = step1.get_mutation_matrix_from_file(input_wmm_file_name)
        bmm = step1.get_mutation_matrix_from_file(input_bmm_file_name)
        am = step1.get_alteration_matrix(bmm)
        logging.info('Input Type 2 or 4 - End of Step 1')
    
    # generating mutation matrices file
    util.create_txt_file_mutation_matrix(wmm, output_folder + "/" + output_wmm_file_name)
    util.create_txt_file_mutation_matrix(bmm, output_folder + "/" + output_bmm_file_name)
    util.create_txt_file_alteration_matrix(am, output_folder + "/" + output_am_file_name)
    # generating reports about samples and genes
    util.get_sample_report(am, output_folder + "/" + output_sample_report_file_name)
    util.get_gene_report(bmm, output_folder + "/" + output_gene_report_file_name)
    
    genes_in_analysis = step1.get_genes_from_bmm(bmm)
    if input_type == 1:
        # Step 2
        gene_network_file = open(input_gene_network_file_name, 'rb')
        gene_network_nx = step2.read_gene_network_nx(gene_network_file)
        logging.info('Input Type 1 - End of Step 2 - Read gene network')
        logging.info('Input Type 1 - End of Step 2 - Getting genes in the analysis')
        rgn = step2.create_related_gene_network(gene_network_nx, rgn_threshold_gamma, genes_in_analysis)    
        gene_network_file.close()
        nx.write_weighted_edgelist(rgn, output_folder + "/" + output_rgn_file_name, delimiter='\t')
        logging.info('End of Step 2')
    elif input_type == 2:
        # Step 2 - Using file generated before
        rgn = nx.read_weighted_edgelist(input_rgn_file_name, delimiter='\t')
        logging.info('Input Type 2 - End of Step 2')
    else: #input_type == 3 or input_type == 4:
        wgn = nx.read_weighted_edgelist(input_wgn_file_name, delimiter='\t')
        rgn = step2.create_related_gene_network_from_weighted(wgn, rgn_threshold_gamma, genes_in_analysis)
        nx.write_weighted_edgelist(rgn, output_folder + "/" + output_rgn_file_name, delimiter='\t')
        logging.info('Input Type 3 or 4 - End of Step 2')
    
    # generating report about the Related Gene Network (RGN)
    #util.get_network_report(rgn, output_folder + "/" + network_report_file_name)
    #diff = set(list(bmm)) - set(list(rgn.nodes))
    #print(diff)
    
    #Step 4
    connected_components = step4.find_connected_components(rgn, component_size_k)
    
    components_file = open(output_folder + "/" + output_components_file_name, "w")
    for c in connected_components:
        components_file.write("\t".join(c) + "\n")
    components_file.close()
    
    logging.info('Step 4 - End of finding components')
    
    ranking = step4.rank_function(wmm, am, connected_components, mutation_relevance_alpha, related_genes_beta, mutual_exclusivity_delta)
    
    ranking_file = open(output_folder + "/" + output_ranking_file_name, "w")
    
    ranking_file.write("component\tnumber and frequency of mutated samples\ttype of mutation\tgenes relation\tmutual exclusivity\tscore\n")
    for component in ranking.index:
        ranking_file.write("[{}]\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(component, ranking.at[component, "number_and_frequency_of_mutated_samples"], ranking.at[component, "mutation_relevance"], ranking.at[component, "genes_relation"], ranking.at[component, "me_weight"], ranking.at[component, "score"]))
        
    ranking_file.close()
    logging.info('End of Step 4')
    #print(ranking[["mutation_relevance_raw", "mutation_relevance", "genes_relation_raw", "genes_relation", "me_weight_raw", "me_weight", "score"]])


# python main.py input_parameters.in 
main()
