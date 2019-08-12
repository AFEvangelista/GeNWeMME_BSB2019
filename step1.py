import sys
import pandas as pd

# constants
NONSENSE = "Nonsense_Mutation"
MISSENSE = "Missense_Mutation"
SPLICE ="Splice_Site"
FRAMESHIFT_DEL = "Frame_Shift_Del"
FRAMESHIFT_INS = "Frame_Shift_Ins"
INFRAME_DEL = "In_Frame_Del"
INFRAME_INS = "In_Frame_Ins"
_3UTR = "3'UTR"
_5UTR = "5'UTR"
NONSTOP_MUTATION = "Nonstop_Mutation"
TRANSLATION_START_SITE = "Translation_Start_Site"

'''
Input: Pandas DataFrame with the MAF file information
OutPut: Binary Mutation Matrix (BMM), in a pandas DataFrame

Binary Mutation matrix format:
Rows are patients ((samples) and columns are genes.
1: gene is mutated in patient; 0: gene is not mutated in patient.

Example:
Hugo_Sumbol g1  g2  g3  ... gm
p1  0   0   1   0
p2  1   0   0   1
p3  1   0   1   0
...
pn  1   0   0   0

'''
def get_binary_mutation_matrix(maf):
    bmm = pd.crosstab(maf.Tumor_Sample_Barcode, maf.Hugo_Symbol).clip_upper(1)
    return bmm


'''
Input: Pandas DataFrame with the MAF file information
Output: Weighted Mutation Matrix (WMM), in a pandas DataFrame

TODO: insert scores for known cancer related gene (likely MSigDB)
'''
def get_weighted_mutation_matrix(maf, mutation_weights):
    patients = list(set(maf.Tumor_Sample_Barcode))
    genes = list(set(maf.Hugo_Symbol))
    mutation_types = list(set(maf.Variant_Classification))
    #print("patients: {}".format(patients))
    #print("genes: {}".format(genes))
    #print("mutation types: {}".format(mutation_types))
    
    # cross tab to get the number of mutations of each pair gene-patient
    mut = pd.crosstab(maf.Hugo_Symbol, [maf.Tumor_Sample_Barcode, maf.Variant_Classification], values=maf.Variant_Classification, aggfunc=len, dropna=False)

    wmm = pd.DataFrame(index=patients, columns=genes) # Creating Weighted Mutation Matrix to keep the scores.
    wmm = wmm.fillna(0.0) # fill with 0s rather than NaNs
    
    # calcuting the scores for each pait gene-patient: s(p_i, g_j)
    for p in patients:
        mut[p] = mut[p].fillna(0) # fill with 0s rather NaNs
        for g in genes:
            nonsense_count = mut[p][NONSENSE][g] if NONSENSE in mutation_types else 0
            nonsense_weight = nonsense_count * mutation_weights[NONSENSE]
            
            missense_count = mut[p][MISSENSE][g] if MISSENSE in mutation_types else 0
            missense_weight = missense_count * mutation_weights[MISSENSE]
            
            splice_count = mut[p][SPLICE][g] if SPLICE in mutation_types else 0
            splice_weight = splice_count * mutation_weights[SPLICE]
            
            frame_shift_ins_count = mut[p][FRAMESHIFT_INS][g] if FRAMESHIFT_INS in mutation_types else 0
            frame_shift_ins_weight = frame_shift_ins_count * mutation_weights[FRAMESHIFT_INS]
            
            frame_shift_del_count = mut[p][FRAMESHIFT_DEL][g] if FRAMESHIFT_DEL in mutation_types else 0
            frame_shift_del_weight = frame_shift_del_count * mutation_weights[FRAMESHIFT_DEL]
            
            inframe_ins_count = mut[p][INFRAME_INS][g] if INFRAME_INS in mutation_types else 0
            inframe_ins_weight = inframe_ins_count * mutation_weights[INFRAME_INS]
            
            inframe_del_count = mut[p][INFRAME_DEL][g] if INFRAME_DEL in mutation_types else 0
            inframe_del_weight = inframe_del_count * mutation_weights[INFRAME_DEL]
            
            _3utr_count = mut[p][_3UTR][g] if _3UTR in mutation_types else 0
            _3utr_weight = _3utr_count * mutation_weights[_3UTR]
            
            _5utr_count = mut[p][_5UTR][g] if _5UTR in mutation_types else 0
            _5utr_weight = _5utr_count * mutation_weights[_5UTR]
            
            nonstop_mutation_count = mut[p][NONSTOP_MUTATION][g] if NONSTOP_MUTATION in mutation_types else 0
            nonstop_mutation_weight = nonstop_mutation_count * mutation_weights[NONSTOP_MUTATION]
            
            translation_start_site_count = mut[p][TRANSLATION_START_SITE][g] if TRANSLATION_START_SITE in mutation_types else 0
            translation_start_site_weight = translation_start_site_count * mutation_weights[TRANSLATION_START_SITE]
            
            score = 0
            count_total = nonsense_count + missense_count + splice_count + frame_shift_ins_count + frame_shift_del_count + inframe_ins_count + inframe_del_count + _3utr_count + _5utr_count + nonstop_mutation_count + translation_start_site_count
            weight_total = nonsense_weight + missense_weight + splice_weight + frame_shift_ins_weight + frame_shift_del_weight + inframe_ins_weight + inframe_del_weight + _3utr_weight + _5utr_weight + nonstop_mutation_weight + translation_start_site_weight
            if count_total > 0:
                score =  round(weight_total / count_total, 3)
            wmm.at[p, g] = score
            
    return wmm

def get_mutation_matrix_from_file(mm_file_name):
    mm = pd.read_csv(mm_file_name, sep='\t', index_col=0)
    return mm

'''
Input: Binary Mutation Matrix (BMM)
Output: An Alteration Matrix, in a dictionary

Each row contais patient (sample) in the first columns. The follow columns have the genes which are mutated in the sample.
Example:
p1  g1  g4  g13
p2  g7  g16 g53 g104
p3  g10 g12
...
pn  g1  g3
'''
def get_alteration_matrix(bmm):
    matrix = dict()
    sample_list = bmm.index.tolist()
    for sample in sample_list:
        matrix[sample] = list(bmm.columns.values[bmm.loc[sample] == 1])

    return matrix

'''
Input: Pandas DataFrame with the MAF file information and a the minimum frequency of a gene in the samples
Output: The input MAF filtered.

The removal is performed when a gene is mutated in less than a specific percentage of samples
'''
def filter_genes_min_freq(maf, min_freq):
    num_patients = len(set(maf.Tumor_Sample_Barcode))
    # keeping only one pait gene-patient
    maf_no_duplicate = maf.drop_duplicates(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])
    genes = list(maf_no_duplicate.Hugo_Symbol) # list of genes
    
    # dictionary with the number of samples in which a gene is mutated
    genes_count = {gene : genes.count(gene) for gene in genes}
    min_quantity = num_patients * min_freq
    genes_to_remove = [gene for gene in genes_count if genes_count[gene] < min_quantity]
    
    for gene in genes_to_remove:
        maf.drop(maf[maf.Hugo_Symbol == gene].index, inplace=True)
    
    return 1

def get_genes_from_bmm(bmm):
    return list(bmm)













