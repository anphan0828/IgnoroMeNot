# This script finds ignorome that are highly-associated with highly-annotated genes
# Ranking based on annotation count, Philip Lord information content, PubMed mentions, Jensen score;
# Association based on STRING co-expression

from __future__ import print_function
from __future__ import division
import networkx as nx
import sys
import matplotlib.pyplot as plt
import argparse
import pickle as cp
import math
import numpy as np
from numpy import percentile
import collections
import Bio
from Bio.UniProt import GOA
from Bio.Seq import Seq
from dateutil import parser
import os
import xlsxwriter
import requests
import pandas as pd

verbose = 0
options = ""

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="ignoromenot.py")

    requiredArguments = parser.add_argument_group("Required arguments")
    requiredArguments.add_argument('--ranktable', '-rank', help="Ranking table with 1st column being gene names, 2nd column being Wyatt Clark IC",
                                   required=True)
    requiredArguments.add_argument('--idtable', '-id', help="ID mapping table from STRING",
                                   required=True)
    requiredArguments.add_argument('--stringppi', '-ppi', help="Protein-protein interaction network from STRING",
                                   required=True)
    parser.add_argument('--threshold_top', '-ttop', help="The threshold level above which most-annotated genes will be selected",
                        default=1000)
    parser.add_argument('--threshold_bot', '-tbot', help="The threshold level under which least-annotated genes will be selected",
                        default=5)
    parser.add_argument('--threshold_ppi', '-tppi', help="The threshold value (0.1) above which STRING interactions will be selected",
                        default=500)
    # parser.add_argument('--output', '-odir', help="Writes the final outputs to the directory in this path.")
    # requiredArguments.add_argument('--input', '-i', nargs="+",
    #                                help="The input file path. Please remember the name of the file must start with goa in front of it, with the name of the species following separated by an underscore",
    #                                required=True)
    # mutex_parser_IC_THRESH.add_argument('--info_threshold_Wyatt_Clark_percentile', '-WCTHRESHp',
    #                                     help="Provide the percentile p. All annotations having information content below p will be discarded")
    # mutex_parser_IC_THRESH.add_argument('--info_threshold_Wyatt_Clark', '-WCTHRESH',
    #                                     help="Provide a threshold value t. All annotations having information content below t will be discarded")
    #
    # mutex_parser_IC_THRESH.add_argument('--info_threshold_Phillip_Lord_percentile', '-PLTHRESHp',
    #                                     help="Provide the percentile p. All annotations having information content below p will be discarded")
    # mutex_parser_IC_THRESH.add_argument('--info_threshold_Phillip_Lord', '-PLTHRESH',
    #                                     help="Provide a threshold value t. All annotations having information content below t will be discarded")
    # parser.add_argument('--verbose', '-v',
    #                     help="Set this argument to 1 if you wish to view the outcome of each operation on console",
    #                     default=0)
    args = parser.parse_args()
    return args

def vprint(*s):
    global verbose
    # print(s,verbose)
    if verbose == 1:
        for string in s:
            print(string, end="")
        print()

def debias_inherit():
    # currently, take WyattClarkIC-perprotein.tsv as input (UniProt protein names and WC_IC)
    # TODO: handle gaf file as input (GOA._gaf20iterator)
    return

# def filter_most_coexp(filename, percentile):
#     all_interaction = pd.read_csv(filename, sep=",", header=0, usecols=['protein1', "protein2", "coexpression"])
#     all_coexp = all_interaction[all_interaction['coexpression'] > 0]
#     MAX_COEXP = max(all_coexp['coexpression'])
#     if percentile != None:
#         top_coexp = all_coexp[all_coexp['coexpression'] > percentile*(1-percentile/100)*MAX_COEXP]
#     return top_coexp

def filter_most_annotated(ranktable, idmatch, threshold):
    rankIC = pd.read_csv(ranktable, sep= "\t", header=0, names= ["alias", "WC_IC"])
    idtable = pd.read_csv(idmatch, sep = "\t", header=0, names=["string_id", "alias", "source"], usecols=["string_id", "alias"])
    rankIC_stringID = pd.merge(rankIC, idtable, on="alias", how="left")
    ranked = rankIC_stringID.drop_duplicates(keep="first").dropna()
    ranked_missing = rankIC_stringID[rankIC_stringID.isna().any(axis=1)]
    # TODO: resolve NaN (add UniProt ID table)
    # TODO: create dict metric["prot"] = ["string_id", "preferred_name", "count", "wc_ic", ...]
    most_annt_list = []
    for index, row in ranked.iterrows():
        if row['WC_IC'] >= threshold:
            most_annt_list.append(row['string_id'])
    return most_annt_list

def filter_least_annotated(ranktable, idmatch, threshold):
    rankIC = pd.read_csv(ranktable, sep="\t", header=0, names=["alias", "WC_IC"])
    idtable = pd.read_csv(idmatch, sep="\t", header=0, names=["string_id", "alias", "source"],
                          usecols=["string_id", "alias"])
    rankIC_stringID = pd.merge(rankIC, idtable, on="alias", how="left")
    ranked = rankIC_stringID.drop_duplicates(keep="first").dropna()
    ranked_missing = rankIC_stringID[rankIC_stringID.isna().any(axis=1)]
    # TODO: resolve NaN (add UniProt ID table)
    # TODO: create dict metric["prot"] = ["string_id", "preferred_name", "count", "wc_ic", ...]
    least_annt_list = []
    for index, row in ranked.iterrows():
        if row['WC_IC'] <= threshold:
            least_annt_list.append(row['string_id'])
    return least_annt_list

def associated_ignorome(all_ppi, tppi, most, least):
    all_ppi = pd.read_csv(all_ppi, header=0, delim_whitespace=True)
    coexp = all_ppi.loc[all_ppi.coexpression > tppi,['protein1','protein2','coexpression']]
    most_least = coexp[coexp['protein1'].isin(most) & coexp['protein2'].isin(least)]
    print(most_least.head(5))
    ignorome_list = most_least.iloc[:,1].tolist()
    print("Total ignorome found:", len(ignorome_list))
    return ignorome_list

    # TODO: separate module: idmapping

# def get_coexpression_api(coexpression_threshold, gene_all_common):
#     string_api_url = "https://version-11-5.string-db.org/api"
#     output_format = "tsv-no-header"
#     method = "network" # network or interaction_partners
#     request_url = "/".join([string_api_url, output_format, method])
#     my_genes = gene_all_common
#     params = {
#         "identifiers": "%0d".join(my_genes),  # list of proteins in common names
#         "species": 9606,  # species NCBI identifier
#         "limit" : 5,
#         "caller_identity": "ignoromenot"  # your app name
#     }
#     response = requests.post(request_url, data=params)
#     # print("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\tnscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore")
#     for line in response.text.strip().split("\n"):
#         l = line.strip().split("\t")
#         p1, p2= l[2], l[3]  # protein preferred name
#         # filter the interaction according to coexpression score (l[9])
#         coexpression_score = float(l[10])
#         if coexpression_score > coexpression_threshold:
#             ## print
#             print("\t".join([p1, p2, "coexpression (score. %.3f)" % coexpression_score]))


def main():
    # global verbose, options, report
    # commandLineArg = sys.argv
    # if len(commandLineArg) == 1:
    #     print("Please use the --help option to get usage information")
    # # Parse command line arguments
    # options = parseCommandLineArguments()
    # most = filter_most_annotated(options.ranktable, options.idtable, options.threshold_top)
    # least = filter_least_annotated(options.ranktable, options.idtable, options.threshold_bot)
    # print("List of ignorome genes based on threshold: ")
    # print(associated_ignorome(options.stringppi, options.threshold_ppi, most, least))

    ranktable = "WyattClarkIC-perprotein.tsv"
    idtable = "UniProtID_9606.protein.aliases.v11.5.txt"
    stringppi = "9606.protein.links.full.v11.5.txt"
    threshold_top = 1000
    threshold_bot = 5
    threshold_ppi = 500
    most = filter_most_annotated(ranktable, idtable,threshold_top)
    least = filter_least_annotated(ranktable, idtable, threshold_bot)
    print("List of ignorome genes based on threshold: ")
    print(associated_ignorome(stringppi, threshold_ppi, most, least))

main()