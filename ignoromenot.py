# This script finds ignorome that are highly-associated with highly-annotated genes
# Ranking based on annotation count, Philip Lord information content, PubMed mentions, Jensen score;
# Association based on STRING co-expression

from __future__ import print_function
from __future__ import division
import networkx as nx
from networkx.algorithms import bipartite
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


def parse_commandlinearguments():
    parser = argparse.ArgumentParser(prog="ignoromenot.py")
    mutex_parser_top_THRESH = parser.add_mutually_exclusive_group()
    mutex_parser_bot_THRESH = parser.add_mutually_exclusive_group()
    mutex_parser_ppi_THRESH = parser.add_mutually_exclusive_group()

    requiredArguments = parser.add_argument_group("Required arguments")
    requiredArguments.add_argument('--ranktable', '-rank', help="Ranking table with 1st column being gene names, 2nd column being Wyatt Clark IC",
                                   required=True)
    requiredArguments.add_argument('--idtable', '-id', help="ID mapping table from STRING",
                                   required=True)
    requiredArguments.add_argument('--stringppi', '-ppi', help="Protein-protein interaction network from STRING",
                                   required=True)
    mutex_parser_top_THRESH.add_argument('--threshold_top', '-ttop', help="The threshold level above which most-annotated genes will be selected",
                        default=100)
    mutex_parser_top_THRESH.add_argument('--percentile_top', '-ptop', help="Genes at k-th percentile and above will be selected. Example: -ptop 95 selects top 5% of genes with highest value")
    mutex_parser_bot_THRESH.add_argument('--threshold_bot', '-tbot', help="The threshold level under which least-annotated genes will be selected",
                        default=5)
    mutex_parser_bot_THRESH.add_argument('--percentile_bot', '-pbot', help="Genes at k-th percentile and below will be selected. Example: -pbot 10 selects top 10% of genes with lowest value")
    mutex_parser_ppi_THRESH.add_argument('--threshold_ppi', '-tppi', help="The threshold value (0-1000) above which STRING interaction scores will be selected",
                        default=500)
    mutex_parser_ppi_THRESH.add_argument('--percentile_ppi', '-pppi', help="STRING interaction pairs at k-th percentile and above will be selected. Example: -pppi 95 selects top 5% of associated pairs with highest score")
    args = parser.parse_args()
    return args


# def vprint(*s):
#     global verbose
#     # print(s,verbose)
#     if verbose == 1:
#         for string in s:
#             print(string, end="")
#         print()


def debias_inherit():
    # currently, take WyattClarkIC-perprotein.tsv as input (UniProt protein names and WC_IC)
    # TODO: handle gaf file as input (GOA._gaf20iterator)
    return


def filter_most_annotated(ranktable, ttop, ptop):
    rankIC = pd.read_csv(ranktable, sep= "\t", header=0, names= ["alias", "WC_IC"])
    # idtable = pd.read_csv(idmatch, sep = "\t", header=0, names=["string_id", "alias", "source"], usecols=["string_id", "alias"])
    # rankIC_stringID = pd.merge(rankIC, idtable, on="alias", how="left")
    # ranked = rankIC_stringID.drop_duplicates(keep="first").dropna()
    # ranked_missing = rankIC_stringID[rankIC_stringID.isna().any(axis=1)]
    # TODO: resolve NaN (add UniProt ID table)
    # TODO: create dict metric["prot"] = ["string_id", "preferred_name", "count", "wc_ic", ...]
    most_annt_list = []
    if ptop is not None:
        threshold = np.percentile(rankIC["WC_IC"], float(ptop))
    else:
        threshold = float(ttop)
    for index, row in rankIC.iterrows():
        if row['WC_IC'] >= threshold:
            most_annt_list.append(row['alias'])
    print("IC threshold for most annotated genes:", threshold)
    return most_annt_list


def filter_least_annotated(ranktable, tbot, pbot):
    rankIC = pd.read_csv(ranktable, sep="\t", header=0, names=["alias", "WC_IC"])

    # ranked_missing = rankIC_stringID[rankIC_stringID.isna().any(axis=1)]
    # TODO: resolve NaN (add UniProt ID table)
    # TODO: create dict metric["prot"] = ["string_id", "preferred_name", "count", "wc_ic", ...]
    least_annt_list = []
    if pbot is not None:
        threshold = np.percentile(rankIC["WC_IC"], float(pbot))
    else:
        threshold = float(tbot)
    for index, row in rankIC.iterrows():
        if row['WC_IC'] <= threshold:
            least_annt_list.append(row['alias'])
    print("IC threshold for least annotated genes:", threshold)
    return least_annt_list


def associated_ignorome(all_ppi, tppi, pppi, most, least, idtable, ranktable):
    all_ppi = pd.read_csv(all_ppi, header=0, delim_whitespace=True)
    coexp = all_ppi.loc[all_ppi.coexpression > 0, ['protein1', 'protein2', 'coexpression']]
    if pppi is not None:
        threshold = np.percentile(coexp['coexpression'], float(pppi))
    else:
        threshold = float(tppi)
    print("Coexpression threshold:", threshold)

    # Filter most coexpressed
    coexp = all_ppi.loc[all_ppi.coexpression > threshold, ['protein1', 'protein2', 'coexpression']]
    string_most, idmapping = common_to_string(idtable, ranktable, most)  # idmatch is path
    # print("Top IC genes: ", string_most)
    string_least, idmapping = common_to_string(idtable, ranktable, least)
    # print("Bot IC genes: ", string_least)
    most_least = coexp[coexp['protein1'].isin(string_most) & coexp['protein2'].isin(string_least)].sort_values(by= 'coexpression', ascending=False)
    known_set = set(most_least.iloc[:, 0].tolist())
    ignorome_set = set(most_least.iloc[:, 1].tolist())
    print("\nStrongest interaction pairs:")
    print(most_least.head(20))

    # Print results with STRING and preferred name
    ignorome_table = idmapping[idmapping['string_id'].isin(ignorome_set)]
    known_table = idmapping[idmapping['string_id'].isin(known_set)]
    print("\nTotal ignorome found:", len(ignorome_set), ignorome_table['alias'].to_list())
    print(ignorome_table)
    return ignorome_set, known_set


def bipartite_graph(most, least):
    G = nx.Graph()
    G.add_nodes_from(most, bipartite=0)
    G.add_nodes_from(least, bipartite=1)
    # TODO: get edges (list of tuple pairs) from all_ppi table


def common_to_string(idtable, ranktable, genelist):
    ranktable = pd.read_csv(ranktable, sep="\t", header=0, names=["alias", "WC_IC"])
    idtable = pd.read_csv(idtable, sep="\t", header=0, names=["string_id", "alias", "source"],
                          usecols=["string_id", "alias"])
    idmapping = pd.merge(ranktable, idtable, on="alias", how="left")
    idmapping = idmapping.drop_duplicates(keep="first").dropna()

    returnlist = []
    for gene in genelist:
        returnlist.append(idmapping[idmapping['alias'] == gene]['string_id'].to_string(index = False))
    return returnlist, idmapping


def get_coexpression_api(coexpression_threshold, gene_all_common):
    # TODO: replace hard code
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network" # network or interaction_partners
    request_url = "/".join([string_api_url, output_format, method])
    my_genes = gene_all_common
    params = {
        "identifiers": "%0d".join(my_genes),  # list of proteins in common names
        "species": 9606,  # species NCBI identifier
        "limit" : 5,
        "caller_identity": "ignoromenot"  # your app name
    }
    response = requests.post(request_url, data=params)
    print("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\t"
          "ncbiTaxonId\tscore\tnscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore")
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        p1, p2= l[2], l[3]  # protein preferred name
        # filter the interaction according to coexpression score (l[9])
        coexpression_score = float(l[10])
        if coexpression_score > coexpression_threshold:
            ## print
            print("\t".join([p1, p2, "coexpression (score. %.3f)" % coexpression_score]))


def main():
    global options
    commandLineArg = sys.argv
    if len(commandLineArg) == 1:
        print("Please use the --help option to get usage information")
    # Parse command line arguments
    options = parse_commandlinearguments()
    # # Execute ranking and find ignoromes
    most = filter_most_annotated(options.ranktable, options.threshold_top, options.percentile_top)
    least = filter_least_annotated(options.ranktable, options.threshold_bot, options.percentile_bot)
    associated_ignorome(options.stringppi, options.threshold_ppi, options.percentile_ppi, most, least, options.idtable, options.ranktable)

    # # Run from IDE
    # ranktable = "WyattClarkIC-perprotein.tsv"
    # idtable = "511145.protein.aliases.v11.5.txt"
    # stringppi = "511145.protein.links.full.v11.5.txt"
    # threshold_top = 100
    # threshold_bot = None
    # threshold_ppi = None
    # percentile_top = None
    # percentile_bot = 0.5
    # percentile_ppi = 95
    # most = filter_most_annotated(ranktable,threshold_top, percentile_top)
    # least = filter_least_annotated(ranktable, threshold_bot, percentile_bot)
    # associated_ignorome(stringppi, threshold_ppi, percentile_ppi, most, least, idtable, ranktable)

if __name__ == "__main__":
    main()