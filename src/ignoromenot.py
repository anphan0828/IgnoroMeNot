# This script finds ignorome that are highly-associated with highly-annotated genes
# Ranking currently based on Wyatt Clark information content
# Association currently based on STRING co-expression

from __future__ import print_function
from __future__ import division
import networkx as nx
from networkx.algorithms import bipartite
import sys
import matplotlib.pyplot as plt
import argparse
import math
import numpy as np
from numpy import percentile
import Bio
from Bio.UniProt import GOA
from Bio.Seq import Seq
import os
import requests
import pandas as pd

# Global variables
verbose = 0
options = ""
STRING_API_URL = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "interaction_partners"  # network or interaction_partners
request_url = "/".join([STRING_API_URL, output_format, method])


def parse_commandlinearguments():
    """Parsing command line arguments: required input files: rank tsv, STRING alias txt and STRING protein.links txt,
    and input threshold (otherwise use default).
    """
    parser = argparse.ArgumentParser(prog="ignoromenot.py")
    mutex_parser_top_thresh = parser.add_mutually_exclusive_group()
    mutex_parser_bot_thresh = parser.add_mutually_exclusive_group()
    mutex_parser_ppi_thresh = parser.add_mutually_exclusive_group()

    required_arguments = parser.add_argument_group("Required arguments")
    required_arguments.add_argument('--ranktable', '-rank', required=True,
                                    help="Ranking table with 1st column being gene names, other columns being metrics")
    required_arguments.add_argument('--idtable', '-id', required=True,
                                    help="ID mapping table from STRING")
    required_arguments.add_argument('--stringppi', '-ppi', help="Protein-protein interaction network from STRING",
                                    required=True)
    mutex_parser_top_thresh.add_argument('--threshold_top', '-ttop', default=100,
                                         help="The threshold level above which most-annotated genes will be selected")
    mutex_parser_top_thresh.add_argument('--percentile_top', '-ptop',
                                         help="Genes at k-th percentile and above will be selected. "
                                              "Example: -ptop 95 selects top 5% of genes with highest value")
    mutex_parser_bot_thresh.add_argument('--threshold_bot', '-tbot', default=5,
                                         help="The threshold level under which least-annotated genes will be selected")
    mutex_parser_bot_thresh.add_argument('--percentile_bot', '-pbot',
                                         help="Genes at k-th percentile and below will be selected. "
                                              "Example: -pbot 10 selects top 10% of genes with lowest value")
    mutex_parser_ppi_thresh.add_argument('--threshold_ppi', '-tppi', default=500,
                                         help="The threshold value (0-1000) above which "
                                              "STRING interaction scores will be selected")
    mutex_parser_ppi_thresh.add_argument('--percentile_ppi', '-pppi',
                                         help="STRING interaction pairs at k-th percentile and above will be selected. "
                                              "Example: -pppi 95 selects top 5% of associated pairs with highest score")
    args = parser.parse_args()
    return args


# def vprint(*s):
#     global verbose
#     # print(s,verbose)
#     if verbose == 1:
#         for string in s:
#             print(string, end="")
#         print()


def rank_gaf_file(gaf_input):
    # TODO: handle multiple input files (enumerate(options.input))
    """Takes 1 GAF file as input and return ranked data frame of 4 per-protein metrics:
    annotation count, article count, information content, propagated information content
    """
    gaf_output = GOA._gaf20iterator(open(gaf_input, "r"))  # generator object, each row is a dict
    # TODO: aggregate metrics per protein and get ranking
    return gaf_output


def filter_most_annotated(ranktable, ttop, ptop):
    """
    This function filters most annotated genes based on threshold
    :param ranktable: tab-separated gene list with one or more metrics
    :param ttop: absolute threshold
    :param ptop: percentile threshold
    :return: list of top genes
    """
    allrank = pd.read_csv(ranktable, sep="\t", header=0)
    # TODO: create dict metric["prot"] = ["string_id", "preferred_name", "count", "wc_ic", ...]
    most_annt_list = []

    # Get threshold and filter for each metric
    if ptop is not None:
        threshold = np.percentile(allrank.iloc[:, 1], float(ptop))
    else:
        threshold = float(ttop)
    for index, row in allrank.iterrows():  # another loop if >1 metric
        if row[1] >= threshold:
            most_annt_list.append(row[0])
    print("IC threshold for most annotated genes:", threshold)
    return most_annt_list


def filter_least_annotated(ranktable, tbot, pbot):
    """
    This function filters least annotated genes based on threshold
    :param ranktable: tab-separated gene list with one or more metrics
    :param tbot: absolute threshold
    :param pbot: percentile threshold
    :return: list of bot genes
    """
    allrank = pd.read_csv(ranktable, sep="\t", header=0)
    least_annt_list = []

    # Get threshold and filter for each metric
    if pbot is not None:
        threshold = np.percentile(allrank.iloc[:, 1], float(pbot))
    else:
        threshold = float(tbot)
    for index, row in allrank.iterrows():  # another loop if >1 metric
        if row[1] <= threshold:
            least_annt_list.append(row[0])
    print("IC threshold for least annotated genes:", threshold)
    return least_annt_list


def associated_ignorome(all_ppi, tppi, pppi, most, least, idtable, ranktable):
    """
    This function finds ignorome genes based on association with high information genes
    :param all_ppi: path to STRING interaction file
    :param tppi: absolute threshold
    :param pppi: percentile threshold
    :param most: list of top genes
    :param least: list of bot genes
    :param idtable: path to STRING alias file (or another id mapping table)
    :param ranktable: path to ranktable
    :return: interaction threshold, ignorome genes, highly-annotated genes interacting with ignorome
    """
    all_ppi = pd.read_csv(all_ppi, header=0, delim_whitespace=True)
    coexp = all_ppi.loc[all_ppi.coexpression > 0, ['protein1', 'protein2', 'coexpression']]  # or other interaction
    if pppi is not None:
        threshold = np.percentile(coexp['coexpression'], float(pppi))
    else:
        threshold = float(tppi)
    print("Coexpression threshold:", threshold)

    # Filter most coexpressed pairs
    coexp = all_ppi.loc[all_ppi.coexpression > threshold, ['protein1', 'protein2', 'coexpression']]
    string_most, idmapping = common_to_string(idtable, ranktable, most)
    # print("Top IC genes: ", string_most)
    string_least, idmapping = common_to_string(idtable, ranktable, least)
    # print("Bot IC genes: ", string_least)

    # Get ignorome genes: have low knowledge/interest but interact strongly with high knowledge/interest genes
    most_least = coexp[coexp['protein1'].isin(string_most) & coexp['protein2'].isin(string_least)].\
        sort_values(by='coexpression', ascending=False)
    known_set = set(most_least.iloc[:, 0].tolist())
    ignorome_set = set(most_least.iloc[:, 1].tolist())
    # print("\nStrongest most-least annotated pairs:\n", most_least.head(10))

    # Print results with STRING ids, information value
    ignorome_table = idmapping[idmapping['string_id'].isin(ignorome_set)]
    known_table = idmapping[idmapping['string_id'].isin(known_set)]
    print("\nTotal ignorome found:", len(ignorome_set), ignorome_table['string_id'].to_list())
    print(ignorome_table)
    return threshold, ignorome_set, known_set


def bipartite_graph(most, least):
    """
    :param most: list of top genes
    :param least: list of bot genes
    :return:
    """
    G = nx.Graph()
    G.add_nodes_from(most, bipartite=0)
    G.add_nodes_from(least, bipartite=1)
    # TODO: get edges (list of tuple pairs) from all_ppi table


def common_to_string(idtable, ranktable, genelist):
    """
    This function takes UniProt protein symbols and returns STRING ids
    :param idtable: path to STRING alias file
    :param ranktable: path to ranktable (common gene names)
    :param genelist: input common gene names
    :return: gene STRING ids and idmapping table
    """
    ranktable = pd.read_csv(ranktable, sep="\t", header=0).rename(columns={"DB_Object_Symbol": "alias"})
    idtable = pd.read_csv(idtable, sep="\t", header=0, names=["string_id", "alias", "source"],
                          usecols=["string_id", "alias"])
    all_idmapping = pd.merge(ranktable, idtable, on="alias", how="left")
    idmapping = all_idmapping.drop_duplicates(keep="first").dropna()
    # na_idmapping = all_idmapping[all_idmapping.isna().any(axis=1)]
    # TODO: resolve NaN (input another mapping table from UniProt)

    returnlist = []
    for gene in genelist:
        returnlist.append(idmapping[idmapping['alias'] == gene]['string_id'].to_string(index=False))
    return returnlist, idmapping


def get_network_api(calculated_threshold, gene_all_common, species):
    """
    This function uses STRING API to get interaction partners of target genes
    :param calculated_threshold: interaction threshold
    :param gene_all_common: list of ignorome genes
    :param species: species ID
    :return
    """
    my_genes = gene_all_common
    params = {
        "identifiers": "%0d".join(my_genes),  # list of proteins in common names
        "species": species,  # species NCBI identifier, obtained from alias file
        "limit": 10,
        "caller_identity": "ignoromenot"  # your app name
    }
    response = requests.post(request_url, data=params)
    print("\nInteraction network of ignorome genes:")
    print("QueryIgnorome\tPartners\tInteraction")
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        p1, p2 = l[0], l[1]  # protein preferred name
        # filter the interaction according to coexpression score (l[9])
        coexpression_score = float(l[10])
        if coexpression_score > float(calculated_threshold)/1000:
            # print
            print("\t".join([p1, p2, "coexpression (score. %.3f)" % coexpression_score]))
    return


def main():
    global options
    command_line_arg = sys.argv
    if len(command_line_arg) == 1:
        print("Please use the --help option to get usage information")
    # Parse command line arguments
    options = parse_commandlinearguments()
    # # Execute ranking from GAF file
    # ranktable = rank_gaf_file(options.input)

    # Find ignorome from ranktable
    most = filter_most_annotated(options.ranktable, options.threshold_top, options.percentile_top)
    least = filter_least_annotated(options.ranktable, options.threshold_bot, options.percentile_bot)
    coexp_threshold, ignorome_list, known_list = associated_ignorome(options.stringppi, options.threshold_ppi,
                                                                     options.percentile_ppi, most, least,
                                                                     options.idtable, options.ranktable)
    species = options.idtable.split("/")[-1].split(".")[0]
    get_network_api(coexp_threshold, ignorome_list, species)

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
