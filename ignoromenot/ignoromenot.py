#!/usr/bin/env python
# This script finds ignorome that are highly-associated with highly-annotated genes
# Ranking currently based on Wyatt Clark information content -> based on GAF-calculated files
# Association currently based on STRING co-expression

from __future__ import print_function
from __future__ import division
import collections
import datetime

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
from dateutil import parser
import debias_supp as debias
import pickle as cp

# Global variables
verbose = 0
options = ""
STRING_API_URL = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "interaction_partners"  # network or interaction_partners
request_url = "/".join([STRING_API_URL, output_format, method])

# Files needed for debias:
DATADIR = "data_for_ignoromenot/" # TODO: run debias_prep instead of providing readily used files
# os.chdir('bar/')
# Some filenames
FILE_ALTERNATE_ID_TO_ID_MAPPING = DATADIR + "alt_to_id.graph"
FILE_MFO_ONTOLOGY_GRAPH = DATADIR + "mf.graph"
FILE_BPO_ONTOLOGY_GRAPH = DATADIR + "bp.graph"
FILE_CCO_ONTOLOGY_GRAPH = DATADIR + "cc.graph"
FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH = DATADIR + "mf_ancestors.map"
FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH = DATADIR + "bp_ancestors.map"
FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH = DATADIR + "cc_ancestors.map"

debias.FILE_ALTERNATE_ID_TO_ID_MAPPING = FILE_ALTERNATE_ID_TO_ID_MAPPING
debias.FILE_MFO_ONTOLOGY_GRAPH = FILE_MFO_ONTOLOGY_GRAPH
debias.FILE_BPO_ONTOLOGY_GRAPH = FILE_BPO_ONTOLOGY_GRAPH
debias.FILE_CCO_ONTOLOGY_GRAPH = FILE_CCO_ONTOLOGY_GRAPH
debias.FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH = FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH
debias.FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH = FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH
debias.FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH = FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH


def parse_commandlinearguments():
    """Parsing command line arguments: required input files: rank tsv, STRING alias txt and STRING protein.links txt,
    and input threshold (otherwise use default).
    """
    parser = argparse.ArgumentParser(prog="ignoromenot.py", description="IgnoroMeNot find ignorome genes",add_help=True)
    mutex_parser_top_thresh = parser.add_mutually_exclusive_group()
    mutex_parser_bot_thresh = parser.add_mutually_exclusive_group()
    mutex_parser_ppi_thresh = parser.add_mutually_exclusive_group()
    specifier_parser = parser.add_argument_group("Specifiers")

    required_arguments = parser.add_argument_group("Required arguments")
    required_arguments.add_argument('--input', '-i', required=True,
                                    help="Input could be a path to a tab-separated file of a list of genes with any single "
                                         "annotation metric (see above sample rank table) or\n"
                                         "Input could be a GO Annotation File (GAF). In this case, filename must end with .gaf")
    required_arguments.add_argument('--idtable', '-id', required=True,
                                    help="The path to a STRING protein alias file of the organism being examined.\n"
                                         "The filename must start with the organism ID (e.g., 9606 for human, 511145 for E.coli)")
    required_arguments.add_argument('--stringppi', '-ppi', required=True,
                                    help="The path to a STRING interaction network of the organism being examined.\n"
                                         "The filename must start with the organism ID (e.g., 9606 for human, 511145 for E.coli)")
    specifier_parser.add_argument('--metric', '-m', default="ct",
                                    help="A single annotation metric obtained from GAF based on which the proteins will be ranked. If not specified,"
                                         "proteins will be ranked based on annotation count.\n"
                                         "Accepted metrics: ct, ic, ia. Default to 'ct' ")
    specifier_parser.add_argument('--aspect', '-a', default="All",
                                    help="If GAF file is provided as input, specify which aspect of GO to rank the proteins. If not specified, "
                                         "proteins will be ranked based on the total value across 3 aspects.\n"
                                         "Accepted aspects: All, MFO, BPO, CCO. Default to 'All' (sum accross 3 aspects).")
    mutex_parser_top_thresh.add_argument('--threshold_top', '-ttop', default=100,
                                         help="Set an absolute upper threshold for most annotated genes based on the given metric.\nDefault to 100")
    mutex_parser_top_thresh.add_argument('--percentile_top', '-ptop',
                                         help="Set a relative upper threshold for most annotated genes at k-th percentile based on the given metric.\n"
                                              "Cannot be provided simultaneously with --threshold_top."
                                              "Example: -ptop 95 selects top 5%% of genes with highest value")
    mutex_parser_bot_thresh.add_argument('--threshold_bot', '-tbot', default=5,
                                         help="Set an absolute lower threshold for least annotated genes based on the given metric.\nDefault to 5")
    mutex_parser_bot_thresh.add_argument('--percentile_bot', '-pbot',
                                         help="Set a relative lower threshold for least annotated genes at k-th percentile based on the given metric.\n"
                                              "Cannot be provided simultaneously with --threshold_bot."
                                              "Example: -pbot 10 selects top 10%% of genes with lowest value")
    mutex_parser_ppi_thresh.add_argument('--threshold_ppi', '-tppi', default=500,
                                         help="Set an absolute upper threshold for STRING protein-protein interaction score.\nDefault to 500")
    mutex_parser_ppi_thresh.add_argument('--percentile_ppi', '-pppi',
                                         help="Set a relative upper threshold for STRING protein-protein interaction score.\n"
                                              "Cannot be provided simultaneously with --threshold_ppi."
                                              "Example: -pppi 95 selects top 5%% of associated pairs with highest score")
    args = parser.parse_args()
    return args


# def vprint(*s):
#     global verbose
#     # print(s,verbose)
#     if verbose == 1:
#         for string in s:
#             print(string, end="")
#         print()

# __init__.py (hash and hash)
# intermediate file keeps all the metrics
# make a new branch on gothresher, modify it
# use different name, gocrusher

def rank_gaf_file(gaf_input, metric, aspect):
    # TODO: handle multiple input files (enumerate(options.input))
    """Takes 1 GAF file as input and return ranked data frame of 4 per-protein metrics:
    annotation count, article count, information content, propagated information content
    """
    # gaf_input = '/home/ahphan/RotationData/Friedberg/gaf_files/Only_Experimental_Annotations/goa_human27Jul22Filtered.gaf'
    # metric = "ct"
    # aspect = "MFO"
    gaf = GOA._gaf20iterator(open(gaf_input, "r"))
    data = debias.convertFromGAFToRequiredFormat(gaf)
    # Prot_to_GO_Map, all_GO_Terms_in_corpus = debias.createProteinToGOMapping(data)
    # Prot_to_GO_Map_propagated = debias.propagateOntologies(Prot_to_GO_Map)
    prev_go_term_freq = debias.freqGO_TERM(data)
    data, goia_f, goia_p, goia_c = debias.calculateWyattClarkInformationContent(data, 1, 0, None, './', 0)
    new_data, discarded_data = debias.chooseProteinsBasedOnPublications(data, 0, None)
    later_go_term_freq = debias.freqGO_TERM(data)
    ct_dict, ic_dict, ia_dict = debias.generateHistogram(data, data, discarded_data, 0, prev_go_term_freq,
                      later_go_term_freq, goia_f, goia_p, goia_c)
    data_dict = str(metric) + "_dict"
    ranktable = pd.DataFrame.from_dict(ct_dict,orient='index')
    ranktable.columns = [str(metric)+"_All", str(metric)+"_F", str(metric)+"_P", str(metric)+"_C"]
    ranktable.reset_index(inplace=True)
    ranktable=ranktable.rename(columns={"index":"Genes"})
    if aspect == "All":
        ranktable=ranktable.iloc[0,[0,1]]
    elif aspect == "MFO":
        ranktable=ranktable.iloc[:,[0,2]]
    elif aspect == "BPO":
        ranktable=ranktable.iloc[:,[0,3]]
    elif aspect == "CCO":
        ranktable=ranktable.iloc[:,[0,4]]
    else:
        print("Error in aspect code. Please choose one of the following: All, MFO, BPO, CCO")
        exit()
    return ranktable


def filter_most_annotated(ranktable, ttop, ptop):
    """
    This function filters most annotated genes based on threshold
    :param ranktable: tab-separated gene list with one or more metrics
    :param ttop: absolute threshold
    :param ptop: percentile threshold
    :return: list of top genes
    """
    # allrank = pd.read_csv(ranktable, sep="\t", header=0) # TODO: parse ranktable here, categorize from main with options.input
    allrank = ranktable
    most_annt_list = []

    # Get threshold and filter for each metric
    if ptop is not None:
        threshold = np.percentile(allrank.iloc[:, 1], float(ptop))
    else:
        threshold = float(ttop)
    for index, row in allrank.iterrows():  # another loop if >1 metric
        if row[1] >= threshold:
            most_annt_list.append(row[0])
    print()
    print("Absolute threshold for most annotated genes:", threshold)
    return most_annt_list


def filter_least_annotated(ranktable, tbot, pbot):
    """
    This function filters least annotated genes based on threshold
    :param ranktable: tab-separated gene list with one or more metrics
    :param tbot: absolute threshold
    :param pbot: percentile threshold
    :return: list of bot genes
    """
    # allrank = pd.read_csv(ranktable, sep="\t", header=0)
    allrank = ranktable
    least_annt_list = []

    # Get threshold and filter for each metric
    if pbot is not None:
        threshold = np.percentile(allrank.iloc[:, 1], float(pbot))
    else:
        threshold = float(tbot)
    for index, row in allrank.iterrows():  # another loop if >1 metric
        if row[1] <= threshold:
            least_annt_list.append(row[0])
    print("Absolute threshold for least annotated genes:", threshold)
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
    :param ranktable: pd.DataFrame ranktable with 2 columns: Genes and metric
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
    ignorome_common = []
    for ignorome in ignorome_set:
        common = idmapping.loc[idmapping['string_id']==ignorome, 'alias']
        ignorome_common.append(common)
    # print("\nStrongest most-least annotated pairs:\n", most_least.head(10))

    # Print results with STRING ids, information value
    ignorome_table = idmapping[idmapping['string_id'].isin(ignorome_set)]
    known_table = idmapping[idmapping['string_id'].isin(known_set)]
    print("\nTotal ignorome found:", len(ignorome_common), ignorome_common)
    print(ignorome_table)
    ignorome_table.to_csv('./Ignorome_Genes_' + str(datetime.datetime.now) + '.tsv', sep="\t", header=True)
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
    :param ranktable: pd.DataFrame ranktable with 2 columns: Genes (common names) and metric
    :param genelist: input common gene names
    :return: gene STRING ids and idmapping table
    """
    ranktable.columns.values[0] = "alias"
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
    This function uses STRING API to get interaction partners of target genes (with internet)
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
    options = parse_commandlinearguments()
    if len(sys.argv) < 4:
        print("Please use the --help option to get usage information")
    # Parse command line arguments


    # # Execute ranking from GAF file
    # ranktable = rank_gaf_file(options.input)
    if options.input.split("/")[-1].split(".")[-1] != "gaf":
        # Find ignorome from ranktable
        ranktable = pd.read_csv(options.input, sep="\t")
    else:
        ranktable = rank_gaf_file(options.input, options.metric, options.aspect)

    most = filter_most_annotated(ranktable, options.threshold_top, options.percentile_top)
    least = filter_least_annotated(ranktable, options.threshold_bot, options.percentile_bot)
    coexp_threshold, ignorome_list, known_list = associated_ignorome(options.stringppi, options.threshold_ppi,
                                                                     options.percentile_ppi, most, least,
                                                                     options.idtable, ranktable)

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
