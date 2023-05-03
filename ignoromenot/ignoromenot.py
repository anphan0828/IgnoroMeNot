#!/usr/bin/env python
# This script finds ignorome that are strongly associated with highly-annotated genes
# Ranking currently based on Wyatt Clark information content
# Association currently based on STRING co-expression, database, experiments

from __future__ import print_function
from __future__ import division
import collections
import datetime
import networkx as nx
from networkx.algorithms import bipartite
import sys
import matplotlib.pyplot as plt
import seaborn as sns
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
import uniprot_api as uniprot
import pickle as cp
import time
import warnings
with warnings.catch_warnings():
    warnings.simplefilter(action='ignore', category=FutureWarning)


# Global variables
verbose = 0
options = ""
STRING_API_URL = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "interaction_partners"  # network or interaction_partners
request_url = "/".join([STRING_API_URL, output_format, method])

# Files needed for debias:
DATADIR = "data_for_ignoromenot/"
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
    required_arguments.add_argument('--stringppi', '-ppi', required=True,
                                    help="The path to a STRING interaction network of the organism being examined.\n"
                                         "The filename must start with the organism ID (e.g., 9606 for human, 511145 for E.coli)")
    specifier_parser.add_argument('--metric', '-m', default="ia",
                                    help="A single annotation metric obtained from GAF based on which the proteins will be ranked. If not specified,"
                                         "proteins will be ranked based on annotation count.\n"
                                         "Accepted metrics: ct, ia, fc. Default to 'ia' ")
    specifier_parser.add_argument('--aspect', '-a', default="MFO",
                                    help="If GAF file is provided as input, specify which aspect of GO to rank the proteins. If not specified, "
                                         "proteins will be ranked based on the total value across 3 aspects.\n"
                                         "Accepted aspects: All, MFO, BPO, CCO. Default to 'MFO' (molecular function).")
    specifier_parser.add_argument('--recal', '-r', default=1,
                                    help="Set recal to 1 if you are running on a new GAF, the program reconstructs the ancestors and recalculates all metrics."
                                         "Set recal to 0 if you only change the threshold values but run on the same GAF, the program uses old files to save time."
                                         "Default to 1.")
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


def rank_gaf_file(gaf_input, metric, aspect, recal, swiss_no_gaf_path=None):
    # TODO: handle multiple input files (enumerate(options.input))
    """Takes 1 GAF file as input and return ranked data frame of 4 per-protein metrics:
    annotation count, article count, information content, propagated information content
    """
    gaf = GOA._gaf20iterator(open(gaf_input, "r"))
    data = debias.convertFromGAFToRequiredFormat(gaf)
    data, Prot_to_GO_Map_propagated, ontology_to_ia_map, Prot_IA_propagated, Prot_Count = \
        debias.calculateWyattClarkInformationContent(data, int(recal), 0, None, './', 0, "All")
    new_data, discarded_data = debias.chooseProteinsBasedOnPublications(data, 100, None)
    debias.calculateWyattClarkInformationContent(new_data, int(recal), 0, None, './', 0, "LTP")
    debias.calculateWyattClarkInformationContent(discarded_data, int(recal), 0, None, './', 0, "HTP")
    fc_dict = debias.calculateFractionalCountPerProtein(data, type="All")  # before filtered HTP, LTP
    if metric == "ct":
        ranktable = pd.DataFrame.from_dict(Prot_Count,orient='index')
    elif metric == "fc":
        ranktable = pd.DataFrame.from_dict(fc_dict,orient='index')
    elif metric == "ia":
        ranktable = pd.DataFrame.from_dict(Prot_IA_propagated,orient='index')
    else:
        print("Please provide a valid metric that can be calculated from GAF: ct, fc, ia")
        exit()
    ranktable.columns = [str(metric)+"_All", str(metric)+"_F", str(metric)+"_P", str(metric)+"_C"]
    ranktable.reset_index(inplace=True)
    ranktable=ranktable.rename(columns={"index":"protein_name"})
    go_accreted_dict = dict()
    for GO_term, value in ontology_to_ia_map.items():
        try:
            go_accreted_dict[GO_term] = value[2]
        except IndexError:
            continue
    specific_term = 6
    ia_blacklist = []
    for prot, GO_terms in Prot_to_GO_Map_propagated.items():
        for term in GO_terms:
            if go_accreted_dict[term] > specific_term:
                ia_blacklist.append(prot.split("_")[-1])
                break
    # print(ia_blacklist)
    if aspect == "All":
        ranktable_asp=ranktable.iloc[:,[0,1]]
    elif aspect == "MFO":
        ranktable_asp=ranktable.iloc[:,[0,2]]
    elif aspect == "BPO":
        ranktable_asp=ranktable.iloc[:,[0,3]]
    elif aspect == "CCO":
        ranktable_asp=ranktable.iloc[:,[0,4]]
    else:
        print("Error in aspect code. Please choose one of the following: All, MFO, BPO, CCO")
        exit()

    # Add SwissProt proteins that are not in GAF:
    if swiss_no_gaf_path is not None:
        no_gaf = open(swiss_no_gaf_path, 'r')
        allLines = no_gaf.readlines()
        for protein in allLines:
            name = protein.strip()
            columnName = ranktable_asp.columns.values[1]
            df2=pd.DataFrame([[name,0]],columns=['protein_name',columnName])
            ranktable_asp=pd.concat([ranktable_asp,df2],ignore_index=True)
    return ranktable_asp, ia_blacklist


def filter_most_annotated(ranktable, ttop, ptop):
    """
    This function filters most annotated genes based on threshold
    :param ranktable: tab-separated gene list with one or more metrics
    :param ttop: absolute threshold
    :param ptop: percentile threshold
    :return: list of top genes
    """
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
    print("Upper absolute threshold for " + str(allrank.columns.values[1]) +": "+ str(threshold)+". Number of genes: " + str(len(most_annt_list)))
    fh = open('HighlyAnnt' + str(threshold) + '.txt', 'w')
    for gene in most_annt_list:
        fh.write(gene + '\n')
    fh.close()
    return most_annt_list


def filter_least_annotated(ranktable, ia_blacklist, tbot, pbot):
    """
    This function filters least annotated genes based on threshold
    :param ranktable: tab-separated gene list with one or more metrics
    :param tbot: absolute threshold
    :param pbot: percentile threshold
    :return: list of bot genes
    """
    allrank = ranktable
    least_annt_list = []
    not_least_annt_list = []
    if pbot is not None:
        threshold = np.percentile(allrank.iloc[:, 1], float(pbot))
    else:
        threshold = float(tbot)
    for index, row in allrank.iterrows():  # another loop if >1 metric
        if row[1] <= threshold:
            if row[0] not in ia_blacklist:
                least_annt_list.append(row[0])
            else:
                not_least_annt_list.append(row[0])
    print("Lower absolute threshold for " + str(allrank.columns.values[1]) +": " + str(threshold)+". Number of genes: "+str(len(least_annt_list)))
    fh = open('LessAnnt'+str(threshold)+'.txt','w')
    for gene in least_annt_list:
        fh.write(gene+'\n')
    fh.close()
    print("Genes that have at least 1 specific experimental annotation not considered for less annotated: ", len(not_least_annt_list))
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
    string_least, idmapping = common_to_string(idtable, ranktable, least)

    most_least = coexp[coexp['protein1'].isin(string_most) & coexp['protein2'].isin(string_least)].\
        sort_values(by='protein2', ascending=True)
    known_set = set(most_least.iloc[:, 0].tolist())
    ignorome_set = set(most_least.iloc[:, 1].tolist())
    ignorome_common = []
    for ignorome in ignorome_set:
        common = idmapping.loc[idmapping['string_id']==ignorome, 'protein_name']
        ignorome_common.append((common.tolist()))

    # Print results with STRING ids, information value
    ignorome_table = idmapping[idmapping['string_id'].isin(ignorome_set)].sort_values(by='string_id',ascending=True)
    known_table = idmapping[idmapping['string_id'].isin(known_set)]
    print("\nIgnorome proteins found:", len(ignorome_common))
    print((ignorome_common))
    print(ignorome_table)
    ignorome_table.to_csv('./IgnoromeGenes_' + time.strftime(r"%Y-%m-%d-%H-%M-%S",time.localtime()) + '.tsv', sep="\t", header=True)
    return threshold, ignorome_set, known_set


def associated_ignorome2(all_ppi, tppi, pppi, most, least, ranktable):
    all_ppi = pd.read_csv(all_ppi, header=0, delim_whitespace=True)
    intppi = all_ppi.loc[((all_ppi.coexpression > 0) | (all_ppi.experiments > 0)) | (all_ppi.database > 0), ['protein1','protein2','coexpression','experiments','database']]
    if pppi is not None:
        threshold_coexp = np.percentile(intppi['coexpression'], float(pppi))
        threshold_exprm = np.percentile(intppi['experiments'], float(pppi))
        threshold_db = np.percentile(intppi['database'], float(pppi))
    else:
        threshold_coexp = float(tppi)  # TODO: remove absolute threshold or triplicate it
        threshold_exprm = 700
        threshold_db = 800
    print("Coexpression threshold:", threshold_coexp)
    print("Coexpression mean and median and std:", np.mean(intppi['coexpression']), np.median(intppi['coexpression']), np.std(intppi['coexpression']))
    print("Experiments threshold:", threshold_exprm)
    print("Experiments mean and median and std:", np.mean(intppi['experiments']), np.median(intppi['experiments']), np.std(intppi['experiments']))
    print("Database threshold:", threshold_db)
    print("Database mean and median and std:", np.mean(intppi['database']), np.median(intppi['database']), np.std(intppi['database']))

    # Filter most coexpressed pairs, or highest exprm score, or highest db score
    intppi_top = intppi.loc[((intppi.coexpression >= threshold_coexp) | (intppi.experiments >= threshold_exprm)) | (intppi.database >= threshold_db), ['protein1','protein2','coexpression','experiments','database']]
    print("Total interaction pairs in STRING: ", (intppi_top.shape[0])/2)
    intppi_unre = intppi.loc[((intppi.coexpression > 200) | (intppi.experiments > 0)) | (intppi.database > 0), ['protein1','protein2','coexpression','experiments','database']]

    try:
        mapping_dict = cp.load(open("data/temp/mapping_dict.txt", "rb"))
        mapping_dict_unre = cp.load(open("data/temp/mapping_dict_unre.txt", "rb"))
    except IOError as e:
        print("\nMapping ID with UniProt:")
        mapping_dict, mapping_dict_unre = uniprot.api_map_STRING_UniProt(all_ppi) # STRING to UniProtKB
        if os.path.isdir("data/temp/") == False:
            os.makedirs("data/temp/")
        cp.dump(mapping_dict, open("data/temp/mapping_dict.txt", "wb"))
        cp.dump(mapping_dict_unre, open("data/temp/mapping_dict_unre.txt", "wb"))

    string_most , id_rank, na_idrank = common_to_string2(mapping_dict, ranktable, most)
    string_least, id_rank, na_idrank = common_to_string2(mapping_dict, ranktable, least)
    print("\nProteins not mapped to STRING (no interaction data): " + str(len(na_idrank)))  # STRING uses synonym ygcB for cas3 gene
    na_idrank.to_csv('./Proteins_no_STRING.tsv',sep='\t',header=True,index=False)

    # Get ignorome genes: have low knowledge/interest but interact strongly with high knowledge/interest genes
    most_least = intppi_top[intppi_top['protein1'].isin(string_most) & intppi_top['protein2'].isin(string_least)]. \
        sort_values(by='protein2', ascending=True)
    known_set = set(most_least.iloc[:, 0].tolist())
    ignorome_set = set(most_least.iloc[:, 1].tolist())
    ignorome_common = []
    top_common = dict()
    for ignorome in ignorome_set:
        top = []
        for known in most_least[most_least['protein2'] == ignorome].loc[:,'protein1'].tolist():  # only print strongest association and over threshold top genes
            common = id_rank[id_rank['string_id'] == known].loc[:,'protein_name'].tolist()  # list of rich genes associated with ignorome genes, regardless of type of interaction
            top.extend(common) # unpack nested list into overall list
        top_common[ignorome] = tuple(top)
        common = id_rank.loc[id_rank['string_id'] == ignorome, 'protein_name']
        ignorome_common.append((common.tolist()))

    # Print results with STRING ids, information value
    ignorome_table = id_rank[id_rank['string_id'].isin(ignorome_set)].sort_values(by='string_id', ascending=True)
    known_table = id_rank[id_rank['string_id'].isin(known_set)]
    ignorome_table['associated_with'] = ignorome_table['string_id'].map(top_common)
    # TODO: how to distinguish
    print("Ignorome proteins found: " + str(len(ignorome_common)) + "\n")
    s = ignorome_table.associated_with.str.len().sort_values(ascending=False).index
    ignorome_table = ignorome_table.reindex(s)

    # Use mapping_dict_unre to find strongest coexp links between unreviewed and most annotated
    unre_set = set(mapping_dict_unre.keys())
    unre_common = []
    unre_coexp_top_common = dict()
    for unre in unre_set:
        coexp_top = []
        for known in intppi_unre[intppi_unre['protein2'] == unre].loc[:,'protein1'].tolist():  # only print strongest association and over threshold top genes
            if known in (known_set):
                common = id_rank[id_rank['string_id'] == known].loc[:,'protein_name'].tolist()
                coexp_top.extend(common)
        unre_coexp_top_common[unre] = tuple(coexp_top)
        common = mapping_dict_unre[unre]
        unre_common.append(common)

    # Print results with STRING ids, information value
    unre_table = intppi_unre[intppi_unre['protein1'].isin(known_set) & intppi_unre['protein2'].isin(unre_set)].sort_values(by='protein2', ascending=True)
    unre_table['associated_with'] = unre_table['protein2'].map(unre_coexp_top_common)
    unre_table['TrEMBL_ID'] = unre_table['protein2'].map(mapping_dict_unre)
    unre_table = unre_table.loc[:,['protein2', 'TrEMBL_ID', 'associated_with']].drop_duplicates(keep='first')
    s = unre_table.associated_with.str.len().sort_values(ascending=False).index
    unre_table = unre_table.reindex(s)
    return threshold_coexp, ignorome_set, known_set, ignorome_table, unre_table

def bipartite_graph(most, least, ignorome_table):
    """
    :param most: list of top genes
    :param least: list of bot genes
    :return:
    """
    G = nx.Graph()
    # G.add_nodes_from(most, bipartite=0)
    # G.add_nodes_from(least, bipartite=1)
    for index,row in ignorome_table.iterrows():
        # print(row)
        G.add_node(row['protein_name'],bipartite=1)
        for rich in row['associated_with']:
            G.add_node(rich,bipartite=0)
            G.add_edge(row['protein_name'],rich)

    edges = G.edges()
    m1=30
    n1=15
    pos = dict()
    pos.update((n, (i - m1/2, 1)) for i, n in enumerate(most))  # put nodes from X at x=1
    pos.update((n, (i-n1/2, 0)) for i, n in enumerate(least))  # put nodes from Y at x=2
    nx.draw_networkx(G, pos=pos, with_labels=True,font_size=6,node_size=100)
    plt.show()
    graph_df = nx.to_pandas_edgelist(G,source='Source',target='Target')
    graph_df.to_csv('Bipartite_Edgelist.tsv',sep="\t", header=True, index=False)
    print(graph_df.head())
    # plt.clf()
    # nx.draw_networkx(G, with_labels=True,font_size=6,node_size=100)
    # plt.show()

    # pos1=nx.bipartite_layout(G,nx.bipartite.sets(G)[0],align="horizontal",scale=1,aspect_ratio=5/2)
    # nx.draw(G, pos=pos1,with_labels=True)
    # plt.show()


def common_to_string2(iddict, ranktable, genelist):
    """
    This function takes UniProt protein symbols and returns STRING ids. This will be deprecated, replaced
    by UniProt ID mapping API
    :param idtable: path to STRING alias file
    :param ranktable: pd.DataFrame ranktable with 2 columns: Genes (common names) and metric
    :param genelist: input common gene names
    :return: gene STRING ids and idmapping table
    """
    idtable = pd.DataFrame.from_dict(iddict, orient='index')
    idtable.reset_index(inplace=True)
    idtable = idtable.rename(columns={"index": "string_id", 0:"protein_name"})
    id_rank_all = pd.merge(ranktable, idtable, on=['protein_name'], how="left")
    id_rank = id_rank_all.dropna()
    na_idrank = id_rank_all[id_rank_all.isna().any(axis=1)].sort_values(by=id_rank_all.columns[1], ascending=False)
    returnlist = []
    for gene in genelist:
        returnlist.append(id_rank[id_rank['protein_name'] == gene]['string_id'].to_string(index=False))
    return returnlist, id_rank, na_idrank


def UniProt_StringToCommon(strings):
    mapping_dict = cp.load(open("data/temp/mapping_dict.txt", "rb"))
    idtable = pd.DataFrame.from_dict(mapping_dict, orient='index')
    idtable.reset_index(inplace=True)
    idtable = idtable.rename(columns={"index": "string_id", 0: "protein_name"})
    protein_counts = pd.read_csv('/home/ahphan/Downloads/protein_counts.tsv', sep="\t", header=None,
                                 names=['string_id', 'year', 'fracct'])
    protein_counts = protein_counts.astype({'year':'int64', 'fracct':'float64'})
    idtable['string_id'] = idtable['string_id'].apply(lambda x: x.split(".")[-1])
    protein_counts_ovt = protein_counts.groupby(['string_id','year']).sum().groupby(level=0).cumsum().reset_index()
    protein_counts['ovt_fracct'] = protein_counts_ovt['fracct']
    protein_counts_1121 = protein_counts.loc[(protein_counts['year'] >= 2011) & (protein_counts['year'] < 2022)].iloc[:,[0,1,3]]
    fc = pd.merge(protein_counts_1121, idtable, on='string_id', how='left').reset_index().iloc[:,[4,2,3]].dropna()
    fc.to_csv('/home/ahphan/RotationData/Friedberg/ProteinFractionCounts_1121.tsv', sep="\t", header=True,index=False)

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
        "limit": 20,
        "caller_identity": "ignoromenot"  # your app name
    }
    response = requests.post(request_url, data=params)
    print("\nInteraction partners of ignorome genes:")
    print("QueryIgnorome\tPartners\tInteraction")
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        p1, p2 = l[0], l[1]  # protein preferred name
        # filter the interaction according to coexpression score (l[9])
        coexpression_score = float(l[9])
        if coexpression_score > float(calculated_threshold)/1000:
            # print
            print("\t".join([p1, p2, "coexpression (score. %.3f)" % coexpression_score]))
    return


def main():
    import cProfile
    import pstats
    import io
    seed = 123456
    pr = cProfile.Profile()
    pr.enable()

    global options
    options = parse_commandlinearguments()
    if len(sys.argv) < 4:
        print("Please use the --help option to get usage information")
    if options.input.split("/")[-1].split(".")[-1] != "gaf":
        # Find ignorome from ranktable
        ranktable = pd.read_csv(options.input, sep="\t")
        specifier = "given-metric"
    else:
        ranktable, ia_blacklist = rank_gaf_file(options.input, options.metric, options.aspect, options.recal,
                                                swiss_no_gaf_path='./swiss-no-gaf.txt')
        specifier = options.input.split("/")[-1].split(".")[0].split("_")[1]

    most = filter_most_annotated(ranktable, options.threshold_top, options.percentile_top)
    # for tbot in [0,6,10,15]:
    #     for pppi in [99,98,95,85]:
    least = filter_least_annotated(ranktable, ia_blacklist, options.threshold_bot, options.percentile_bot)
    # coexp_threshold, ignorome_list, known_list = associated_ignorome(options.stringppi, options.threshold_ppi,
    #                                                                  options.percentile_ppi, most, least,
    #                                                                  options.idtable, ranktable)
    coexp_threshold, ignorome_list, known_list, ignorome_table, unre_table = associated_ignorome2(options.stringppi, options.threshold_ppi,
                                                                     options.percentile_ppi, most, least, ranktable)
    # bipartite_graph(most, least, ignorome_table)
    # species = options.stringppi.split("/")[-1].split(".")[0]
    # get_network_api(coexp_threshold, ignorome_list, species)
    ignorome_table.to_csv('./NewIgnoromeGenes_' + str(specifier) +'_'+ str(options.metric) + str(options.aspect) +'_string' + str(int(coexp_threshold)) + '.tsv',
                          sep="\t", header=True, index=False)
            # ignorome_table.to_csv('./IgnoromeGenes_' + str(options.metric) + str(options.aspect) +
            #                           '_t' + str(options.percentile_top) + '_b' + str(tbot) + '_string' + str(int(coexp_threshold)) + '.tsv', sep="\t", header=True, index=False)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
    ps.print_stats()
    with open('./cProfile_ign.txt', 'w') as f:
        f.write(s.getvalue())

if __name__ == "__main__":
    main()
