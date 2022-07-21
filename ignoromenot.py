# This script finds ignorome that are highly-associated with highly-annotated genes
# Ranking based on annotation count, Philip Lord information content, PubMed mentions, Jensen score;
# Association based on STRING co-expression

from __future__ import print_function
from __future__ import division
import networkx as nx
import sys
from networkx.algorithms import bipartite
from operator import itemgetter
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
report = 0

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="ignoromenot.py")
    mutex_parser_evidence = parser.add_mutually_exclusive_group()
    mutex_parser_assigned_by = parser.add_mutually_exclusive_group()
    mutex_parser_IC_THRESH = parser.add_mutually_exclusive_group()
    # mutex_parser_PL_THRESH=parser.add_mutually_exclusive_group()
    mutex_parser_select_references = parser.add_mutually_exclusive_group()
    requiredArguments = parser.add_argument_group("Required arguments")

    parser.add_argument('--prefix', '-pref', help="Add a prefix to the name of your output files.")
    parser.add_argument('--cutoff_prot', '-cprot',
                        help="The threshold level for deciding to eliminate annotations which come from references that annotate more than the given 'threshold' number of PROTEINS")
    parser.add_argument('--cutoff_attn', '-cattn',
                        help="The threshold level for deciding to eliminate annotations which come from references that annotate more than the given 'threshold' number of ANNOTATIONS")

    parser.add_argument('--output', '-odir', help="Writes the final outputs to the directory in this path.")

    mutex_parser_evidence.add_argument('--evidence', '-e', nargs="+",
                                       help="Accepts Standard Evidence Codes outlined in http://geneontology.org/page/guide-go-evidence-codes. All 3 letter code for each standard evidence is acceptable. In addition to that EXPEC is accepted which will pull out all annotations which are made experimentally. COMPEC will extract all annotations which have been done computationally. Similarly, AUTHEC and CUREC are also accepted. Cannot be provided if -einv is provided ")
    mutex_parser_evidence.add_argument('--evidence_inverse', '-einv', nargs="+",
                                       help="Leaves out the provided Evidence Codes. Cannot be provided if -e is provided")

    requiredArguments.add_argument('--input', '-i', nargs="+",
                                   help="The input file path. Please remember the name of the file must start with goa in front of it, with the name of the species following separated by an underscore",
                                   required=True)
    parser.add_argument('--aspect', '-a', nargs="+",
                        help="Enter P, C or F for Biological Process, Cellular Component or Molecular Function respectively")

    mutex_parser_assigned_by.add_argument('--assigned_by', '-assgn', nargs="+",
                                          help="Choose only those annotations which have been annotated by the provided list of databases. Cannot be provided if -assgninv is provided")
    mutex_parser_assigned_by.add_argument('--assigned_by_inverse', '-assgninv', nargs="+",
                                          help="Choose only those annotations which have NOT been annotated by the provided list of databases. Cannot be provided if -assgn is provided")
    parser.add_argument('--recalculate', '-recal',
                        help="Set this to 1 if you wish to enforce the recalculation of the Information Accretion for every GO term. Calculation of the information accretion is time consuming. Therefore keep it to zero if you are performing rerun on old data. The program will then read the information accretion values from a file which it wrote to in the previous run of the program",
                        default=0)

    mutex_parser_IC_THRESH.add_argument('--info_threshold_Wyatt_Clark_percentile', '-WCTHRESHp',
                                        help="Provide the percentile p. All annotations having information content below p will be discarded")
    mutex_parser_IC_THRESH.add_argument('--info_threshold_Wyatt_Clark', '-WCTHRESH',
                                        help="Provide a threshold value t. All annotations having information content below t will be discarded")

    mutex_parser_IC_THRESH.add_argument('--info_threshold_Phillip_Lord_percentile', '-PLTHRESHp',
                                        help="Provide the percentile p. All annotations having information content below p will be discarded")
    mutex_parser_IC_THRESH.add_argument('--info_threshold_Phillip_Lord', '-PLTHRESH',
                                        help="Provide a threshold value t. All annotations having information content below t will be discarded")

    parser.add_argument('--verbose', '-v',
                        help="Set this argument to 1 if you wish to view the outcome of each operation on console",
                        default=0)
    parser.add_argument('--date_before', '-dbfr',
                        help="The date entered here will be parsed by the parser from dateutil package. For more information on acceptable date formats please visit https://github.com/dateutil/dateutil/. All annotations made prior to this date will be picked up")
    parser.add_argument('--date_after', '-daftr',
                        help="The date entered here will be parsed by the parser from dateutil package. For more information on acceptable date formats please visit https://github.com/dateutil/dateutil/. All annotations made after this date will be picked up")
    parser.add_argument('--single_file', '-single', default=0,
                        help="Set to 1 in order to output the results of each individual species in a single file.")

    mutex_parser_select_references.add_argument('--select_references', '-selref', nargs='+',
                                                help='Provide the paths to files which contain references you wish to select. It is possible to include references in case you wish to select annotations made by a few references. This will prompt the program to interpret string which have the keywords \'GO_REF\',\'PMID\' and \'Reactome\' as a GO reference. Strings which do not contain that keyword will be interpreted as a file path which the program will except to contain a list of GO references. The program will accept a mixture of GO_REF and file names. It is also possible to choose all references of a particular category and a handful of references from another. For example if you wish to choose all PMID references, just put PMID. The program will then select all PMID references. Currently the program can accept PMID, GO_REF and Reactome')
    mutex_parser_select_references.add_argument('--select_references_inverse', '-selrefinv', nargs='+',
                                                help='Works like -selref but does not select the references which have been provided as input')
    parser.add_argument('--report', '-r',
                        help="Provide the path where the report file will be stored. If you are providing a path please make sure your path ends with a '/'. Otherwise the program will assume the last string after the final '/' as the name of the report file. A single report file will be generated. Information for each species will be put into individual worksheets.")
    parser.add_argument('--histogram', '-hist',
                        help="Set this option to 1 if you wish to view the histogram of GO_TERM frequency before and after debiasing is performed with respect to cutoffs based on number of proteins or annotations. If you wish to save the file then please enter a filepath. If you are providing a path please make sure your path ends with a '/'. Otherwise the program will assume the last string after the final '/' as the name of the image file. Separate histograms will be generated for each species.")

    args = parser.parse_args()
    return args

def vprint(*s):
    global verbose
    # print(s,verbose)
    if verbose == 1:
        for string in s:
            print(string, end="")
        print()

def filter_most_coexp(filename, percentile):
    all_interaction = pd.read_csv(filename, sep=",", header=0, usecols=['protein1', "protein2", "coexpression"])
    all_coexp = all_interaction[all_interaction['coexpression'] > 0]
    MAX_COEXP = max(all_coexp['coexpression'])
    if percentile != None:
        top_coexp = all_coexp[all_coexp['coexpression'] > percentile*(1-percentile/100)*MAX_COEXP]
    return top_coexp

def filter_most_annotated(ranktable, idmatch, threshold):
    rankIC = pd.read_csv(ranktable, sep= "\t", header=0, names= ["alias", "C", "F", "P"])
    idtable = pd.read_csv(idmatch, sep = "\t", header=0, names=["string_id", "alias", "source"], usecols=["string_id", "alias"])
    rankIC_stringID = pd.merge(rankIC, idtable, on = "alias", how="left")

    most_annt_dict = dict()
    for row in rankIC_stringID:
        if row["C"] > threshold or row["F"] > threshold or row["P"] > threshold:
            most_annt_dict[row["alias"]] = {row["C"], row["F"], row["P"]}

    return most_annt_dict
def filter_least_annotated(ranktable, idmatch, threshold):
    rankIC = pd.read_csv(ranktable, sep="\t", header=0)

    return least_annt_dict

def associated_ignorome(top_coexp, ):

def get_coexpression(coexpression_threshold, gene_all_common):
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
    # print("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\tnscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore")
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        p1, p2= l[2], l[3]  # protein preferred name
        # filter the interaction according to coexpression score (l[9])
        coexpression_score = float(l[10])
        if coexpression_score > coexpression_threshold:
            ## print
            print("\t".join([p1, p2, "coexpression (score. %.3f)" % coexpression_score]))

gene_list = ["ENSP00000395234"]
get_coexpression(0, gene_list)

def main():
    global verbose, options, report
    commandLineArg = sys.argv
    if len(commandLineArg) == 1:
        print("Please use the --help option to get usage information")
    # Parse command line arguments
    options = parseCommandLineArguments()

if __name__ == "main":
    main()