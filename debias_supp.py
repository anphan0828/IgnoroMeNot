#!/usr/bin/
"""
"""
from __future__ import print_function
from __future__ import division
import networkx as nx
import sys

import pandas as pd
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

# DATADIR = "data_for_debias_0914/"
# # Some filenames
# FILE_ALTERNATE_ID_TO_ID_MAPPING = DATADIR + "alt_to_id.graph"
# FILE_CAFA_ID_TO_UNIPROT_ID_MAP = DATADIR + "CAFAIDTOUniprotIDMap.txt"
# FILE_MFO_ONTOLOGY_GRAPH = DATADIR + "mf.graph"
# FILE_BPO_ONTOLOGY_GRAPH = DATADIR + "bp.graph"
# FILE_CCO_ONTOLOGY_GRAPH = DATADIR + "cc.graph"
# FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH = DATADIR + "mf_ancestors.map"
# FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH = DATADIR + "bp_ancestors.map"
# FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH = DATADIR + "cc_ancestors.map"

verbose = 0
options = ""
report = 0
GAF21FIELDS = [
    'DB',
    'DB_Object_ID',
    'DB_Object_Symbol',
    'Qualifier',
    'GO_ID',
    'DB:Reference',
    'Evidence',
    'With',
    'Aspect',
    'DB_Object_Name',
    'Synonym',
    'DB_Object_Type',
    'Taxon_ID',
    'Date',
    'Assigned_By',
    'Annotation_Extension',
    'Gene_Product_Form_ID'
]

EXPEC = [
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP"
]

COMPEC = [
    "ISS",
    "ISO",
    "ISA",
    "ISM",
    "IGC",
    "IBA",
    "IBD",
    "IKR",
    "IRD",
    "RCA"
]

AUTHEC = [
    "TAS",
    "NAS"
]

CUREC = [
    "IC",
    "ND"
]

IEA = ["IEA"]


def column(matrix, i):
    f = itemgetter(i)
    return map(f, matrix)


def chooseProteinsBasedOnPublications(data, cutoff_prot, cutoff_attn):
    """
    This function will read the set of proteins and will choose only those proteins which have been probed in those publications which
    deal with less than -d <number> of proteins
    This program creates a bipartite graph with one set as the GO terms and the other set as the references and cross links them with the 
    GO_TERMS as weights to the edges.
    This function can be used to select a cut off based on number of annotations by a particular reference or even number of proteins
    annotated by a reference. It is recommended that the protein cut-off, i.e. -cprot, be used instead of the annotations cutoff. Since
    it is relevant for a reference to provide more annotations to fewer proteins than to work with a lot of proteins.
    """
    mapping = []
    for attnid in data:
        per_annotation = data[attnid]
        go = per_annotation['GO_ID']  # Extracting the Gene Ontology
        protein = per_annotation['DB'] + "_" + per_annotation['DB_Object_ID']
        ref = per_annotation['DB:Reference']  # The reference
        mapping.append([protein, ref, go, attnid])  # Appending the annotation id for later identification
    g = nx.MultiGraph()
    g.add_nodes_from(column(mapping, 0), bipartite=0)
    g.add_nodes_from(column(mapping, 1), bipartite=1)
    for triplet in mapping:
        g.add_edge(triplet[0], triplet[1], weight=triplet[2] + "_" + triplet[3])

    simple_g = nx.Graph(g)  # Converting the multi graph to a simple graph without parallel edges

    no_of_prot_annotations_by_each_ref = []
    for ref in list(set(column(mapping, 1))):
        no_of_prot_annotations_by_each_ref.append(simple_g.degree(ref))

    if cutoff_attn == None:
        graph = simple_g
        threshold = int(cutoff_prot)
    elif cutoff_prot == None:
        graph = g
        threshold = int(cutoff_attn)
    list_of_chosen_attn = []
    list_of_discarded_attn = []
    # Looping through each GO term and selecting those for which there is at least one reference which probes fewer proteins than threshold
    for protein in list(set(column(mapping, 0))):
        references = g.neighbors(protein)
        for ref in references:
            # Condition for inclusion
            if graph.degree(ref) <= threshold:
                for key in g.get_edge_data(protein, ref):
                    weight = g.get_edge_data(protein, ref)[key]
                    list_of_chosen_attn.append(weight['weight'].split("_")[1])
            else:
                for key in g.get_edge_data(protein, ref):
                    weight = g.get_edge_data(protein, ref)[key]
                    list_of_discarded_attn.append(weight['weight'].split("_")[1])
    new_data = dict()
    discarded_data = dict()
    for attid in list_of_chosen_attn:
        new_data[attid] = data[attid]
    for nid in list_of_discarded_attn:
        discarded_data[nid] = data[nid]
    return new_data, discarded_data


#

def convertToDictionary(filename):
    """
    This function reads from the input gaf file and converts it to a dictionary. This function is deprecated and will be removed in further releases.
    Instead of using this function the program now makes use of the gaf iterator function from biopython.
    """
    alt_id_to_id_map = cp.load(open(FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb"))
    fhr = open(filename, "r")
    data = dict()
    counter = 1
    for line in fhr:
        if "!" not in line:
            line = line.split("\t")
            id = "anntn" + str(counter)
            per_annotation = dict()
            for f_no, field in enumerate(GAF21FIELDS):
                if field == "GO_ID":
                    if line[f_no].strip() in alt_id_to_id_map:
                        # print(line[f_no].strip())
                        line[f_no] = alt_id_to_id_map[line[f_no].strip()]
                per_annotation[field] = line[f_no]
            data[id] = per_annotation
            counter += 1
            """if(len(data)==10):
                break"""
    fhr.close()
    return data


# def convertFromGAFToRequiredFormat(gaf):  # commented on Sep 14
#     """
#     This function takes the data input which is created by gaf iterator and then makes few changes
#     in the annotations which is relevant to this program.
#     """
#     alt_id_to_id_map = cp.load(open(FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb"))
#     counter = 1
#     data = dict()
#     for annotation in gaf:
#         id = "anntn" + str(counter)
#         if annotation['GO_ID'] in alt_id_to_id_map:
#             annotation['GO_ID'] = alt_id_to_id_map[annotation['GO_ID']]
#         annotation['DB:Reference'] = annotation['DB:Reference'][0]
#         annotation['Date'] = parser.parse(annotation['Date']).date()
#         # annotation['Qualifier']='|'.join(annotation['Qualifier'])
#         # print(annotation['Evidence'])
#         data[id] = annotation
#         counter += 1
#     return data

def convertFromGAFToRequiredFormat(gaf):  # newly rewritten Sep 14
    """
    This function takes the data input which is created by gaf iterator and then makes few changes
    in the annotations which is relevant to this program.
    """
    alt_id_to_id_map = cp.load(open(FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb"))
    #print(alt_id_to_id_map)
    counter = 1
    data = dict()
    for annotation in gaf:
        id = "anntn" + str(counter)
        # if annotation['GO_ID'] == 'GO:0002506':
        #     print(annotation['GO_ID'], annotation['Aspect'])
        if annotation['GO_ID'] in alt_id_to_id_map:
            #print(annotation['GO_ID'], alt_id_to_id_map[annotation['GO_ID']][0], alt_id_to_id_map[annotation['GO_ID']][1])
            if alt_id_to_id_map[annotation['GO_ID']][1] == 'mf':
                annotation['Aspect'] = "F"
            elif alt_id_to_id_map[annotation['GO_ID']][1] == 'bp':
                annotation['Aspect'] = "P"
            elif alt_id_to_id_map[annotation['GO_ID']][1] == 'cc':
                annotation['Aspect'] = "C"
            annotation['GO_ID'] = alt_id_to_id_map[annotation['GO_ID']][0]
        # if annotation['GO_ID'] == 'GO:0030247':
        #     print("again", annotation['GO_ID'], annotation['Aspect'])
        annotation['DB:Reference'] = annotation['DB:Reference'][0]
        annotation['Date'] = parser.parse(annotation['Date']).date()
        # annotation['Qualifier']='|'.join(annotation['Qualifier'])
        # print(annotation['Evidence'])
        data[id] = annotation
        counter += 1
    return data

def writeToFile(data, filename, input_filename):
    """
    This function will write the content of the data structure 'data' to the output file.
    It requires the input file to read the header. Inclusion of the header is mandatory.
    """
    vprint("Writing to file ", filename)
    # print(filename)
    filepath = "/".join(filename.split("/")[:-1])
    try:
        if os.path.isdir(filepath) == False:
            os.makedirs(filepath)
    except OSError:
        print(
            "You do not have sufficient Permissions to create the folder. Please alter the permissions or provide a different path.")
        sys.exit()
    fhr = open(input_filename, "r")
    header = ""
    for line in fhr:
        if line[0] == '!':
            header += line
    fhr.close()
    fhw = open(filename + ".gaf", "w")
    fhw.write(header)
    for key in data:
        per_annotation = data[key]
        per_annotation['Qualifier'] = '|'.join(per_annotation['Qualifier'])
        per_annotation['With'] = '|'.join(per_annotation['With'])
        per_annotation['Synonym'] = '|'.join(per_annotation['Synonym'])
        per_annotation['Taxon_ID'] = '|'.join(per_annotation['Taxon_ID'])
        per_annotation['Date'] = ''.join(str(per_annotation['Date']).split("-"))
        # vprint(per_annotation)
        string = ""
        for field in GAF21FIELDS:
            try:
                string += per_annotation[field] + "\t"
            except TypeError:
                print("Exception has occurred in function writeToFile")
                print(per_annotation)
                print(field)
                print(per_annotation[field])
                exit()
        string += '\n'
        fhw.write(string)
    fhw.close()


def checkEvidenceCodeForCorrectness(codes):
    """
    This function checks whether the Evidence Codes provided by the user.
    It will return false if any incorrect evidence code is provided.
    """
    for evidence in codes:
        if (
                evidence != "COMPEC" and evidence != "EXPEC" and evidence != "AUTHEC" and evidence != "CUREC" and evidence != "IEA"):
            evidence = [evidence]
            if True not in set([set(evidence).issubset(set(COMPEC)), set(evidence).issubset(set(EXPEC)),
                                set(evidence).issubset(set(CUREC)), set(evidence).issubset(set(IEA))]):
                return False
    return True


def chooseProteinsBasedOnEvidenceCodes(data, evidence_list, evidence_inverse_list):
    """
    This function will select only those annotations which have been annotated by the provided Evidence Codes
    """
    # Checking whether the provided Evidence codes are correct or not
    if evidence_list is not None:
        if checkEvidenceCodeForCorrectness(evidence_list) == False:
            vprint(
                "Invalid arguments for Evidence Codes provided please check http://geneontology.org/page/guide-go-evidence-codes")
            sys.exit()
    else:
        if checkEvidenceCodeForCorrectness(evidence_inverse_list) == False:
            vprint(
                "Invalid arguments for Evidence Codes provided please check http://geneontology.org/page/guide-go-evidence-codes")
            sys.exit()

    select_these_EC = []
    if (evidence_list is not None):
        for evidence in evidence_list:
            if (
                    evidence != "COMPEC" and evidence != "EXPEC" and evidence != "AUTHEC" and evidence != "CUREC" and evidence != "IEA"):
                select_these_EC.append(evidence)
            else:
                EC_set = ""
                if (evidence == "COMPEC"):
                    EC_set = COMPEC
                elif (evidence == "EXPEC"):
                    EC_set = EXPEC
                elif (evidence == "AUTHEC"):
                    EC_set = AUTHEC
                elif (evidence == "CUREC"):
                    EC_set = CUREC
                elif (evidence == "IEA"):
                    EC_set = IEA
                for ele in EC_set:
                    select_these_EC.append(ele)
    else:
        select_these_EC.extend(COMPEC)
        select_these_EC.extend(EXPEC)
        select_these_EC.extend(AUTHEC)
        select_these_EC.extend(CUREC)
        select_these_EC.extend(IEA)

        for evidence in evidence_inverse_list:
            if (
                    evidence != "COMPEC" and evidence != "EXPEC" and evidence != "AUTHEC" and evidence != "CUREC" and evidence != "IEA"):
                select_these_EC.remove(evidence)
            else:
                EC_set = ""
                if (evidence == "COMPEC"):
                    EC_set = COMPEC
                elif (evidence == "EXPEC"):
                    EC_set = EXPEC
                elif (evidence == "AUTHEC"):
                    EC_set = AUTHEC
                elif (evidence == "CUREC"):
                    EC_set = CUREC
                elif (evidence == "IEA"):
                    EC_set = IEA
                for ele in EC_set:
                    select_these_EC.remove(ele)

    new_data = dict()
    vprint(select_these_EC)
    for attnid in data:
        per_annotation = data[attnid]
        if per_annotation['Evidence'] in select_these_EC:
            new_data[attnid] = per_annotation
    return new_data


# This function reads the data entered via command line and returns a dictionary with all relevant options

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="debias.py")
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


def createProteinToGOMapping(data):
    """
    This function creates a dictionary where key is a protein. Each protein refers to a list where the list consists of GO_TERMS.
    """
    prot_to_go = dict()
    all_GO = []
    alt_id_to_id_map = cp.load(open(FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb"))
    for attnid in data:
        annotation = data[attnid]
        prot_id = annotation['DB'] + '_' + annotation['DB_Object_ID']
        GO_term = annotation['GO_ID']
        if GO_term in alt_id_to_id_map:
            GO_term = alt_id_to_id_map[GO_term]
        all_GO.append(GO_term)
        if prot_id not in prot_to_go:
            prot_to_go[prot_id] = []
            if [GO_term, annotation['Aspect']] not in prot_to_go[prot_id]:
                prot_to_go[prot_id].append([GO_term, annotation['Aspect']])
        else:
            if [GO_term, annotation['Aspect']] not in prot_to_go[prot_id]:
                prot_to_go[prot_id].append([GO_term, annotation['Aspect']])
        # vprint(prot_to_go[prot_id])
    return prot_to_go, list(set(all_GO))


def propagateOntologies(Prot_to_GO_Map):
    """
    This function takes in each annotation and constructs the ancestors of that term from their respective Aspect
    """
    mf_g = cp.load(open(FILE_MFO_ONTOLOGY_GRAPH, "rb"))
    bp_g = cp.load(open(FILE_BPO_ONTOLOGY_GRAPH, "rb"))
    cc_g = cp.load(open(FILE_CCO_ONTOLOGY_GRAPH, "rb"))
    alt_id_to_id_map = cp.load(open(FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb"))
    # vprint(alt_id_to_id_map)
    Prot_to_GO_Map_new = dict()
    mf_ancestors = cp.load(open(FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH, "rb"))
    bp_ancestors = cp.load(open(FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH, "rb"))
    cc_ancestors = cp.load(open(FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH, "rb"))
    for eachprotein in Prot_to_GO_Map:
        ancestors = []
        annotations = Prot_to_GO_Map[eachprotein]
        for annotation in annotations:
            aspect = annotation[1]
            GO_term = annotation[0]
            # print(aspect, GO_term)
            if aspect == 'F':
                ancestors.extend(mf_ancestors[GO_term])
            if aspect == 'P':
                ancestors.extend(bp_ancestors[GO_term])
            if aspect == 'C':
                ancestors.extend(cc_ancestors[GO_term])
        ancestors = list(set(ancestors))
        Prot_to_GO_Map_new[eachprotein] = ancestors
    return Prot_to_GO_Map_new


def findFrequency(annotations, Prot_to_GO_Map):
    count = 0
    if annotations == None:
        return 0
    for prot in Prot_to_GO_Map:
        if set(annotations).issubset(set(Prot_to_GO_Map[prot])):
            count += 1
    return count


def assignProbabilitiesToOntologyTree(g, Prot_to_GO_Map, all_GO_Terms, ontology_to_ia_map, aspect):  # conditional prob of each GO term
    fh = open("WC_IC_"+str(aspect)+".txt", "w")
    go_to_ia_dict = {}  # return a dict for each aspect {GO_ID: IA}
    for node_num, node in enumerate(g.nodes()):
        if (node not in all_GO_Terms):
            ontology_to_ia_map[node] = [0, 0]
            continue
        if node_num % 100 == 0:
            vprint(node_num, " proteins processed for ", aspect)
        predecessor = g.successors(node)
        # vprint(node,predecessor)
        predecessor_with_node = []
        predecessor_with_node.extend(predecessor)
        predecessor_with_node.append(node)
        denom = findFrequency(predecessor, Prot_to_GO_Map)
        num = findFrequency(predecessor_with_node, Prot_to_GO_Map)
        # vprint(node,g.successors(node))
        """vprint(predecessor_with_node,num)
        vprint(predecessor,denom)"""
        if (denom == 0):
            prob = 0
        else:
            prob = num / denom
        ontology_to_ia_map[node] = [prob, -math.log(prob, 2)]
        # print(node, num, denom)
        fh.write(str(node) + "\t" + str(num) + "\t" + str(denom) + "\t" + str(prob) + "\t" + str(-math.log(prob, 2)) + "\n")
        go_to_ia_dict[node] = -math.log(prob,2) # dictionary of IA per term
    fh.close()
    return go_to_ia_dict


def assignProbabilitiesToOntologyGraphs(Prot_to_GO_Map, all_GO_Terms):  # first step: graph before tree,
    # deleted last argument aspects because not used
    mf_g = cp.load(open(FILE_MFO_ONTOLOGY_GRAPH, "rb"))
    bp_g = cp.load(open(FILE_BPO_ONTOLOGY_GRAPH, "rb"))
    cc_g = cp.load(open(FILE_CCO_ONTOLOGY_GRAPH, "rb"))
    ontology_to_ia_map = dict()
    go_to_ia_dict_f = assignProbabilitiesToOntologyTree(mf_g, Prot_to_GO_Map, all_GO_Terms, ontology_to_ia_map, 'MFO')
    # fdict = open("WC_dictf", "w")
    # for GO in go_to_ia_dict_f.keys():
    #     fdict.write(GO + "\t" + str(go_to_ia_dict_f[GO]) + "\n")
    # fdict.close()
    go_to_ia_dict_p = assignProbabilitiesToOntologyTree(bp_g, Prot_to_GO_Map, all_GO_Terms, ontology_to_ia_map, 'BPO')
    # fdict = open("WC_dictp", "w")
    # for GO in go_to_ia_dict_p.keys():
    #     fdict.write(GO + "\t" + str(go_to_ia_dict_p[GO]) + "\n")
    # fdict.close()
    go_to_ia_dict_c = assignProbabilitiesToOntologyTree(cc_g, Prot_to_GO_Map, all_GO_Terms, ontology_to_ia_map, 'CCO')
    # fdict = open("WC_dictc", "w")
    # for GO in go_to_ia_dict_c.keys():
    #     fdict.write(GO + "\t" + str(go_to_ia_dict_c[GO]) + "\n")
    # fdict.close()

    """for GO in ontology_to_ia_map:
        vprint(ontology_to_ia_map[GO])"""
    fileTemp1 = open("IC.txt", "w")
    for GO in ontology_to_ia_map:
        fileTemp1.write(GO + "\t")
        # print(ontology_to_ia_map[GO])
        # print(ontology_to_ia_map[GO][1])
        fileTemp1.write(str(ontology_to_ia_map[GO][1]) + "\n")
    fileTemp1.close()
    fileTemp1 = open("FreqIC.txt", "w")
    for GO in ontology_to_ia_map:
        fileTemp1.write(GO + "\t" + str(ontology_to_ia_map[GO][0]) + "\t")
        # print(ontology_to_ia_map[GO])
        # print(ontology_to_ia_map[GO][1])
        fileTemp1.write(str(ontology_to_ia_map[GO][1]) + "\n")
    fileTemp1.close()
    return ontology_to_ia_map, go_to_ia_dict_f, go_to_ia_dict_p, go_to_ia_dict_c


def calculateInformationAccretionForEachProtein(Prot_to_GO_Map, ontology_to_ia_map):
    vprint("Starting calculation of ia")
    infoAccr = dict()
    alt_id_to_id_map = cp.load(open(FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb"))
    for prot in Prot_to_GO_Map:
        annotations = Prot_to_GO_Map[prot]
        ia = 0
        for annotation in annotations:
            if len(annotation) == 2:
                GO_term = annotation[0]
            else:
                GO_term = annotation
            if GO_term not in ontology_to_ia_map:
                GO_term = alt_id_to_id_map[GO_term]
            # vprint(prot,annotation[0])
            ia += ontology_to_ia_map[GO_term][1]
        infoAccr[prot] = ia
    print(len(infoAccr))
    return infoAccr # dictionary of IA per protein


def chooseGOBasedOnAspect(data, aspect):
    new_data = dict()
    # vprint( aspect )
    for attnid in data:
        if data[attnid]['Aspect'] in aspect:
            new_data[attnid] = data[attnid]
    return new_data


def chooseGOBasedOnAssignedBy(data, assigned_by, assigned_by_inverse):
    new_data = dict()
    for attnid in data:
        if (assigned_by != None):
            if data[attnid]['Assigned_By'] in assigned_by:
                new_data[attnid] = data[attnid]
        else:
            if data[attnid]['Assigned_By'] not in assigned_by_inverse:
                new_data[attnid] = data[attnid]
    return new_data


def calculatePhillipLordInformationContent(data, crisp, percentile_val):
    go_terms = []
    """Prot_to_GO_Map, all_GO_Terms_in_corpus = createProteinToGOMapping( data )
    Prot_to_GO_Map_propagated = propagateOntologies( Prot_to_GO_Map )"""
    # alt_id_to_id_map = cp.load( open( FILE_ALTERNATE_ID_TO_ID_MAPPING, "rb" ) )
    """for eachprot in Prot_to_GO_Map:
        go_terms.extend([annotation[0] for annotation in Prot_to_GO_Map[eachprot]])"""
    """for eachprot in Prot_to_GO_Map_propagated:
        go_terms.extend(Prot_to_GO_Map_propagated[eachprot])"""

    for attnid in data:
        go_terms.append(data[attnid]["GO_ID"])
    GO_term_to_PL_info = collections.Counter(go_terms)

    ic = []
    for term in GO_term_to_PL_info:
        # vprint(term,x[term],x[term]/len(go_terms))
        GO_term_to_PL_info[term] = -math.log(GO_term_to_PL_info[term] / len(go_terms), 2)
        ic.append(GO_term_to_PL_info[term])
        # print(term, GO_term_to_PL_info[term], len(go_terms), -math.log(GO_term_to_PL_info[term]/len(go_terms),2))

    if crisp == None:
        threshold = ((max(ic) - min(ic)) * float(percentile_val) / 100) + min(ic)
    else:
        threshold = float(crisp)
        num_bins = 10
        # the histogram of the data
        # n, bins, patches = plt.hist(ic, num_bins, facecolor='green', alpha=0.9)
        # plt.show()

    print("The maximum value of information content is ", max(ic))
    print("The minimum value of information content is ", min(ic))
    print("The chosen threshold is ", threshold)
    new_data = dict()
    for attnid in data:
        annotation = data[attnid]
        if GO_term_to_PL_info[annotation["GO_ID"]] >= threshold:
            new_data[attnid] = data[attnid]
    return new_data
    # vprint(collections.Counter(go_terms))


def calculateWyattClarkInformationContent(data, recal, crisp, percentile_val, outputfiles, input_num):
    """
    This function will display some essential statistics when the value of threshold
    is crisp and a percentile is not provided.
    """
    # vprint(outputfiles[0].split("_"))
    ontology_to_ia_map_filename = "ontology_to_ia_map_" + "_".join(
        outputfiles[input_num].split("/")[-1].split("_")[:-2]) + ".txt"
    goia_f_filename = "goia_f"+"_".join(
        outputfiles[input_num].split("/")[-1].split("_")[:-2]) + ".txt"
    goia_p_filename = "goia_p" + "_".join(
        outputfiles[input_num].split("/")[-1].split("_")[:-2]) + ".txt"
    goia_c_filename = "goia_c" + "_".join(
        outputfiles[input_num].split("/")[-1].split("_")[:-2]) + ".txt"
    # vprint(ontology_to_ia_map_filename)
    # exit()
    Prot_to_GO_Map, all_GO_Terms_in_corpus = createProteinToGOMapping(data)
    Prot_to_GO_Map_propagated = propagateOntologies(Prot_to_GO_Map)
    if (recal == 1):
        vprint(
            "Recalculating Information Accretion for Wyatt Clark Information Content. This may take a long time depending on the size of input")
        ontology_to_ia_map, goia_f, goia_p, goia_c = assignProbabilitiesToOntologyGraphs(Prot_to_GO_Map_propagated, all_GO_Terms_in_corpus)  # 3 dicts returned from Graph function call will be used in generateHistogram
        if os.path.isdir("data/temp/") == False:
            os.makedirs("data/temp/")
        cp.dump(ontology_to_ia_map, open("data/temp/" + ontology_to_ia_map_filename, "wb"))  #TODO: dump 3 goia dicts, too
        cp.dump(goia_f, open("data/temp/"+ goia_f_filename, "wb"))
        cp.dump(goia_p, open("data/temp/"+ goia_p_filename, "wb"))
        cp.dump(goia_c, open("data/temp/"+ goia_c_filename, "wb"))
    else:
        vprint("Skipping recalculation of Information Accretion for Wyatt Clark")
    try:
        ontology_to_ia_map = cp.load(open("data/temp/" + ontology_to_ia_map_filename, "rb")) #TODO: load 3 goia dicts
        goia_f = cp.load(open("data/temp/" + ontology_to_ia_map_filename, "rb"))
        goia_p = cp.load(open("data/temp/" + ontology_to_ia_map_filename, "rb"))
        goia_c = cp.load(open("data/temp/" + ontology_to_ia_map_filename, "rb"))
    except IOError as e:
        print("File for GO_Term to ia NOT FOUND. Please rerun the program with the argument -recal 1")
        exit()
    protInfoAccretion = calculateInformationAccretionForEachProtein(Prot_to_GO_Map_propagated, ontology_to_ia_map )
    print("Printing Wyatt Clark IA per protein: ")
    fhwc = open("WyattClarkIC-perprotein.tsv", "w") # just for IgnoroMeNot use
    fhwc.write("DB_Object_Symbol\tWyattClarkIC\n")
    for prot in protInfoAccretion:  # a dict of IA per protein
        # print(protInfoAccretion[prot], ia)
        prot_noDB = prot.split("_")[-1]
        fhwc.write(prot_noDB + "\t" + str(protInfoAccretion[prot]) + "\n")
    ia = []
    for mapping in ontology_to_ia_map:
        if ontology_to_ia_map[mapping][0] != 0:
            ia.append(ontology_to_ia_map[mapping][1])

    vprint(sorted(ia))
    vprint("Total proteins:" + str(len(ia)))
    # Doing Some statistical analysis with the distribution of information content

    if crisp == None:
        threshold = (max(ia) - min(ia)) * float(percentile_val) / 100 + min(ia)
    else:
        threshold = float(crisp)
    # vprint("Wyatt Clark Threshold",threshold,min(ia),max(ia))
    new_data = dict()

    if crisp is not None:
        num_bins = 10
        # the histogram of the data
        # n, bins, patches = plt.hist(ia, num_bins, facecolor='green', alpha=0.9)
        # plt.show()
    print("Printing Wyatt Clark IC per term:") # not doing filter to save time, term not in later.txt (after -cprot 100 is not included -> KeyError when loop over dict)
    for attnid in data:
        annotation = data[attnid]
        print(annotation["GO_ID"], ontology_to_ia_map[annotation["GO_ID"]][1])
        if ontology_to_ia_map[annotation["GO_ID"]][1] >= threshold:
            new_data[attnid] = data[attnid]

    # print("Printing Wyatt Clark IC per protein: ")
    # print(calculateInformationAccretionForEachProtein(Prot_to_GO_Map_propagated, ontology_to_ia_map))
    # vprint(threshold)
    return new_data, goia_f, goia_p, goia_c


def chooseProteinsBasedOnReferences(data, select, inverse_select):
    group = []
    references = []
    if select is not None:
        ptr = select
    else:
        ptr = inverse_select

    for item in ptr:
        vprint(item)
        if item in ("GO_REF", "PMID", "Reactome"):
            group.append(item)
        elif "GO_REF" in item or "PMID" in item or "Reactome" in item:
            references.append(item)
        else:
            for line in open(item, "r"):
                references.append(line.strip())

    new_data = dict()
    vprint(group)
    vprint(references)
    for attnid in data:
        for item in group:
            if item in data[attnid]['DB:Reference']:
                new_data[attnid] = data[attnid]
        if data[attnid]['DB:Reference'] in references:
            new_data[attnid] = data[attnid]

    if inverse_select is not None:
        newer_data = dict()
        for key in set(data.keys()) - set(new_data.keys()):
            newer_data[key] = data[key]
        # vprint(key)
        # vprint(newer_data[key])
        new_data = newer_data

    return new_data


def vprint(*s):
    global verbose
    # print(s,verbose)
    if verbose == 1:
        for string in s:
            print(string, end="")
        print()
#
#
# def printDetailsAboutData(data):
#     print("Total number of annotations in the provided Database ", len(data))
#     prots = []
#     ref = []
#     for attnid in data:
#         annotation = data[attnid]
#         prots.append(annotation['DB'] + "_" + annotation['DB_Object_ID'])
#         ref.append(annotation['DB:Reference'])
#     print("Total number of unique proteins in the provided Database ", len(set(prots)))
#     print("Total number of unique references in the provided Database ", len(set(ref)))
#
#
def chooseAnnotationsBasedOnDate(data, before, after):
    if before != None:
        before = parser.parse(before).date()
    if after != None:
        after = parser.parse(after).date()
    new_data = dict()
    for attnid in data:
        annotation = data[attnid]
        if before != None and after != None:
            if annotation['Date'] <= before and annotation['Date'] >= after:
                new_data[attnid] = annotation
        elif before != None:
            if annotation['Date'] <= before:
                new_data[attnid] = annotation
        elif after != None:
            if annotation['Date'] >= after:
                new_data[attnid] = annotation
    return new_data


def changeNameofOutputFiles(options):
    longAspect = {'P': 'BPO', 'C': 'CCO', 'F': 'MFO'}
    if options.prefix is not None:
        prefix = options.prefix
    else:
        prefix = ""
    if options.output == None:
        options.output = []
        path = "./" + prefix
    else:
        path = options.output + "/" + prefix
        options.output = []

    # vprint("Output Options ",options.output)
    for num, inputfile in enumerate(options.input):
        final_outputfilename = ""
        # vprint(inputfile)

        vprint(options.output)
        file = inputfile.split("/")[-1]
        species = file.split(".gaf")[0].split("_")[1]
        vprint("Species: " + species)
        final_outputfilename = path + species
        aspect = ""
        if options.aspect:
            for i in range(len(options.aspect)):
                aspect += longAspect[options.aspect[i]] + "_"
            final_outputfilename += "_" + aspect[:-1]
        if options.cutoff_prot:
            final_outputfilename += '_REF_' + options.cutoff_prot
        if options.evidence:
            final_outputfilename += "_" + "_".join(options.evidence)
        if options.info_threshold_Phillip_Lord:
            final_outputfilename += '_PL_' + options.info_threshold_Phillip_Lord
        elif options.info_threshold_Phillip_Lord_percentile:
            final_outputfilename += '_PLP_' + options.info_threshold_Phillip_Lord_percentile
        if options.info_threshold_Wyatt_Clark:
            final_outputfilename += '_WC_' + options.info_threshold_Wyatt_Clark
        elif options.info_threshold_Wyatt_Clark_percentile:
            final_outputfilename += '_WCP_' + options.info_threshold_Wyatt_Clark_percentile
        options.output.append(final_outputfilename)
        vprint(options.output[num])
        vprint()
        # exit()
    # vprint(options.output)
    return options.output


def combineOutputFiles(outputfiles, options):
    # print(outputfiles)
    path = "/".join(outputfiles[0].split("/")[:-1])
    file = outputfiles[0].split("/")[-1]
    """if ("./" in outputfiles[0]):
        finaloutputfilename="all_"+"_".join(file.split("_")[1:])
    else:"""
    finaloutputfilename = "all_" + "_".join(outputfiles[0].split("/")[-1].split("_")[1:])
    if options.prefix != None:
        finaloutputfilename = options.prefix + finaloutputfilename
    finaloutputfilename = path + "/" + finaloutputfilename
    print(finaloutputfilename)
    # Combine the gaf files
    header = ""
    for line in open(outputfiles[0] + ".gaf", "r"):
        if "!" in line:
            header += line
    d = ""
    for filename in outputfiles:
        for line in open(filename + ".gaf", "r"):
            if "!" not in line:
                d += line
    open(finaloutputfilename + ".gaf", "w").write(header + d)


def deleteTemporaryFiles(options):
    """
    This function deletes all the files for each organism
    """
    print("Inside delete temporary files")
    for filename in options.output:
        os.remove(filename + ".gaf")


def createReportFile(filepath, outputfilename):
    # print("Outputfilename ",outputfilename)
    # print("Filepath ",filepath)
    repfilepath = ""
    if (filepath[-1] == '/'):
        if os.path.isdir(filepath) == False:
            os.makedirs(filepath)
        filepath += "report_" + "_".join(outputfilename.split("/")[-1].split("_")[1:]) + ".xlsx"
    elif (filepath == "."):
        filepath = "report_" + "_".join(outputfilename.split("/")[-1].split("_")[1:]) + ".xlsx"
    else:
        if "/" in filepath and os.path.isdir('/'.join(filepath.split("/")[:-1])) == False:
            os.makedirs('/'.join(filepath.split("/")[:-1]))
        if ("." in filepath.split("/")[-1]):
            filepath = ".".join(filepath.split(".")[:-1]) + ".xlsx"
        else:
            filepath += ".xlsx"
    # print("Report Filepath ",filepath)
    return filepath


def countProteins(data):
    allprots = []
    for annotation in data:
        allprots.append(data[annotation]['DB'] + "_" + data[annotation]['DB_Object_ID'])
    return len(set(allprots))


def writeReport(filename, report):
    # fhw=open(filename,"w")
    all_filenames = []
    if ("/" not in filename):
        for species in report:
            all_filenames.append(species + "_" + filename[:-4] + "tsv")
    else:
        for species in report:
            all_filenames.append(
                "/".join(filename.split("/")[:-1]) + "/" + species + "_" + filename.split("/")[1][:-4] + "tsv")
    # print(all_filenames)

    # print(report)
    for sp_num, species in enumerate(report):
        fhw = open(all_filenames[sp_num], "w")
        for element in report[species]:
            # print(element)
            # fhw.write("\t".join(element)+"\n")
            for col, ele in enumerate(element):
                # print(species,element,ele)
                fhw.write(str(ele) + "\t")
            fhw.write("\n")
        # print(chr(column),chr(column+1))

    """workbook = xlsxwriter.Workbook(filename)
    bold = workbook.add_format({'bold': True,'text_wrap':True,'font_name':'monaco'})
    font = workbook.add_format({'font_name':'monaco'})
    for species in report:
        worksheet = workbook.add_worksheet(species.split("/")[-1])
        column=0
        row=0
        worksheet.write(row,column,"Operation",bold)
        worksheet.set_column(column,column, 30)
        worksheet.write(row,column+1,"Number of proteins before",bold)
        worksheet.write(row,column+2,"Number of proteins after",bold)
        worksheet.write(row,column+3,"Number of annotations before",bold)
        worksheet.write(row,column+4,"Number of annotations after",bold)
        row+=1
        #print(report[species])
        for element in report[species]:
            for col,ele in enumerate(element):
                if col==0:
                    worksheet.write(row,column+col,ele,bold)
                else:
                    worksheet.write(row,column+col,ele,font)
            row+=1
        #print(chr(column),chr(column+1))
    workbook.close()"""


def generateHistogram(data, originalData, discarded_data, species, prev, lat, go_to_ia_dict_f, go_to_ia_dict_p, go_to_ia_dict_c):
    # first argument options is removed and files are not written
    proteinDict = {}
    newproteinDict = {}
    completelyRemoved = []
    discardedDict = {}
    debiasedDict = {}
    discardedDict_F = {}
    discardedDict_P = {}
    discardedDict_C = {}
    # print(originalData)

    for new_keys in discarded_data.keys():
        if discarded_data[new_keys]['DB_Object_Symbol'] in discardedDict:
            discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][0] = discardedDict[
                                                                                 discarded_data[new_keys][
                                                                                     'DB_Object_Symbol']][0] + 1
            if discarded_data[new_keys]['Aspect'] == "F":
                discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][1] = \
                    discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][1] + 1
            elif discarded_data[new_keys]['Aspect'] == "P":
                discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][2] = \
                    discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][2] + 1
            elif discarded_data[new_keys]['Aspect'] == "C":
                discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][3] = \
                    discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][3] + 1
            # print("yes")
        else:
            discardedDict[discarded_data[new_keys]['DB_Object_Symbol']] = [1, 0, 0, 0]
            if discarded_data[new_keys]['Aspect'] == "F":
                discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][1] = \
                    discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][1] + 1
            elif discarded_data[new_keys]['Aspect'] == "P":
                discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][2] = \
                    discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][2] + 1
            elif discarded_data[new_keys]['Aspect'] == "C":
                discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][3] = \
                    discardedDict[discarded_data[new_keys]['DB_Object_Symbol']][3] + 1
        # print(discarded_data[new_keys]['DB_Object_Symbol']+"\t"+str(discardedDict[discarded_data[new_keys]['DB_Object_Symbol']]))

    discFH1 = open("HTPProteins_New.tsv", "w")
    for i in discardedDict:
        discFH1.write(str(i) + "\t" + str(discardedDict[i]) + "\n")
    discFH1.close()

    removedProt = 0
    for keys in originalData:
        # print(keys)

        if keys not in data.keys():
            # print(originalData[keys]['DB_Object_Symbol']+"\t"+originalData[keys]['GO_ID'])
            removedProt += 1
            completelyRemoved.append(originalData[keys]['DB_Object_Symbol'])
        else:
            if data[keys]['DB_Object_Symbol'] in debiasedDict:
                debiasedDict[data[keys]['DB_Object_Symbol']][0] = \
                    debiasedDict[data[keys]['DB_Object_Symbol']][0] + 1
                if data[keys]['Aspect'] == "F":
                    debiasedDict[data[keys]['DB_Object_Symbol']][1] = \
                        debiasedDict[data[keys]['DB_Object_Symbol']][1] + 1
                elif data[keys]['Aspect'] == "P":
                    debiasedDict[data[keys]['DB_Object_Symbol']][2] = \
                        debiasedDict[data[keys]['DB_Object_Symbol']][2] + 1
                elif data[keys]['Aspect'] == "C":
                    debiasedDict[data[keys]['DB_Object_Symbol']][3] = \
                        debiasedDict[data[keys]['DB_Object_Symbol']][3] + 1
            else:
                debiasedDict[data[keys]['DB_Object_Symbol']] = [1, 0, 0, 0]
                if data[keys]['Aspect'] == "F":
                    debiasedDict[data[keys]['DB_Object_Symbol']][1] = \
                        debiasedDict[data[keys]['DB_Object_Symbol']][1] + 1
                elif data[keys]['Aspect'] == "P":
                    debiasedDict[data[keys]['DB_Object_Symbol']][2] = \
                        debiasedDict[data[keys]['DB_Object_Symbol']][2] + 1
                elif data[keys]['Aspect'] == "C":
                    debiasedDict[data[keys]['DB_Object_Symbol']][3] = \
                        debiasedDict[data[keys]['DB_Object_Symbol']][3] + 1

            if data[keys]['DB_Object_Symbol'] in newproteinDict:
                newproteinDict[data[keys]['DB_Object_Symbol']] = newproteinDict[data[keys]['DB_Object_Symbol']] + 1
                # print("yes")
            else:
                newproteinDict[data[keys]['DB_Object_Symbol']] = 1
        if originalData[keys]['DB_Object_Symbol'] in proteinDict:
            proteinDict[originalData[keys]['DB_Object_Symbol']][0] = \
                proteinDict[originalData[keys]['DB_Object_Symbol']][0] + 1
            if originalData[keys]['Aspect'] == "F":
                proteinDict[originalData[keys]['DB_Object_Symbol']][1] = \
                    proteinDict[originalData[keys]['DB_Object_Symbol']][1] + 1
            elif originalData[keys]['Aspect'] == "P":
                proteinDict[originalData[keys]['DB_Object_Symbol']][2] = \
                    proteinDict[originalData[keys]['DB_Object_Symbol']][2] + 1
            elif originalData[keys]['Aspect'] == "C":
                proteinDict[originalData[keys]['DB_Object_Symbol']][3] = \
                    proteinDict[originalData[keys]['DB_Object_Symbol']][3] + 1
            # print("yes")
        else:
            proteinDict[originalData[keys]['DB_Object_Symbol']] = [1, 0, 0, 0]
            if originalData[keys]['Aspect'] == "F":
                proteinDict[originalData[keys]['DB_Object_Symbol']][1] = \
                    proteinDict[originalData[keys]['DB_Object_Symbol']][1] + 1
            elif originalData[keys]['Aspect'] == "P":
                proteinDict[originalData[keys]['DB_Object_Symbol']][2] = \
                    proteinDict[originalData[keys]['DB_Object_Symbol']][2] + 1
            elif originalData[keys]['Aspect'] == "C":
                proteinDict[originalData[keys]['DB_Object_Symbol']][3] = \
                    proteinDict[originalData[keys]['DB_Object_Symbol']][3] + 1
            # print(proteinDict[originalData[keys]['DB_Object_Symbol']])
    partialCounter = 0
    compRemoved = 0
    print("Removed annotations " + str(removedProt))
    fFH1 = open("PartiallyRemoved.txt", "w")
    cFH1 = open("CompletelyRemovedProteinList.txt", "w")
    dFH1 = open("AllProteins_New.tsv", "w")  # previously: ProteinListWithAnnotations_New.tsv
    for i in proteinDict:
        # if newproteinDict[i] < proteinDict[i]:
        dFH1.write(str(i) + "\t" + str(proteinDict[i]) + "\n")  # i is protein common name, proteinDict[i] = [All, F, P, C] annotation count
        if i not in newproteinDict:
            cFH1.write(i + "\n")
            compRemoved += 1
        else:
            diff = proteinDict[i][0] - newproteinDict[i]
            partialCounter += 1
            fFH1.write(i + "\n")
            # print(str(i)+"\t"+str(diff))

    lFH1 = open("LTPProteins_New.tsv", "w")
    for k in debiasedDict:
        lFH1.write(str(k) + "\t" + str(debiasedDict[k]) + "\n")
    lFH1.close()
    fFH1.close()
    cFH1.close()
    dFH1.close()
    # compRemoved = np.unique(completelyRemoved)
    # for c in compRemoved:
    #    cFH1.write(c+"\n")
    # cFH1.close()
    print("Total proteins " + str(len(proteinDict)))
    print("Completely removed " + str(compRemoved))
    print("Partially removed " + str(partialCounter))
    # print(removedProt)
    # filepath = options.histogram
    # outputfilename = "/."
    fileTemp1 = open("prev.txt", "w")
    for i in prev:
        fileTemp1.write(str(i) + "\t" + str(prev[i]) + "\n")
    fileTemp1.close()
    fileTemp2 = open("later.txt", "w")
    for j in lat:
        fileTemp2.write(str(j) + "\t" + str(lat[j]) + "\n")
    fileTemp2.close()
    # print(prev)
    # if (filepath[-1] == '/'):
    #     if os.path.isdir(filepath) == False:
    #         os.makedirs(filepath)
    #     filepath += species.split(".gaf")[0].split("_")[-1] + "_" + "histogram_" + "_".join(
    #         outputfilename.split("/")[-1].split("_")[1:]) + ".png"
    # elif (filepath == "."):
    #     # print("I am here")
    #     filepath = species.split(".gaf")[0].split("_")[-1] + "_" + "histogram_" + "_".join(
    #         outputfilename.split("/")[-1].split("_")[1:]) + ".png"
    #     # print(filepath)
    # else:
    #     if "/" in filepath and os.path.isdir('/'.join(filepath.split("/")[:-1])) == False:
    #         os.makedirs('/'.join(filepath.split("/")[:-1]))
    #     temp = species.split(".gaf")[0].split("_")[-1] + "_"
    #     if ("." not in filepath and "/" not in filepath):
    #         filepath = species.split(".gaf")[0].split("_")[-1] + "_" + "histogram_" + "_".join(
    #             outputfilename.split("/")[-1].split("_")[1:]) + ".png"
    #     elif ("." in filepath.split("/")[-1]):
    #         filepath = "/".join((".".join(filepath.split(".")[:-1]) + ".png").split("/")[:-1]) + "/" + temp + \
    #                    (".".join(filepath.split(".")[:-1]) + ".png").split("/")[-1]
    #     else:
    #         filepath = "/".join(filepath.split("/")[:-1]) + "/" + temp + filepath.split("/")[-1] + ".png"
    # print("IMG filepath",filepath)

    prev_P = freqGO_TERM(data)
    prev_val_p = [prev_P[key] for key in prev_P]

    prev_val = [prev[key] for key in prev]

    # print(sum(prev_val),np.mean(prev_val))
    new_prev_val = [-math.log(val / sum(prev_val), 2) for val in prev_val]
    # print(sum(prev_val))
    # prev_val=new_prev_val
    lat_val = [lat[key] for key in lat]
    # print(sum(lat_val),np.mean(lat_val))
    new_lat_val = [-math.log(val / sum(lat_val), 2) for val in lat_val]
    # print(sum(lat_val))
    # lat_val=new_lat_val

    ##############################################################################################################

    #############################################################################################################

    go_freq, go_count = freqGO_TERM_Protein(discarded_data)

    # Calculate number of C, P, F
    num_F = 0
    num_P = 0
    num_C = 0
    aspect_dict = {}
    for keys in originalData:
        if originalData[keys]['Aspect'] == "F":
            num_F += 1
            aspect_dict[originalData[keys]['GO_ID']] = originalData[keys]['Aspect']
        elif originalData[keys]['Aspect'] == "P":
            num_P += 1
            aspect_dict[originalData[keys]['GO_ID']] = originalData[keys]['Aspect']
        elif originalData[keys]['Aspect'] == "C":
            num_C += 1
            aspect_dict[originalData[keys]['GO_ID']] = originalData[keys]['Aspect']
    print("F", num_F)
    print("P", num_P)
    print("C", num_C)

    print(len(prev))
    print("Aspect dict", len(aspect_dict))

    ############################################################################################################

    prev_P = freqGO_TERM(originalData)
    prev_val_p = [prev_P[key] for key in prev_P]

    prev_val = [prev[key] for key in prev]

    # print(sum(prev_val),np.mean(prev_val))
    new_prev_val = [-math.log(val / sum(prev_val), 2) for val in prev_val]
    # print(sum(prev_val))
    # prev_val=new_prev_val
    lat_val = [lat[key] for key in lat]
    # print(sum(lat_val),np.mean(lat_val))
    new_lat_val = [-math.log(val / sum(lat_val), 2) for val in lat_val]

    go_to_ic_dict_f = {}
    go_to_ic_dict_p = {}
    go_to_ic_dict_c = {}
    ft_m = open("all_deltaIC.txt", "w")

    # ia = pd.read_csv("WC_IC_MFO.txt",header=None,sep="\t",names=["GO","nom","denom","prob","WCIC"])
    # print("Sample IA per term dict:\n", ia.head(5))
    # go_to_ia_dict_f = dict(zip(ia.GO, ia.WCIC))
    # ia = pd.read_csv("WC_IC_BPO.txt",header=None,sep="\t",names=["GO","nom","denom","prob","WCIC"])
    # go_to_ia_dict_p = dict(zip(ia.GO, ia.WCIC))
    # ia = pd.read_csv("WC_IC_CCO.txt",header=None,sep="\t",names=["GO","nom","denom","prob","WCIC"])
    # go_to_ia_dict_c = dict(zip(ia.GO, ia.WCIC))


    for i in prev_P:
        if aspect_dict[i] == "F":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(prev_P[i] / num_F) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_F, 2)) + "\t" + str(
                (prev_P[i] / num_F) * (-math.log(prev_P[i] / num_F, 2))) + "\n")  # Shannon entropy -P*log2(P)
            go_to_ic_dict_f[i] = -math.log(prev_P[i] / num_F, 2)
            # go_to_ia_dict_f[i] = go_to_ia_dict[i]
        if aspect_dict[i] == "P":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(
                prev_P[i] / num_P) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_P, 2)) + "\t" + str(
                (prev_P[i] / num_P) * (-math.log(prev_P[i] / num_P, 2))) + "\n")
            go_to_ic_dict_p[i] = -math.log(prev_P[i] / num_P, 2)
            # go_to_ia_dict_p[i] = go_to_ia_dict[i]
        if aspect_dict[i] == "C":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(
                prev_P[i] / num_C) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_C, 2)) + "\t" + str(
                (prev_P[i] / num_C) * (-math.log(prev_P[i] / num_C, 2))) + "\n")
            go_to_ic_dict_c[i] = -math.log(prev_P[i] / num_C, 2)
            # go_to_ia_dict_c[i] = go_to_ia_dict[i]

        # print(i, prev_P[i], sum(prev_val), prev_P[i]/sum(prev_val))
        # else:
        #     ft_m.write(str(lat[i]) + "\t" + str(-math.log(prev[i] / sum(prev_val), 2)) + "\t" + str(
        #         -math.log(lat[i] / sum(prev_val), 2)) + "\n")

    # for k in range(0,len(tempArrayPrev)):
    # ft1.write(str(tempArrayPrev[k]-tempArrayLat[k])+"\n")
    ft_m.close()

    afhi = open("AllProteins_IC.tsv", "w")
    # fhi = open("All_Proteins_With_IC.tsv", "w")   # with IC
    # fhi.write("Protein\tAnnotations\tAnnotation_Set\tMFO\tMFO_IC\tBPO\tBPO_IC\tCCO\tCCO_IC" + "\n") # with IC

    prot_go_dict = {}
    ic_dict = {}
    for keys in originalData:
        if originalData[keys]['DB_Object_Symbol'] in prot_go_dict:
            prot_go_dict[originalData[keys]['DB_Object_Symbol']].append(originalData[keys]['GO_ID'])
        else:
            prot_go_dict[originalData[keys]['DB_Object_Symbol']] = [originalData[keys]['GO_ID']]

    # Print Phillip Lord IC to file "AllProteins_IC.tsv"
    for i in prot_go_dict:
        # fhi.write(str(i) + "\t" + str(prot_go_dict[i]) + "\t" + str(set(prot_go_dict[i])) + "\t")   # with IC
        afhi.write(str(i) + "\t")
        prot_go_dict[i] = set(prot_go_dict[i]) # Remove for with IC stuff
        print(i, prot_go_dict[i])
        mfo_ic = 0
        bpo_ic = 0
        cco_ic = 0

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_f.keys():
                # print(go_to_ic_dict_f[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_f[j]
        mfo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t") # with IC
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_p.keys():
                # print(go_to_ic_dict_p[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_p[j]
        bpo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t") # with IC
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_c.keys():
                # print(go_to_ic_dict_c[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_c[j]
        cco_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\n") # with IC
        afhi.write(
            "[" + str(mfo_ic + bpo_ic + cco_ic) + ", " + str(mfo_ic) + ", " + str(bpo_ic) + ", " + str(cco_ic) + "]\n")
        ic_dict[i] =[(mfo_ic + bpo_ic + cco_ic),(mfo_ic), (bpo_ic), (cco_ic)]

    # fhi.close()
    afhi.close()

    # Print Wyatt Clark information accretion to file "AllProteins_IA.tsv, tab-separated not dict format
    bfhi = open("AllProteins_IA.tsv", "w")
    ia_dict = {}

    for i in prot_go_dict:
        # fhi.write(str(i) + "\t" + str(prot_go_dict[i]) + "\t" + str(set(prot_go_dict[i])) + "\t")   # with IC
        bfhi.write(str(i) + "\t")
        prot_go_dict[i] = set(prot_go_dict[i]) # Remove for with IC stuff
        print(i, prot_go_dict[i])
        mfo_ic = 0
        bpo_ic = 0
        cco_ic = 0

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_f.keys():
                # print(go_to_ia_dict_f[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_f[j]
        mfo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t") # with IC
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_p.keys():
                # print(go_to_ia_dict_p[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_p[j]
        bpo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t") # with IC
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_c.keys():
                # print(go_to_ia_dict_c[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_c[j]
        cco_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\n") # with IC
        bfhi.write(
            "[" + str(mfo_ic + bpo_ic + cco_ic) + ", " + str(mfo_ic) + ", " + str(bpo_ic) + ", " + str(cco_ic) + "]\n")
        ia_dict[i] =[(mfo_ic + bpo_ic + cco_ic),(mfo_ic), (bpo_ic), (cco_ic)]
    # fhi.close()
    bfhi.close()

    ##############################################################################################################

    # ft1 = open("deltaIC.txt", "w")
    # sigmaDelta = 0
    # # tempArrayPrev = []
    # # tempArrayLat = []
    # fileTemp1 = open("prevIC.txt", "w")
    # for i in lat:
    #     ft1.write(str(i) + "\t" + str(prev[i]) + "\t" + str(lat[i]) + "\t" + str(
    #         -math.log(prev[i] / sum(prev_val), 2)) + "\t" + str(-math.log(lat[i] / sum(lat_val), 2)) + "\t" + str(
    #         -math.log(prev[i] / sum(prev_val), 2) - (-math.log(lat[i] / sum(lat_val), 2))) + "\t" + str(
    #         (-math.log(prev[i] / sum(prev_val), 2) - (-math.log(lat[i] / sum(lat_val), 2))) ** 2) + "\n")
    #     fileTemp1.write(str(i) + "\t" + str(prev[i]) + "\n")
    #     sigmaDelta += (-math.log(prev[i] / sum(prev_val), 2) - (-math.log(lat[i] / sum(lat_val), 2))) ** 2
    #     # tempArrayPrev.append(i)
    # fileTemp2 = open("laterIC.txt", "w")
    # print(sigmaDelta)
    # for j in new_lat_val:
    #     fileTemp2.write(str(j) + "\n")
    # tempArrayLat.append(j)

    # ft_m = open("modified_deltaIC.txt", "w")
    # for i in prev:
    #     ft_m.write(str(i) + "\t" + str(prev[i]) + "\t")
    #     if i not in lat.keys():
    #         ft_m.write(str(0) + "\t" + str(-math.log(prev[i] / sum(prev_val), 2)) + "\t" + " " + "\n")
    #     else:
    #         ft_m.write(str(lat[i]) + "\t" + str(-math.log(prev[i] / sum(prev_val), 2)) + "\t" + str(
    #             -math.log(lat[i] / sum(prev_val), 2)) + "\n")

    ''' commenting to line 1309 for corrected calculation of IC between ltp and htp 28 Jun 2022

    prev_P = freqGO_TERM(data)
    prev_val_p = [prev_P[key] for key in prev_P]

    prev_val = [prev[key] for key in prev]

    # print(sum(prev_val),np.mean(prev_val))
    new_prev_val = [-math.log(val / sum(prev_val), 2) for val in prev_val]
    # print(sum(prev_val))
    # prev_val=new_prev_val
    lat_val = [lat[key] for key in lat]
    # print(sum(lat_val),np.mean(lat_val))
    new_lat_val = [-math.log(val / sum(lat_val), 2) for val in lat_val]

    go_to_ic_dict_f = {}
    go_to_ic_dict_p = {}
    go_to_ic_dict_c = {}
    '''

    ft_m = open("ltp_deltaIC.txt", "w")
    for i in prev_P:
        if aspect_dict[i] == "F":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(prev_P[i] / num_F) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_F, 2)) + "\t" + str(
                (prev_P[i] / num_F) * (-math.log(prev_P[i] / num_F, 2))) + "\n")
            # go_to_ic_dict_f[i] = -math.log(prev_P[i] / num_F, 2) # 28 Jun 22 fix
        if aspect_dict[i] == "P":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(
                prev_P[i] / num_P) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_P, 2)) + "\t" + str(
                (prev_P[i] / num_P) * (-math.log(prev_P[i] / num_P, 2))) + "\n")
            # go_to_ic_dict_p[i] = -math.log(prev_P[i] / num_P, 2)  # 28 Jun 22 fix
        if aspect_dict[i] == "C":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(
                prev_P[i] / num_C) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_C, 2)) + "\t" + str(
                (prev_P[i] / num_C) * (-math.log(prev_P[i] / num_C, 2))) + "\n")
            # go_to_ic_dict_c[i] = -math.log(prev_P[i] / num_C, 2)  # 28 Jun 22 fix

        # print(i, prev_P[i], sum(prev_val), prev_P[i]/sum(prev_val))
        # else:
        #     ft_m.write(str(lat[i]) + "\t" + str(-math.log(prev[i] / sum(prev_val), 2)) + "\t" + str(
        #         -math.log(lat[i] / sum(prev_val), 2)) + "\n")

    # for k in range(0,len(tempArrayPrev)):
    # ft1.write(str(tempArrayPrev[k]-tempArrayLat[k])+"\n")
    ft_m.close()


    # fhi = open("LTP_Proteins_With_IC.tsv", "w")
    afhi = open("LTPProteins_IC.tsv", "w")
    # fhi.write("Protein\tAnnotations\tAnnotation_Set\tMFO\tMFO_IC\tBPO\tBPO_IC\tCCO\tCCO_IC" + "\n")   # For with IC stuff

    prot_go_dict = {}
    for keys in data:
        if data[keys]['DB_Object_Symbol'] in prot_go_dict:
            prot_go_dict[data[keys]['DB_Object_Symbol']].append(data[keys]['GO_ID'])
        else:
            prot_go_dict[data[keys]['DB_Object_Symbol']] = [data[keys]['GO_ID']]

    for i in prot_go_dict:
        # fhi.write(str(i) + "\t" + str(prot_go_dict[i]) + "\t" + str(set(prot_go_dict[i])) + "\t") # For with IC stuff
        afhi.write(str(i) + "\t")
        prot_go_dict[i] = set(prot_go_dict[i])    # Remove for with IC
        # print(i, prot_go_dict[i])
        mfo_ic = 0
        bpo_ic = 0
        cco_ic = 0

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_f.keys():
                # print(go_to_ic_dict_f[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_f[j]
        mfo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For with IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_p.keys():
                # print(go_to_ic_dict_p[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_p[j]
        bpo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For with IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_c.keys():
                # print(go_to_ic_dict_c[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_c[j]
        cco_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\n")   # For with IC stuff
        afhi.write(
           "[" + str(mfo_ic + bpo_ic + cco_ic) + ", " + str(mfo_ic) + ", " + str(bpo_ic) + ", " + str(cco_ic) + "]\n")

    # fhi.close()
    afhi.close()


    # fhi = open("LTP_Proteins_With_IA.tsv", "w")
    bfhi = open("LTPProteins_IA.tsv", "w")
    # fhi.write("Protein\tAnnotations\tAnnotation_Set\tMFO\tMFO_IC\tBPO\tBPO_IC\tCCO\tCCO_IC" + "\n")   # For with IC stuff

    prot_go_dict = {}
    for keys in data:
        if data[keys]['DB_Object_Symbol'] in prot_go_dict:
            prot_go_dict[data[keys]['DB_Object_Symbol']].append(data[keys]['GO_ID'])
        else:
            prot_go_dict[data[keys]['DB_Object_Symbol']] = [data[keys]['GO_ID']]

    for i in prot_go_dict:
        # fhi.write(str(i) + "\t" + str(prot_go_dict[i]) + "\t" + str(set(prot_go_dict[i])) + "\t") # For with IC stuff
        bfhi.write(str(i) + "\t")
        prot_go_dict[i] = set(prot_go_dict[i])  # Remove for with IC
        # print(i, prot_go_dict[i])
        mfo_ic = 0
        bpo_ic = 0
        cco_ic = 0

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_f.keys():
                # print(go_to_ia_dict_f[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_f[j]
        mfo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For with IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_p.keys():
                # print(go_to_ic_dict_p[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_p[j]
        bpo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For with IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # For without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_c.keys():
                # print(go_to_ic_dict_c[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_c[j]
        cco_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\n")   # For with IC stuff
        bfhi.write(
            "[" + str(mfo_ic + bpo_ic + cco_ic) + ", " + str(mfo_ic) + ", " + str(bpo_ic) + ", " + str(cco_ic) + "]\n")

    # fhi.close()
    bfhi.close()

    #########################################################################################################


    ##########################################################################################################

    ''' commenting to line 1389 for corrected calculation of IC between ltp and htp 28 Jun 2022
    prev_P = freqGO_TERM(discarded_data)
    prev_val_p = [prev_P[key] for key in prev_P]

    prev_val = [prev[key] for key in prev]

    # print(sum(prev_val),np.mean(prev_val))
    new_prev_val = [-math.log(val / sum(prev_val), 2) for val in prev_val]
    # print(sum(prev_val))
    # prev_val=new_prev_val
    lat_val = [lat[key] for key in lat]
    # print(sum(lat_val),np.mean(lat_val))
    new_lat_val = [-math.log(val / sum(lat_val), 2) for val in lat_val]

    go_to_ic_dict_f = {}
    go_to_ic_dict_p = {}
    go_to_ic_dict_c = {}
    '''

    ft_m = open("htp_deltaIC.txt", "w")
    for i in prev_P:
        if aspect_dict[i] == "F":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(prev_P[i] / num_F) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_F, 2)) + "\t" + str(
                (prev_P[i] / num_F) * (-math.log(prev_P[i] / num_F, 2))) + "\n")
            # go_to_ic_dict_f[i] = -math.log(prev_P[i] / num_F, 2) # 28 Jun 22 fix
        if aspect_dict[i] == "P":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(
                prev_P[i] / num_P) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_P, 2)) + "\t" + str(
                (prev_P[i] / num_P) * (-math.log(prev_P[i] / num_P, 2))) + "\n")
            # go_to_ic_dict_p[i] = -math.log(prev_P[i] / num_P, 2) # 28 Jun 22 fix
        if aspect_dict[i] == "C":
            ft_m.write(str(i) + "\t" + str(prev_P[i]) + "\t" + aspect_dict[i] + "\t" + str(
                prev_P[i] / num_C) + "\t")
            ft_m.write(str(-math.log(prev_P[i] / num_C, 2)) + "\t" + str(
                (prev_P[i] / num_C) * (-math.log(prev_P[i] / num_C, 2))) + "\n")
            # go_to_ic_dict_c[i] = -math.log(prev_P[i] / num_C, 2) # 28 Jun 22 fix

        # print(i, prev_P[i], sum(prev_val), prev_P[i]/sum(prev_val))
        # else:
        #     ft_m.write(str(lat[i]) + "\t" + str(-math.log(prev[i] / sum(prev_val), 2)) + "\t" + str(
        #         -math.log(lat[i] / sum(prev_val), 2)) + "\n")

    # for k in range(0,len(tempArrayPrev)):
    # ft1.write(str(tempArrayPrev[k]-tempArrayLat[k])+"\n")
    ft_m.close()
    
    ##########################################################################################################

    # ft1.close()
    fileTemp1.close()
    fileTemp2.close()

    # fhi = open("HTP_Proteins_With_IC.tsv", "w")
    afhi = open("HTPProteins_IC.tsv", "w")
    # fhi.write("Protein\tAnnotations\tAnnotation_Set\tMFO\tMFO_IC\tBPO\tBPO_IC\tCCO\tCCO_IC" + "\n")   # With IC stuff

    prot_go_dict = {}
    for keys in discarded_data:
        if discarded_data[keys]['DB_Object_Symbol'] in prot_go_dict:
            prot_go_dict[discarded_data[keys]['DB_Object_Symbol']].append(discarded_data[keys]['GO_ID'])
        else:
            prot_go_dict[discarded_data[keys]['DB_Object_Symbol']] = [discarded_data[keys]['GO_ID']]

    for i in prot_go_dict:
        # fhi.write(str(i) + "\t" + str(prot_go_dict[i]) + "\t" + str(set(prot_go_dict[i])) + "\t") # With IC stuff
        afhi.write(str(i) + "\t")
        prot_go_dict[i] = set(prot_go_dict[i])    # removed for with IC
        print(i, prot_go_dict[i])
        mfo_ic = 0
        bpo_ic = 0
        cco_ic = 0

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_f.keys():
                # print(go_to_ic_dict_f[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_f[j]
        mfo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # With IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # Without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_p.keys():
                # print(go_to_ic_dict_p[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_p[j]
        bpo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # With IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # Without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ic_dict_c.keys():
                # print(go_to_ic_dict_c[j])
                all_j.append(j)
                all_j_ic += go_to_ic_dict_c[j]
        cco_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\n")   # With IC stuff
        afhi.write(
            "[" + str(mfo_ic + bpo_ic + cco_ic) + ", " + str(mfo_ic) + ", " + str(bpo_ic) + ", " + str(cco_ic) + "]\n")


    # fhi.close()
    afhi.close()


    # fhi = open("HTP_Proteins_With_IC.tsv", "w")
    bfhi = open("HTPProteins_IA.tsv", "w")
    # fhi.write("Protein\tAnnotations\tAnnotation_Set\tMFO\tMFO_IC\tBPO\tBPO_IC\tCCO\tCCO_IC" + "\n")   # With IC stuff

    prot_go_dict = {}
    for keys in discarded_data:
        if discarded_data[keys]['DB_Object_Symbol'] in prot_go_dict:
            prot_go_dict[discarded_data[keys]['DB_Object_Symbol']].append(discarded_data[keys]['GO_ID'])
        else:
            prot_go_dict[discarded_data[keys]['DB_Object_Symbol']] = [discarded_data[keys]['GO_ID']]

    for i in prot_go_dict:
        # fhi.write(str(i) + "\t" + str(prot_go_dict[i]) + "\t" + str(set(prot_go_dict[i])) + "\t") # With IC stuff
        bfhi.write(str(i) + "\t")
        prot_go_dict[i] = set(prot_go_dict[i])  # removed for with IC
        print(i, prot_go_dict[i])
        mfo_ic = 0
        bpo_ic = 0
        cco_ic = 0

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_f.keys():
                # print(go_to_ic_dict_f[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_f[j]
        mfo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # With IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # Without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_p.keys():
                # print(go_to_ic_dict_p[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_p[j]
        bpo_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # With IC stuff
        ## fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\t")   # Without IC stuff

        all_j = []
        all_j_ic = 0
        for j in prot_go_dict[i]:
            if j in go_to_ia_dict_c.keys():
                # print(go_to_ic_dict_c[j])
                all_j.append(j)
                all_j_ic += go_to_ia_dict_c[j]
        cco_ic = all_j_ic
        # fhi.write(str(all_j) + "\t" + str(all_j_ic) + "\n")   # With IC stuff
        bfhi.write(
            "[" + str(mfo_ic + bpo_ic + cco_ic) + ", " + str(mfo_ic) + ", " + str(bpo_ic) + ", " + str(cco_ic) + "]\n")

    # fhi.close()
    bfhi.close()

    """prev_val=[]
    for key in prev:
        prev_val.append(prev[key])"""


    # binsize = 100
    # plt.hist(new_prev_val, bins=binsize, color='r', label="Before Debiasing")
    # plt.hist(new_lat_val, bins=binsize, color='b', label="After Debiasing")
    # # plt.xlim((0,max(max(new_prev_val)+0.1,max(new_lat_val)+0.1)))
    # plt.xlabel("Information Content")
    # plt.ylabel("Frequency")
    # plt.title("Histogram for " + species.split(".gaf")[0].split("_")[-1] + "\n" + msg)
    #
    # if options.histogram == "1":
    #     plt.legend(bbox_to_anchor=(0.4, 1))
    #     # plt.show()
    # else:
    #     # plt.figure(figsize=(70,70))
    #     plt.legend(bbox_to_anchor=(0.45, 1))
    #     plt.savefig(filepath, dpi=900)
    #     fhw = open(filepath[:-3] + ".txt", "w")
    #     fhw.write(str(new_prev_val) + "\n" + str(new_lat_val))
    # plt.close()
    return proteinDict, ic_dict, ia_dict

def freqGO_TERM(data):
    go_to_freq = dict()
    for annotation in data:
        if data[annotation]['GO_ID'] in go_to_freq:
            go_to_freq[data[annotation]['GO_ID']] += 1
        else:
            go_to_freq[data[annotation]['GO_ID']] = 1
    return go_to_freq


def freqGO_TERM_Protein(data):
    go_to_freq = dict()
    go_count = dict()
    for annotation in data:
        # print(annotation)
        if data[annotation]['DB_Object_Symbol'] in go_to_freq:
            go_to_freq[data[annotation]['DB_Object_Symbol']].append(data[annotation]['GO_ID'])
            go_count[data[annotation]['DB_Object_Symbol']] += 1
        else:
            # print(data[annotation]['DB_Object_Symbol'], data[annotation]['GO_ID'])
            go_to_freq[data[annotation]['DB_Object_Symbol']] = [data[annotation]['GO_ID']]
            go_count[data[annotation]['DB_Object_Symbol']] = 1
    return go_to_freq, go_count


def writeContentsToFile(dataDict, fileName):
    F = open(fileName + ".txt", "w")
    for key in dataDict:
        F.write(str(key) + "\t" + str(dataDict[key]) + "\n")
    F.close()


def main():
    global verbose, options, report
    commandLineArg = sys.argv
    if len(commandLineArg) == 1:
        print("Please use the --help option to get usage information")
    # Parse command line arguments
    options = parseCommandLineArguments()
    if options.recalculate != 0 and (
            options.info_threshold_Wyatt_Clark == None and options.info_threshold_Wyatt_Clark_percentile == None):
        print("Error in arguments. You must provide Wyatt Clark in order to recalculate")
        # print(options.recalculate,options.info_threshold_Wyatt_Clark,options.info_threshold_Wyatt_Clark_percentile )
        exit()
    if (options.aspect == None):
        options.aspect = ["F", "P", "C"]
    if options.verbose != 0:
        verbose = 1

    # vprint( options )

    options.output = changeNameofOutputFiles(options)
    if options.report != None:
        report = 1
        options.report = createReportFile(options.report, options.output[0])
        reportdict = dict()
        for input in options.input:
            reportdict[input.split(".gaf")[0].split("_")[1]] = []
    # exit()
    data = dict()
    originalData = dict()
    # data = convertToDictionary( options.input )
    for file_num, eachinputfile in enumerate(options.input):
        species = eachinputfile.split(".gaf")[0].split("_")[1]
        gafoutput = GOA._gaf20iterator(open(eachinputfile, "r"))
        data = convertFromGAFToRequiredFormat(gafoutput)

        # data.update(temp_data)
        # vprint("Processing file ",eachInputFile,file_num)
        # if options.verbose != 0:
        #     printDetailsAboutData(data)
        # print(data)
        originalData = data
        # print(originalData)
        vprint()
        if (options.select_references != None or options.select_references_inverse != None):
            vprint("Number of annotations before choosing proteins based on references ", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                report_row.append("References")
                report_row.append(countProteins(data))
            data = chooseProteinsBasedOnReferences(data, options.select_references, options.select_references_inverse)
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[eachinputfile].append(report_row)
            vprint("Number of annotations before choosing proteins based on references ", len(data))
            vprint("Data discarded ", prev_len - len(data))

            vprint()
        if (options.evidence != None or options.evidence_inverse != None):
            vprint("Number of annotations before choosing proteins based on Evidence Codes ", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                report_row.append("Evidence")
                report_row.append(countProteins(data))
            data = chooseProteinsBasedOnEvidenceCodes(data, evidence_list=options.evidence,
                                                      evidence_inverse_list=options.evidence_inverse)
            vprint("Number of annotations after choosing proteins based on Evidence Codes ", len(data))
            vprint("Data discarded ", prev_len - len(data))
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[eachinputfile].append(report_row)
            vprint()
        if (options.aspect != None):
            vprint("Number of annotations before choosing proteins based on aspect ", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                report_row.append("Aspect")
                report_row.append(countProteins(data))
            data = chooseGOBasedOnAspect(data, aspect=options.aspect)
            vprint("Number of annotations after choosing proteins based on aspect ", len(data))
            vprint("Data discarded ", prev_len - len(data))
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[species].append(report_row)
            vprint()
        if (options.assigned_by != None or options.assigned_by_inverse != None):
            vprint("Number of annotations before choosing proteins based on assigned by", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                report_row.append("Assigned By")
                report_row.append(countProteins(data))
            data = chooseGOBasedOnAssignedBy(data, options.assigned_by, options.assigned_by_inverse)
            vprint("Number of annotations after choosing proteins based on assigned_by ", len(data))
            vprint("Data discarded ", prev_len - len(data))
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[eachinputfile].append(report_row)
            vprint()
        if (options.date_before != None or options.date_after != None):
            vprint("Number of annotations before choosing proteins based on provided range of dates ", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                report_row.append("Date")
                report_row.append(countProteins(data))
            data = chooseAnnotationsBasedOnDate(data, options.date_before, options.date_after)
            vprint("Number of annotations after choosing proteins based on provided range of dates ", len(data))
            vprint("Data discarded ", prev_len - len(data))
        #     if (report == 1):
        #         report_row.append(countProteins(data))
        #         report_row.append(prev_len)
        #         report_row.append(len(data))
        #         reportdict[eachinputfile].append(report_row)
        #     vprint()

        if (options.info_threshold_Phillip_Lord != None or options.info_threshold_Phillip_Lord_percentile != None):
            vprint("Number of annotations before choosing proteins based on Phillip Lord Threshold ", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                report_row.append("Phillip Lord Threshold")
                report_row.append(countProteins(data))
            data = calculatePhillipLordInformationContent(data, options.info_threshold_Phillip_Lord,
                                                          options.info_threshold_Phillip_Lord_percentile)
            vprint("Number of annotations after choosing proteins based on Phillip Lord Threshold ", len(data))
            vprint("Data discarded ", prev_len - len(data))
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[eachinputfile].append(report_row)
            vprint()
        if (options.info_threshold_Wyatt_Clark or options.info_threshold_Wyatt_Clark_percentile != None):
            vprint("Number of annotations before choosing proteins based on Wyatt Clark Threshold ", len(data))
            prev_len = len(data)
            # for attnid in data:
            #     print(data[attnid])
            if (report == 1):
                report_row = []
                report_row.append("Wyatt Clark Threshold")
                report_row.append(countProteins(data))
            data, goia_f, goia_p, goia_c = calculateWyattClarkInformationContent(data, int(options.recalculate),
                                                         options.info_threshold_Wyatt_Clark,
                                                         options.info_threshold_Wyatt_Clark_percentile,
                                                         options.output, file_num)
            vprint("Number of annotations after choosing proteins based on Wyatt Clark Threshold ", len(data))
            vprint("Data discarded ", prev_len - len(data))
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[eachinputfile].append(report_row)
            vprint()

        if (options.cutoff_prot != None or options.cutoff_attn != None):  # move down after calling calculateWC so that 3 go dicts are present for generateHistogram
            vprint("Number of annotations before choosing proteins based on Publications ", len(data))
            prev_len = len(data)
            if (report == 1):
                report_row = []
                if (options.cutoff_prot != None):
                    report_row.append("Protein Cut off")
                else:
                    report_row.append("Annotation Cut off")
                report_row.append(countProteins(data))
            if options.histogram != None:
                prev_go_term_freq = freqGO_TERM(data)
                writeContentsToFile(prev_go_term_freq, "PrevInterimFile")
            data, discarded_data = chooseProteinsBasedOnPublications(data, options.cutoff_prot, options.cutoff_attn)
            if options.histogram != None:
                later_go_term_freq = freqGO_TERM(data)
                writeContentsToFile(later_go_term_freq, "laterInterimFile")
                generateHistogram(data, originalData, discarded_data, eachinputfile, prev_go_term_freq,
                                  later_go_term_freq, goia_f, goia_p, goia_c)

            vprint("Number of annotations after choosing proteins based on Publications ", len(data))
            vprint("Data discarded ", prev_len - len(data))
            if (report == 1):
                report_row.append(countProteins(data))
                report_row.append(prev_len)
                report_row.append(len(data))
                reportdict[species].append(report_row)
            vprint()

        writeToFile(data, options.output[file_num], options.input[file_num])
    if len(options.input) > 1:
        combineOutputFiles(options.output, options)
        if (options.single_file != 0):
            deleteTemporaryFiles(options)
    if report == 1:
        writeReport(options.report, reportdict)


if __name__ == "__main__":
    main()

