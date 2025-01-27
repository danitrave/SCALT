#!/usr/bin/python

""" These are the libraries required  """

import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
import concurrent.futures
import itertools
import json
from operator import itemgetter
start_time = datetime.now()

os.system("mkdir ./INTER_benchmark")

def set_train_sample(S):
    if S == "SM2_pbmc1":
        train_sample = "SM2_pbmc1_pbmc1_SM2"
    elif S == "CL2_pbmc1":
        train_sample = "CL2_pbmc1_pbmc1_Celseq2"
    elif S == "10Xv2A_pbmc1":
        train_sample = "10Xv2A_pbmc1_pbmc1_10x_v2_A"
    elif S == "10Xv2B_pbmc1":
        train_sample = "10Xv2B_pbmc1_pbmc1_10x_v2_B"
    elif S == "DR_pbmc1":
        train_sample = "DR_pbmc1_pbmc1_Drop"
    elif S == "SW_pbmc1":
        train_sample = "SW_pbmc1_pbmc1_Seqwell"
    elif S == "10Xv3_pbmc1":
        train_sample = "10Xv3_pbmc1_pbmc1_10x_v3"
    elif S == "iD_pbmc1":
        train_sample = "iD_pbmc1_pbmc1_inDrops"
    return train_sample

ls = {"10Xv2A_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_10X_V2"],"CL2_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_Celseq2"],"SM2_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_SM2"],
        "10Xv2B_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_10X_V2"],"DR_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_Drop"],
        "SW_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_Seqwell"],"10Xv3_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell"],
        "iD_pbmc1":["pbmc1_SM2","pbmc1_10x_v2_A","pbmc1_10x_v2_B","pbmc1_10x_v3","pbmc1_Celseq2","pbmc1_Drop","pbmc1_inDrops","pbmc1_Seqwell","pbmc2_inDrops"]}

#ls = {"10Xv2A_pbmc1":["pbmc2_10X_V2"],"CL2_pbmc1":["pbmc2_Celseq2"],"SM2_pbmc1":["pbmc2_SM2"],"10Xv2B_pbmc1":["pbmc2_10X_V2"],"DR_pbmc1":["pbmc2_Drop"],
        #"SW_pbmc1":["pbmc2_Seqwell"],"iD_pbmc1":["pbmc2_inDrops"]}
for s in ls:
    os.system("mkdir INTER_benchmark/"+s)
    f = open("processedData/inter/"+s+"/read_counts.tsv","r")
    d = {k:[] for k in ls[s]}
    LINES = f.readlines()
    l = LINES[0].strip("\n").split("\t")[1:]    
    f.close()
    indexes = [x for x in range(len(l))]
    z = list(zip(l,indexes))
    for ini in d:
        for i in z:
            if i[0].startswith(ini):
                d[ini]+=[i,]
            else:
                continue
    df = pd.read_csv("processedData/inter/"+s+"/read_counts.tsv",sep="\t",header=0,index_col=0)
    labs = pd.read_csv("processedData/inter/"+s+"/labels.tsv",sep="\t",header=0,index_col=False)
    train_sample = set_train_sample(s)
    test_samples = []
    for part in d:
        os.system("mkdir INTER_benchmark/"+s+"/"+s+"_"+part)
        test_samples += [s+"_"+part]
        counts_index = [t[1] for t in d[part]]
        ddf = df.iloc[:,counts_index]
        llabs = labs.iloc[counts_index,:]
        ddf.to_csv("INTER_benchmark/"+s+"/"+s+"_"+part+"/read_counts.tsv",sep="\t",header=True,index=True)
        llabs.to_csv("INTER_benchmark/"+s+"/"+s+"_"+part+"/labels.tsv",sep="\t",header=True,index=False)
    
    #### Run benchmark ####
    os.system("mkdir INTER_benchmark/"+s+"/results")
    os.system("cp INTER_benchmark/"+s+"/"+train_sample+"/read_counts.tsv ./train_counts.tsv")
    os.system("cp INTER_benchmark/"+s+"/"+train_sample+"/labels.tsv ./train_labels.tsv")
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    for TEST in test_samples:
        os.system("mkdir INTER_benchmark/"+s+"/results/"+TEST+"_testing")
        os.system("cp INTER_benchmark/"+s+"/"+TEST+"/read_counts.tsv ./benchmark_counts.tsv")
        os.system("cp INTER_benchmark/"+s+"/"+TEST+"/labels.tsv ./benchmark_labels.tsv")
        os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
        os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
        os.system("mv results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTER_benchmark/"+s+"/results/"+TEST+"_testing")

    os.system("mv AnnolistsBuilder_results/ custom/ INTER_benchmark/"+s+"/results")

end_time = datetime.now()
print('Duration inter-experimental set up benchmark: {}'.format(end_time - start_time))
