#!/usr/bin/python

""" Libraries required  """

import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
import time
import concurrent.futures
import json
import warnings
import statistics
from operator import itemgetter
import itertools
import argparse

start_time = datetime.now()

def intra():
    os.system("mkdir ./intra")
    #### liu ####
    os.system("mkdir intra/liu")
    os.system("cp FuEtAll_benchmark/Datasets/Liu/Liu_dataset_intra_validation.R intra/liu")
    df = pd.read_csv("FuEtAll_benchmark/Datasets/Liu/Liu_datasets.csv",header=0,index_col=0)
    labs = pd.read_csv("FuEtAll_benchmark/Datasets/Liu/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("intra/liu/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("intra/liu/Labels.csv",sep=",",header=True,index=False)

    #### Zheng 68K ####
    os.system("mkdir intra/zheng68k")
    df = pd.read_csv("FuEtAll_benchmark/Datasets/Zheng68K/Zheng68K.csv",header=0,index_col=0)
    labs = pd.read_csv("FuEtAll_benchmark/Datasets/Zheng68K/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("intra/zheng68k/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("intra/zheng68k/Labels.csv",sep=",",header=True,index=False)

    #### Zheng sorted ####
    os.system("mkdir intra/zheng_sorted")
    os.system("cp FuEtAll_benchmark/Datasets/Zhengsort/Zhengsort_intra_validation.R intra/zheng_sorted")
    df = pd.read_csv("FuEtAll_benchmark/Datasets/Zhengsort/Filtered_DownSampled_SortedPBMC_data.csv",header=0,index_col=0)
    labs = pd.read_csv("FuEtAll_benchmark/Datasets/Zhengsort/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("intra/zheng_sorted/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("intra/zheng_sorted/Labels.csv",sep=",",header=True,index=False)

def inter():
    ls = ["10Xv2A_pbmc1","CL2_pbmc1","SM2_pbmc1","10Xv2B_pbmc1","DR_pbmc1","SW_pbmc1","10Xv3_pbmc1","iD_pbmc1"]
    os.system("mkdir ./inter")
    for e in ls:
        os.system("mkdir inter/"+e)
        F_counts = e+".csv"
        F_labels = e+"Labels.csv"
        df = pd.read_csv("FuEtAll_benchmark/Datasets/PBMCbench/"+e+"/"+F_counts,header=0,index_col=0)
        labs = pd.read_csv("FuEtAll_benchmark/Datasets/PBMCbench/"+e+"/"+F_labels,header=0,index_col=False)
        tdf = df.T
        labs.columns = ["CELL_ANNOTATION"]
        tdf.to_csv("inter/"+e+"/read_counts.tsv",sep="\t",header=True,index=True)
        labs.to_csv("inter/"+e+"/labels.tsv",sep="\t",header=True,index=False)


if __name__ == "__main__":
    #transpositionIntra = intra()
    transpositionInter = inter()

    #os.system("mkdir ./processedData")
    os.system("mv inter/ processedData/")
    #os.system("mv intra/ processedData/")

end_time = datetime.now()
print('Duration files transposition: {}'.format(end_time - start_time))
