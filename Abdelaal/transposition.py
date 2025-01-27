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
    #### Cell lines ####
    os.system("mkdir intra/cell_lines")
    paths = ["10x_5cl","CelSeq2_5cl"]
    for p in paths:
        os.system("mkdir intra/cell_lines/"+p)
        df = pd.read_csv("originalFiles/Intra-dataset/CellBench/"+p+"/"+p+"_data.csv",header=0,index_col=0)
        labs = pd.read_csv("originalFiles/Intra-dataset/CellBench/"+p+"/Labels.csv",header=0,index_col=False)
        tdf = df.T
        labs.columns = ["CELL_ANNOTATION"]
        tdf.to_csv("intra/cell_lines/"+p+"/read_counts.tsv",sep="\t",header=True,index=True)
        labs.to_csv("intra/cell_lines/"+p+"/labels.tsv",sep="\t",header=True,index=False)

    #### Pancreas ####
    os.system("mkdir intra/pancreas")
    paths = [("originalFiles/Intra-dataset/Pancreatic_data/Baron Human/Filtered_Baron_HumanPancreas_data.csv","originalFiles/Intra-dataset/Pancreatic_data/Baron Human/Labels.csv"),("originalFiles/Intra-dataset/Pancreatic_data/Muraro/Filtered_Muraro_HumanPancreas_data.csv","originalFiles/Intra-dataset/Pancreatic_data/Muraro/Labels.csv"),("originalFiles/Intra-dataset/Pancreatic_data/Segerstolpe/Filtered_Segerstolpe_HumanPancreas_data.csv","originalFiles/Intra-dataset/Pancreatic_data/Segerstolpe/Labels.csv")] 
    for p in paths:
        df = pd.read_csv(p[0],header=0,index_col=0)
        labs = pd.read_csv(p[1],header=0,index_col=False)
        tdf = df.T
        labs.columns = ["CELL_ANNOTATION"]
        nameDir = p[0].split("/")[3].replace(" ","_")
        os.system("mkdir intra/pancreas/"+nameDir)
        tdf.to_csv("intra/pancreas/"+nameDir+"/read_counts.tsv",sep="\t",header=True,index=True)
        labs.to_csv("intra/pancreas/"+nameDir+"/labels.tsv",sep="\t",header=True,index=False)
    
    #### Zheng 68k ####
    os.system("mkdir intra/zheng68k")
    df = pd.read_csv("originalFiles/Intra-dataset/Zheng 68K/Filtered_68K_PBMC_data.csv",header=0,index_col=0)
    labs = pd.read_csv("originalFiles/Intra-dataset/Zheng 68K/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("intra/zheng68k/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("intra/zheng68k/labels.tsv",sep="\t",header=True,index=False)

    #### Zheng sorted ####
    os.system("mkdir intra/zheng_sorted")
    df = pd.read_csv("originalFiles/Intra-dataset/Zheng sorted/Filtered_DownSampled_SortedPBMC_data.csv",header=0,index_col=0)
    labs = pd.read_csv("originalFiles/Intra-dataset/Zheng sorted/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("intra/zheng_sorted/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("intra/zheng_sorted/labels.tsv",sep="\t",header=True,index=False)

def inter():
    os.system("mkdir ./inter")
    #### cell line ####
    os.system("mkdir inter/cell_lines")
    df = pd.read_csv("originalFiles/Inter-dataset/CellBench/Combined_10x_CelSeq2_5cl_data.csv",header=0,index_col=0)
    labs = pd.read_csv("originalFiles/Inter-dataset/CellBench/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("inter/cell_lines/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("inter/cell_lines/labels.tsv",sep="\t",header=True,index=False)
    os.system("cp originalFiles/Inter-dataset/CellBench/CV_folds.RData inter/cell_lines")

    #### Pancreas ####
    os.system("mkdir inter/pancreas")
    df = pd.read_csv("originalFiles/Inter-dataset/Pancreatic/Combined_HumanPancreas_data.csv",header=0,index_col=0)
    labs = pd.read_csv("originalFiles/Inter-dataset/Pancreatic/Labels.csv",header=0,index_col=False)
    tdf = df.T
    labs.columns = ["CELL_ANNOTATION"]
    tdf.to_csv("inter/pancreas/read_counts.tsv",sep="\t",header=True,index=True)
    labs.to_csv("inter/pancreas/labels.tsv",sep="\t",header=True,index=False)
    os.system("cp originalFiles/Inter-dataset/Pancreatic/CV_folds.RData inter/pancreas")

    #### PBMCs ####
    os.system("mkdir inter/pbmc")
    platf = ["10Xv2","Drop-Seq","Smart-Seq2","10Xv3","inDrop","CEL-Seq","Seq-Well"]
    for x in platf:
        os.system("mkdir inter/pbmc/"+x)
        files = os.listdir("originalFiles/Inter-dataset/PbmcBench/"+x)
        counts = ""
        labels = ""
        Rdata = ""
        for i in files:
            if i.endswith("pbmc1.csv"):
                counts = i
            elif i.endswith("pbmc1Labels.csv"):
                labels = i
            elif i.endswith("folds.RData"):
                Rdata = i

        df = pd.read_csv("originalFiles/Inter-dataset/PbmcBench/"+x+"/"+counts,header=0,index_col=0)
        labs = pd.read_csv("originalFiles/Inter-dataset/PbmcBench/"+x+"/"+labels,header=0,index_col=False)
        tdf = df.T
        labs.columns = ["CELL_ANNOTATION"]
        tdf.to_csv("inter/pbmc/"+x+"/read_counts.tsv",sep="\t",header=True,index=True)
        labs.to_csv("inter/pbmc/"+x+"/labels.tsv",sep="\t",header=True,index=False)
        os.system("cp originalFiles/Inter-dataset/PbmcBench/"+x+"/"+Rdata+" inter/pbmc/"+x)

if __name__ == "__main__":
    transpositionIntra = intra()
    transpositionInter = inter()

    os.system("mkdir ./processedData")
    os.system("mv inter/ processedData/")
    os.system("mv intra/ processedData/")

end_time = datetime.now()
print('Duration files transposition: {}'.format(end_time - start_time))
