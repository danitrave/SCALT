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

def significance_validator(R):
    L = []
    for j in R:
        if j[1] < 0.05:     #p-value significance validation
            L += [j,]
    return L

def classification(T1,F1):
    annotation = []
    df = pd.read_csv(T1,sep="\t",header=0,index_col=0)    #table reporting the p-values from the likelihood-ratio test
    INDEXES = list(df.index)
    COLS = list(df.columns)
    surv = pd.read_csv(F1,sep="\t",header=0,index_col=0)   #table reporting the PASS-EXCLUDE notation based on the survival or not to the genes-expressed filter
    cellTypes_counts = {"Unclassified":0}
    for i in range(len(INDEXES)):
        if surv.iloc[i,0]=="EXCLUDE":              #cells not passing the genes-expressed filer are annotated as "unclassified"
            cellTypes_counts["Unclassified"]+=1
            annotation+=["Unclassified",]
            continue
        else:
            alt = list(df.iloc[i,:])
            minimum = min(alt)
            if minimum > 0.05:
                cellTypes_counts["Unclassified"]+=1     #cells having the lowest p-value greater than 0.05 are annotated as "unclassified"
                annotation+=["Unclassified",]
                continue
            zipped = list(zip(COLS,alt))
            sortedZipped = sorted(zipped,key=itemgetter(1),reverse=False)  #sort annotations by p-values in increasing order
            upperSortedZipped = [(e[0].replace("_",".").replace(" ",".").replace("-","."),e[1]) for e in sortedZipped]
            retain_significant = significance_validator(upperSortedZipped)       #collect only significant annotations
            retain_significant_annotation = [k[0] for k in retain_significant]
            annotation+=[retain_significant_annotation[0].replace("."," "),]
            if retain_significant_annotation[0].replace("."," ") not in cellTypes_counts:   #count the number of cells annotated to a specific cell type and collect the numbers into a dictionary having the cell type cathegory as keys
                cellTypes_counts[retain_significant_annotation[0].replace("."," ")]=1
            else:
                cellTypes_counts[retain_significant_annotation[0].replace("."," ")]+=1
    df1 = pd.DataFrame.from_dict({"SCALTclassification":annotation})
    df2 = pd.DataFrame.from_dict(cellTypes_counts,orient='index')
    df2.columns = ["Counts"]
    df1.to_csv("SCALT_classification.tsv",sep="\t",header=True,index=False)
    df2.to_csv("grouppedCellTypes.tsv",sep="\t",header=True,index=True)

os.system("mkdir ./classify_withSCALT")

#### Cell lines ####
samples = ["10x_5cl","CelSeq2_5cl"]
for i in samples:
    os.system("mkdir classify_withSCALT/"+i)
    os.system("cp processedData/intra/cell_lines/"+i+"/read_counts.tsv ./")
    os.system("cp processedData/intra/cell_lines/"+i+"/labels.tsv ./")
    if i == "CelSeq2_5cl":
        os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation ensembl_id -Granularity high")
        dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
        os.system("mv labels.tsv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/"+i)
    else:
        os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation gene_symbol -Granularity high") 
        dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
        os.system("mv labels.tsv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/"+i)

#### Pancreas ####
samples = ["Baron_Human","Muraro","Segerstolpe"]
for i in samples:
    os.system("mkdir classify_withSCALT/"+i)
    os.system("cp processedData/intra/pancreas/"+i+"/read_counts.tsv ./")
    os.system("cp processedData/intra/pancreas/"+i+"/labels.tsv ./")
    if i == "Muraro":  
        df = pd.read_csv("read_counts.tsv",sep="\t",header=0,index_col=0)
        genes = [e.split("__")[0] for e in list(df.index)]
        df.index=genes
        os.system("rm read_counts.tsv")
        df.to_csv("read_counts.tsv",sep="\t",header=True,index=True)
        os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation gene_symbol -Granularity high")
        dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
        os.system("mv labels.tsv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/"+i)
    else:
        os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation gene_symbol -Granularity high")
        dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
        os.system("mv labels.tsv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/"+i)

#### Zheng 68K ###
os.system("mkdir classify_withSCALT/zheng68k")
os.system("cp processedData/intra/zheng68k/read_counts.tsv ./")
os.system("cp processedData/intra/zheng68k/labels.tsv ./")
os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation gene_symbol -Granularity high")
dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
os.system("mv labels.tsv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/zheng68k")

#### Zheng sorted ###
os.system("mkdir classify_withSCALT/zheng_sorted")
os.system("cp processedData/intra/zheng_sorted/read_counts.tsv ./")
os.system("cp processedData/intra/zheng_sorted/labels.tsv ./")
os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation gene_symbol -Granularity high")
dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
os.system("mv labels.tsv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/zheng_sorted")

#### PBMCs ####
#os.system("mkdir classify_withSCALT/pbmc")
#samples = ["10Xv2","10Xv3","CEL-Seq","Drop-Seq","inDrop","Seq-Well","Smart-Seq2"]
#for i in samples:
    #os.system("mkdir classify_withSCALT/pbmc/"+i)
    #os.system("cp processedData/inter/pbmc/"+i+"/read_counts.tsv ./")
    #os.system("python3 SCALT.py read_counts.tsv -CPUs 10 -Notation gene_symbol -Granularity low")
    #dfs = classification("results_directory/p_values.tsv","results_directory/read_counts_adj_genesExpressed_filter.tsv")
    #os.system("mv SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/pbmc/"+i)

#### PBMCs version two ####
os.system("cp originalFiles/Inter-dataset/PbmcBench/Statisitics.xlsx ./")

df = pd.read_excel("Statisitics.xlsx", sheet_name="Sheet2")
sliced = {"Smart-Seq2":df.iloc[0:7,:],"10Xv2":df.iloc[8:15,:],"10Xv3":df.iloc[16:22,:],"CEL-Seq":df.iloc[23:30,:],
"Drop-Seq":df.iloc[31:38,:],"inDrop":df.iloc[39:46,:],"Seq-Well":df.iloc[47:54,:]}
datasets = {}

for i in sliced:
    os.system("cp processedData/inter/pbmc/"+i+"/read_counts.tsv ./")
    os.system("cp processedData/inter/pbmc/"+i+"/labels.tsv ./")
    data = pd.read_csv("read_counts.tsv",sep="\t",header=0,index_col=0)
    labs = pd.read_csv("labels.tsv",sep="\t",header=0,index_col=None)
    startTrain = int(sliced[i].iloc[0,0].split("pbmc")[1])-1
    endTrain = int(sliced[i].iloc[1,0])-1
    trainCounts = data.iloc[:,startTrain:endTrain]
    trainLabs = labs.iloc[startTrain:endTrain,:]
    trainCounts.to_csv(i+"pbmc1_counts.tsv",sep="\t",header=True,index=True)
    trainLabs.to_csv(i+"pbmc1_labels.tsv",sep="\t",header=True,index=False)
    datasets[i+"pbmc1"]=(i+"pbmc1_counts.tsv",i+"pbmc1_labels.tsv")
    
    #### pbmc2 of the same sequencing protocol ####
    testNames = list(sliced[i].iloc[:,1])
    testIndex = [(int(x.split(":")[0])-1,int(x.split(":")[1])-1) for x in list(sliced[i].iloc[:,3])]
    z = list(zip(testNames,testIndex))
    for t in z:
        if "pbmc2" in t[0]:
            testCounts = data.iloc[:,t[1][0]:t[1][1]]
            testLabs = labs.iloc[t[1][0]:t[1][1],:]
            testCounts.to_csv(i+"pbmc2_counts.tsv",sep="\t",header=True,index=True)
            testLabs.to_csv(i+"pbmc2_labels.tsv",sep="\t",header=True,index=False)
            datasets[i+"pbmc2"]=(i+"pbmc2_counts.tsv",i+"pbmc2_labels.tsv")

    os.system("rm read_counts.tsv labels.tsv")
os.system("rm Statisitics.xlsx")

os.system("mkdir classify_withSCALT/pbmc")

for sample in datasets:
    os.system("mkdir classify_withSCALT/pbmc/"+sample)
    os.system("python3 SCALT.py "+datasets[sample][0]+" -CPUs 10 -Notation gene_symbol -Granularity high")
    dfs = classification("results_directory/p_values.tsv","results_directory/"+datasets[sample][0].split(".tsv")[0]+"_adj_genesExpressed_filter.tsv")
    os.system("mv "+datasets[sample][1]+" SCALT_classification.tsv grouppedCellTypes.tsv results_directory/ REPORT.html classify_withSCALT/pbmc/"+sample+"/")

end_time = datetime.now()
print('Duration classification with SCALT: {}'.format(end_time - start_time))
