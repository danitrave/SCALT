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

#### cell lines ####
os.system("mkdir INTER_benchmark/cell_lines")
os.system("cp processedData/inter/cell_lines/* ./")
os.system("Rscript --vanilla RdataParser.R CV_folds.RData")
dfcoord = pd.read_csv("coordinates.tsv",sep="\t",header=0,index_col=0)
counts = pd.read_csv("read_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("labels.tsv",sep="\t",header=0,index_col=None)
cl = {"10x_5cl":(dfcoord.iloc[0,0],dfcoord.iloc[0,1],dfcoord.iloc[0,2],dfcoord.iloc[0,3]),
    "CelSeq2_5cl":(dfcoord.iloc[1,0],dfcoord.iloc[1,1],dfcoord.iloc[1,2],dfcoord.iloc[1,3])}
for c in cl:
    os.system("mkdir INTER_benchmark/cell_lines/training_"+c)
    trainDS = counts.iloc[:,cl[c][0]:cl[c][1]]
    trainLABS = labels.iloc[cl[c][0]:cl[c][1],:]
    testDS = counts.iloc[:,cl[c][2]:cl[c][3]]
    testLABS = labels.iloc[cl[c][2]:cl[c][3],:]
    trainDS.to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    trainLABS.to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    testDS.to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    testLABS.to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTER_benchmark/cell_lines/training_"+c)

os.system("rm coordinates.tsv")
os.system("rm read_counts.tsv")
os.system("rm labels.tsv")
os.system("rm CV_folds.RData")

#### Pancreas ####
os.system("mkdir INTER_benchmark/pancreas")
os.system("cp processedData/inter/pancreas/* ./")
total_counts = pd.read_csv("read_counts.tsv",sep="\t",header=0,index_col=0)
total_labels = pd.read_csv("labels.tsv",sep="\t",header=0,index_col=None)
dts = ["Baron_Human","Muraro","Segerstolpe"]
cells = {}
for d in dts:
    counts = pd.read_csv("processedData/intra/pancreas/"+d+"/read_counts.tsv",sep="\t",header=0,index_col=0)
    labels = pd.read_csv("processedData/intra/pancreas/"+d+"/labels.tsv",sep="\t",header=0,index_col=None)
    desideredNames = []
    names = list(counts.columns)
    for i in range(len(list(labels["CELL_ANNOTATION"]))):
        if labels.iloc[i,0] in ["alpha", "beta", "delta", "gamma"]:
            desideredNames += [names[i],]
    indexes = [total_counts.columns.get_loc(col) for col in desideredNames]
    extracted_counts = total_counts.loc[:,desideredNames]
    extracted_labels = total_labels.iloc[indexes,:]
    extracted_counts.to_csv(d+"_counts.tsv",sep="\t",header=True,index=True)
    extracted_labels.to_csv(d+"_labels.tsv",sep="\t",header=True,index=False)

os.system("rm read_counts.tsv")
os.system("rm labels.tsv")
os.system("rm CV_folds.RData")

combinations = {("Baron_Human_counts.tsv","Baron_Human_labels.tsv"):("Muraro_counts.tsv","Muraro_labels.tsv","Segerstolpe_counts.tsv","Segerstolpe_labels.tsv"),("Muraro_counts.tsv","Muraro_labels.tsv"):("Baron_Human_counts.tsv","Baron_Human_labels.tsv","Segerstolpe_counts.tsv","Segerstolpe_labels.tsv"),("Segerstolpe_counts.tsv","Segerstolpe_labels.tsv"):("Baron_Human_counts.tsv","Baron_Human_labels.tsv","Muraro_counts.tsv","Muraro_labels.tsv")}

for combi in combinations:
    os.system("mkdir INTER_benchmark/pancreas/testing_"+combi[0].split("_counts")[0])    
    A = pd.read_csv(combinations[combi][0],sep="\t",header=0,index_col=0)
    B = pd.read_csv(combinations[combi][2],sep="\t",header=0,index_col=0)
    LA = pd.read_csv(combinations[combi][1],sep="\t",header=0,index_col=None)
    LB = pd.read_csv(combinations[combi][3],sep="\t",header=0,index_col=None)
    trainCounts = pd.concat([A,B],axis=1)
    trainLabels = pd.concat([LA,LB],axis=0)
    trainCounts.to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    trainLabels.to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    os.system("cp "+combi[0]+" benchmark_counts.tsv")
    os.system("cp "+combi[1]+" benchmark_labels.tsv")
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTER_benchmark/pancreas/testing_"+combi[0].split("_counts")[0])

os.system("rm Baron_Human_counts.tsv Baron_Human_labels.tsv Muraro_counts.tsv Muraro_labels.tsv Segerstolpe_counts.tsv Segerstolpe_labels.tsv")

#### PBMCs - experiment one ####
#datasets = sorted(os.listdir("processedData/inter/pbmc"))
#os.system("mkdir INTER_benchmark/pbmc")
#for i in datasets:
    #os.system("cp processedData/inter/pbmc/"+i+"/* ./") 
    #files = os.listdir("processedData/inter/pbmc/"+i)
    #os.system("mkdir INTER_benchmark/pbmc/"+i+"_training")
    #rdata = ""
    #for f in files:
        #if f.endswith(".RData"):
            #rdata = f
    #os.system("Rscript --vanilla RdataParser.R "+rdata)
    #dfcoord = pd.read_csv("coordinates.tsv",sep="\t",header=0,index_col=0)
    #counts = pd.read_csv("read_counts.tsv",sep="\t",header=0,index_col=0)
    #labels = pd.read_csv("labels.tsv",sep="\t",header=0,index_col=None)
    #os.system("rm coordinates.tsv read_counts.tsv labels.tsv "+rdata)
    #trainCounts = counts.iloc[:,dfcoord.iloc[0,0]:dfcoord.iloc[0,1]]
    #trainLabels = labels.iloc[dfcoord.iloc[0,0]:dfcoord.iloc[0,1],:]
    #testCounts = counts.iloc[:,dfcoord.iloc[0,2]:dfcoord.iloc[0,3]]
    #testLabels = labels.iloc[dfcoord.iloc[0,2]:dfcoord.iloc[0,3],:]
    #trainCounts.to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    #trainLabels.to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    #testCounts.to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    #testLabels.to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    #os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    #os.system("mkdir INTER_benchmark/pbmc/"+i+"_training/testDataset")
    #os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    #os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    #os.system("mv results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTER_benchmark/pbmc/"+i+"_training/testDataset")
    #for j in datasets:
        #if j==i:
            #continue
        #else:
            #os.system("cp processedData/inter/pbmc/"+j+"/read_counts.tsv benchmark_counts.tsv")
            #os.system("cp processedData/inter/pbmc/"+j+"/labels.tsv benchmark_labels.tsv")
            #os.system("mkdir INTER_benchmark/pbmc/"+i+"_training/"+j+"_testing")
            #os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
            #os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
            #os.system("mv results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTER_benchmark/pbmc/"+i+"_training/"+j+"_testing")
    #os.system("mv AnnolistsBuilder_results/ custom/ INTER_benchmark/pbmc/"+i+"_training")

#### PBMCs - experiment two ####
os.system("mkdir INTER_benchmark/pbmc")
os.system("cp originalFiles/Inter-dataset/PbmcBench/Statisitics.xlsx ./")

df = pd.read_excel("Statisitics.xlsx", sheet_name="Sheet2")
sliced = {"Smart-Seq2":df.iloc[0:7,:],"10Xv2":df.iloc[8:15,:],"10Xv3":df.iloc[16:22,:],"CEL-Seq":df.iloc[23:30,:],
"Drop-Seq":df.iloc[31:38,:],"inDrop":df.iloc[39:46,:],"Seq-Well":df.iloc[47:54,:]}

#sliced = {"inDrop":df.iloc[39:46,:]}

for i in sliced:
    os.system("cp processedData/inter/pbmc/"+i+"/read_counts.tsv ./")
    os.system("cp processedData/inter/pbmc/"+i+"/labels.tsv ./")
    data = pd.read_csv("read_counts.tsv",sep="\t",header=0,index_col=0)
    labs = pd.read_csv("labels.tsv",sep="\t",header=0,index_col=None)

    #### Training #####
    os.system("mkdir INTER_benchmark/pbmc/"+i+"_training")
    startTrain = int(sliced[i].iloc[0,0].split("pbmc")[1])-1
    endTrain = int(sliced[i].iloc[1,0])-1
    trainCounts = data.iloc[:,startTrain:endTrain]
    trainLabs = labs.iloc[startTrain:endTrain,:]
    trainCounts.to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    trainLabs.to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")

	#### Testing ####
    testNames = list(sliced[i].iloc[:,1])
    testIndex = [(int(x.split(":")[0])-1,int(x.split(":")[1])-1) for x in list(sliced[i].iloc[:,3])]
    z = list(zip(testNames,testIndex))
    for t in z:
        os.system("mkdir INTER_benchmark/pbmc/"+i+"_training/"+t[0].replace(".","_")+"_testing")
        testCounts = data.iloc[:,t[1][0]:t[1][1]]
        testLabs = labs.iloc[t[1][0]:t[1][1],:]
        testCounts.to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
        testLabs.to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
        os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
        os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
        os.system("mv results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTER_benchmark/pbmc/"+i+"_training/"+t[0].replace(".","_")+"_testing")

    os.system("rm read_counts.tsv labels.tsv")
    os.system("mv AnnolistsBuilder_results/ custom/ INTER_benchmark/pbmc/"+i+"_training")

os.system("rm Statisitics.xlsx")

end_time = datetime.now()
print('Duration inter-experimental set up benchmark: {}'.format(end_time - start_time))
