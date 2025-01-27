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

def removeSmallPop(counts,labels,toremove):
    keep = []
    for c in range(len(list(labels["CELL_ANNOTATION"]))):
        if labels.iloc[c,0] not in toremove:
            keep += [c,]
        else:
            continue

    labels_keep = labels.iloc[keep,:]
    counts_keep = counts.iloc[:,keep]
    return counts_keep,labels_keep

def identifySmallPop(A):
    occurences = {}
    for i in list(A["CELL_ANNOTATION"]):
        if i not in occurences:
            occurences[i]=1
        else:
            occurences[i]+=1
    annoToeliminate = []
    for x in occurences:
        if occurences[x]<100:
            annoToeliminate+=[x,]
    return annoToeliminate

os.system("mkdir ./INTRA_benchmark")

#### cell lines ####
os.system("mkdir INTRA_benchmark/cellLines")
## 10x_5cl ##
os.system("mkdir INTRA_benchmark/cellLines/10x_5cl")
os.system("cp processedData/intra/cell_lines/10x_5cl/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/cell_lines/10x_5cl/labels.tsv ./benchmark_labels.tsv")
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/cellLines/10x_5cl/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/cellLines/10x_5cl/"+fold)

os.system("mv CV_folds.RData INTRA_benchmark/cellLines/10x_5cl")

## CelSeq2_5cl ##
os.system("mkdir INTRA_benchmark/cellLines/CelSeq2_5cl")
os.system("cp processedData/intra/cell_lines/CelSeq2_5cl/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/cell_lines/CelSeq2_5cl/labels.tsv ./benchmark_labels.tsv")
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/cellLines/CelSeq2_5cl/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation ensembl_id")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation ensembl_id")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/cellLines/CelSeq2_5cl/"+fold)
os.system("mv CV_folds.RData INTRA_benchmark/cellLines/CelSeq2_5cl")

#### Pancreas ####
os.system("mkdir INTRA_benchmark/pancreas")
## Baron_Human  ##
os.system("mkdir INTRA_benchmark/pancreas/Baron_Human")
os.system("cp processedData/intra/pancreas/Baron_Human/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/pancreas/Baron_Human/labels.tsv ./benchmark_labels.tsv")
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/pancreas/Baron_Human/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/pancreas/Baron_Human/"+fold)

os.system("mv CV_folds.RData INTRA_benchmark/pancreas/Baron_Human")

## Segerstolpe  ##
os.system("mkdir INTRA_benchmark/pancreas/segerstolpe")
os.system("cp processedData/intra/pancreas/Segerstolpe/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/pancreas/Segerstolpe/labels.tsv ./benchmark_labels.tsv")
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/pancreas/segerstolpe/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/pancreas/segerstolpe/"+fold)
os.system("mv CV_folds.RData INTRA_benchmark/pancreas/segerstolpe")

## Muraro  ##
os.system("mkdir INTRA_benchmark/pancreas/muraro")
os.system("cp processedData/intra/pancreas/Muraro/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/pancreas/Muraro/labels.tsv ./benchmark_labels.tsv")
df = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
newGenes = []
oldGenes = list(df.index)
for g in oldGenes:
    newGenes += [g.split("__")[0],]
df.index = newGenes
os.system("rm benchmark_counts.tsv")
df.to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/pancreas/muraro/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/pancreas/muraro/"+fold)
os.system("mv CV_folds.RData INTRA_benchmark/pancreas/muraro")

#### Zheng 68k ####
os.system("mkdir INTRA_benchmark/zheng68k")
os.system("cp processedData/intra/zheng68k/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/zheng68k/labels.tsv ./benchmark_labels.tsv")
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/zheng68k/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/zheng68k/"+fold)
os.system("mv CV_folds.RData INTRA_benchmark/zheng68k")

### Zheng sorted ####
os.system("mkdir INTRA_benchmark/zheng_sorted")
os.system("cp processedData/intra/zheng_sorted/read_counts.tsv ./benchmark_counts.tsv")
os.system("cp processedData/intra/zheng_sorted/labels.tsv ./benchmark_labels.tsv")
os.system("Rscript --vanilla cv_folds.R benchmark_labels.tsv 1 ./")
os.system("Rscript --vanilla RdataParserCrossValidation.R CV_folds.RData")
data = pd.read_csv("benchmark_counts.tsv",sep="\t",header=0,index_col=0)
labels = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
whichToEliminate = identifySmallPop(labels)
os.system("rm benchmark_counts.tsv benchmark_labels.tsv")
foldcrossValidation = {"CV_1":("1_crossValidationTrain.tsv","1_crossValidationTest.tsv"),"CV_2":("2_crossValidationTrain.tsv","2_crossValidationTest.tsv"),
        "CV_3":("3_crossValidationTrain.tsv","3_crossValidationTest.tsv"),"CV_4":("4_crossValidationTrain.tsv","4_crossValidationTest.tsv"),
        "CV_5":("5_crossValidationTrain.tsv","5_crossValidationTest.tsv")}
for fold in foldcrossValidation:
    os.system("mkdir INTRA_benchmark/zheng_sorted/"+fold)
    trainCells = pd.read_csv(foldcrossValidation[fold][0],sep="\t",header=0,index_col=0)
    testCells = pd.read_csv(foldcrossValidation[fold][1],sep="\t",header=0,index_col=0)
    trainDF = data.iloc[:,list(trainCells["train"])]
    trainLABELS = labels.iloc[list(trainCells["train"]),:]
    filterTrain = removeSmallPop(trainDF,trainLABELS,whichToEliminate)
    testDF = data.iloc[:,list(testCells["test"])]
    testLABELS = labels.iloc[list(testCells["test"]),:]
    filterTest = removeSmallPop(testDF,testLABELS,whichToEliminate)
    filterTrain[0].to_csv("train_counts.tsv",sep="\t",header=True,index=True)
    filterTrain[1].to_csv("train_labels.tsv",sep="\t",header=True,index=False)
    filterTest[0].to_csv("benchmark_counts.tsv",sep="\t",header=True,index=True)
    filterTest[1].to_csv("benchmark_labels.tsv",sep="\t",header=True,index=False)
    os.system("python3 SCALT_AnnotaionListsBuilder.py train_counts.tsv train_labels.tsv -Notation gene_symbol")
    os.system("python3 SCALT.py benchmark_counts.tsv -Types custom -Notation gene_symbol")
    os.system("python3 medianF1_score.py results_directory/p_values.tsv results_directory/benchmark_counts_adj_genesExpressed_filter.tsv benchmark_labels.tsv")
    os.system("mv "+foldcrossValidation[fold][0]+" "+foldcrossValidation[fold][1]+" AnnolistsBuilder_results/ custom/ results_directory/ REPORT.html benchmark_labels.tsv F1_scoreReport.txt INTRA_benchmark/zheng_sorted/"+fold)
os.system("mv CV_folds.RData INTRA_benchmark/zheng_sorted")

end_time = datetime.now()
print('Duration intra-experimental set up benchmark: {}'.format(end_time - start_time))
