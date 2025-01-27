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
        if j[1] < 0.05:
            L += [j,]
    return L

def uncertainty_validator(ds,t):
    deltaOfdeltas = abs((ds[0])-(ds[1]))
    #print(ds[0],ds[1],deltaOfdeltas)
    if deltaOfdeltas > t: #9.23:
        return "unbiased"
    else:
        return "uncertain"

def annotationComparator(T1,ref,T3):
    ORIGINAL = []
    ANNOTATED = []
    df = pd.read_csv(T1,sep="\t",header=0,index_col=0)
    deltas = pd.read_csv(T3,sep="\t",header=0,index_col=0) ### deltas
    INDEXES = list(df.index)
    dfref = pd.read_csv(ref,sep="\t",header=0,index_col=False)
    references_classes = [x.replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper() for x in list(dfref["CELL_ANNOTATION"].unique())] 
    annotation_classes = [x.replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper() for x in list(df.columns.unique())]
    F1_score = {}
    for y in annotation_classes:
        F1_score[y]={"TP_best":0,"TP_sign":0,"FN_best":0,"FN_sign":0,"FP_best":0,"FP_sign":0}
    COLS = list(df.columns)
    total_cells = 0
    unassigned = 0

    for i in range(len(INDEXES)):
        total_cells += 1
        rf = dfref.iloc[i,0].replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper()
        ORIGINAL += [rf,]
        #if rf not in F1_score:
            #continue
        alt = list(df.iloc[i,:])
        deltas_cell = sorted(list(deltas.iloc[i,:])) ### deltas
        minimum = min(alt)
        if minimum > 0.05:
            ANNOTATED += ["unassigned",]
            unassigned += 1
            continue
        uncertainty = uncertainty_validator(deltas_cell,6.03)
        if uncertainty == "uncertain":
            ANNOTATED += ["multiassigned",]
            #multianno += 1
            continue
        zipped = list(zip(COLS,alt))
        sortedZipped = sorted(zipped,key=itemgetter(1),reverse=False)  #False if it is the p-value matrix or likelihood
        upperSortedZipped = [(e[0].replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper(),e[1]) for e in sortedZipped]
        retain_significant = significance_validator(upperSortedZipped)
        retain_significant_annotation = [k[0] for k in retain_significant]

        ANNOTATED += [retain_significant_annotation[0],]
        if rf == retain_significant_annotation[0]:
            F1_score[rf]["TP_best"] += 1
            if rf in retain_significant_annotation:
                F1_score[rf]["TP_sign"] += 1
        if rf != retain_significant_annotation[0]:
            F1_score[rf]["FN_best"] +=1
            F1_score[retain_significant_annotation[0]]["FP_best"]+=1
            if rf not in retain_significant_annotation:
                F1_score[rf]["FN_sign"] += 1
                for l in retain_significant_annotation[1:]:
                    F1_score[l]["FP_sign"] += 1
    return ORIGINAL,ANNOTATED

if __name__ == "__main__":

    path = sys.argv[1]   #directory path
    cv_folds = ["CV_1","CV_2","CV_3","CV_4","CV_5"]
    collected_original = []
    collected_annotated = []
    for cv in cv_folds:
        print("RESULTS FROM "+cv)
        path_original = path+"/"+cv+"/benchmark_labels.tsv"
        path_annotated = path+"/"+cv+"/results_directory/p_values.tsv"
        path_deltas = path+"/"+cv+"/results_directory/deltas.tsv"
        comparison = annotationComparator(path_annotated,path_original,path_deltas)
        #collected_original += comparison[0]
        #collected_annotated += comparison[1]
        collected_original = comparison[0]
        collected_annotated = comparison[1]
        DFO = pd.DataFrame.from_dict({"original":collected_original})
        DFA = pd.DataFrame.from_dict({"scalt":collected_annotated})
        DFO.to_csv("originalAnno.csv",header=True,index=False)
        DFA.to_csv("scaltAnno.csv",header=True,index=False)

        os.system("Rscript --vanilla evaluate.r "+path+" "+cv)
        os.system("rm originalAnno.csv scaltAnno.csv")
        #os.system("mv originalAnno.csv scaltAnno.csv "+path+"/")
    
end_time = datetime.now()
print('Duration accuracy calculation: {}'.format(end_time - start_time))

