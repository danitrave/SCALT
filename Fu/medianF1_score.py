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

def table_parser(T1,F1,ref):
    df = pd.read_csv(T1,sep="\t",header=0,index_col=0)
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
    surv = pd.read_csv(F1,sep="\t",header=0,index_col=0)
    for i in range(len(INDEXES)):
        total_cells += 1
        rf = dfref.iloc[i,0].replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper()
        if rf not in F1_score:
            continue
        alt = list(df.iloc[i,:])
        minimum = min(alt)
        if minimum > 0.05:
            unassigned += 1
            continue
        zipped = list(zip(COLS,alt))
        sortedZipped = sorted(zipped,key=itemgetter(1),reverse=False)  #False if it is the p-value matrix or likelihood
        upperSortedZipped = [(e[0].replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper(),e[1]) for e in sortedZipped]
        retain_significant = significance_validator(upperSortedZipped)
        retain_significant_annotation = [k[0] for k in retain_significant]
        #print(rf,retain_significant_annotation)
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

    median_top = []
    median_sign = []
    up_F1_score = {}
    ### check coherence between references and annotatable classes available 
    for k in F1_score:
        if k not in references_classes:
            continue
        else:
            up_F1_score[k]=F1_score[k]

    report = open("F1_scoreReport.txt","w")
    report.write("The total number of cells is: "+str(total_cells)+"\n")
    report.write("The number of cells unclassified is: "+str(unassigned)+"\n")

    for cell in up_F1_score: 
        if ((2*up_F1_score[cell]["TP_best"])+up_F1_score[cell]["FP_best"]+up_F1_score[cell]["FN_best"])==0:
            F1B = 0.0
        else:
            F1B = (2*up_F1_score[cell]["TP_best"])/((2*up_F1_score[cell]["TP_best"])+up_F1_score[cell]["FP_best"]+up_F1_score[cell]["FN_best"])
        if ((2*up_F1_score[cell]["TP_sign"])+up_F1_score[cell]["FP_sign"]+up_F1_score[cell]["FN_sign"])==0:
            F1T = 0.0
        else:
            F1T = (2*up_F1_score[cell]["TP_sign"])/((2*up_F1_score[cell]["TP_sign"])+up_F1_score[cell]["FP_sign"]+up_F1_score[cell]["FN_sign"])
        median_top += [F1B,]
        median_sign += [F1T,]
        report.write("Cell type: "+cell+" F1 score best: "+str(F1B)+"; F1 score significant: "+str(F1T)+"\n")

    report.write("The median best F1 score is: "+str(np.median(median_top))+"\n")
    report.write("The median significant F1 score is: "+ str(np.median(median_sign))+"\n")

    report.close()

if __name__ == "__main__":
    t = sys.argv[1]     #table with all scores
    f = sys.argv[2]     #table that says if a cell has more than 500 expressed or not
    a = sys.argv[3]     #reference annotation
    COUNTS = table_parser(t,f,a)

end_time = datetime.now()
print('Duration accuracy calculation: {}'.format(end_time - start_time))

