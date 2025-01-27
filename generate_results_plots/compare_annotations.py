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

def annotation_tables(T1,ref):
    ORIGINAL = []
    ANNOTATED = []
    df = pd.read_csv(T1,sep="\t",header=0,index_col=0)
    INDEXES = list(df.index)
    dfref = pd.read_csv(ref,sep="\t",header=0,index_col=False)
    #references_classes = [x.replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper() for x in list(dfref["CELL_ANNOTATION"].unique())] 
    #annotation_classes = [x.replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper() for x in list(df.columns.unique())]
    COLS = list(df.columns)

    for i in range(len(INDEXES)):
        rf = dfref.iloc[i,0].replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper()
        ORIGINAL += [rf,]
        alt = list(df.iloc[i,:])
        minimum = min(alt)
        if minimum > 0.05:
            ANNOTATED += ["unassigned",]
            continue
        zipped = list(zip(COLS,alt))
        sortedZipped = sorted(zipped,key=itemgetter(1),reverse=False)  #False if it is the p-value matrix or likelihood
        upperSortedZipped = [(e[0].replace("_",".").replace(" ",".").replace("–","-").replace("/",".").replace("(",".").replace(")",".").upper(),e[1]) for e in sortedZipped]
        retain_significant = significance_validator(upperSortedZipped)
        retain_significant_annotation = [k[0] for k in retain_significant]

        ANNOTATED += [retain_significant_annotation[0],]
    
    return ORIGINAL,ANNOTATED

if __name__ == "__main__":
    t1 = sys.argv[1]
    t2 = sys.argv[2]
    res = annotation_tables(t1,t2)
    df1 = pd.DataFrame(res[0])
    df1.to_csv("original.csv",header=False, index=False)
    df2 = pd.DataFrame(res[1])
    df2.to_csv("annotated.csv",header=False, index=False)
    os.system("Rscript --vanilla evaluate.r")
    
end_time = datetime.now()
print('Duration accuracy calculation: {}'.format(end_time - start_time))

