#!/usr/bin/python

""" Libraries required  """

import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
import concurrent.futures
import json
import warnings
import statistics
from operator import itemgetter
import itertools
import seaborn as sns
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function probability_calculation(L) calculates the mean of the probabilite of a gene to be expressed in a cell type starting from the collection
	of probabilites tables calculated from the boostrap pipeline.'''

def probability_calculation(L):
    m = sum(L)/len(L)                #mean
    m = m.replace(1.0, 0.9999)          #adjust if the probability is 1 because it is an estimate of the probability
    m.to_csv("genesCellTypes_probabilities.tsv",sep="\t")   #save the probability table
    return m
    
''' The function general_probabilities(R) calculates the mean of the probabilite of a gene to be expressed in general starting from the collection
	of probabilites tables calculated from the boostrap pipeline.'''

def general_probabilities(R):
    g = sum(R)/len(R)                 #mean
    g.to_csv("genesGeneral_probabilities.tsv",sep="\t")   #save the probability table
    return g
    
''' The function probabilities_ratio(ctp,gp) calculates the mean of the ratio of proabilities previously described starting from the collection
	of probabilites tables calculated from the boostrap pipeline.'''

def probabilities_ratio(ctp,gp):
    cells = list(ctp.columns)
    for ct in cells:
        DIV = ctp.loc[:,ct].div(gp["probs"],axis=0,fill_value=float(0))  #P(G|CT)/G(G)
        ctp[ct]=DIV
    ctp = ctp.fillna(float(0))       #change NA to 0
    ctp.to_csv("genesProbabilities_ratios.tsv",sep="\t",index=True,header=True)   #save the probabilites ratio table

if __name__ == "__main__":
    samples = os.listdir("./boostraps_samples")
    dfs = []
    general = []
    for i in samples:     #collect the tables from each boostrap sample
        p = pd.read_csv("./boostraps_samples/"+i+"/probabilities_tables/cell_type_probabilities.tsv",sep="\t",header=0,index_col=0)  #cell type probability table
        dfs += [p,]
        g = pd.read_csv("./boostraps_samples/"+i+"/probabilities_tables/global_probabilities.tsv",sep="\t",header=0,index_col=0)	#general probability table
        general += [g,]
    cell_type_probabilites=probability_calculation(dfs)
    general_probabilities=general_probabilities(general)
    probabilities_ratio(cell_type_probabilites,general_probabilities)

    #### Extract the total list of genes needed for further analysys ####

    genes = list(general_probabilities.index)
    dfg = pd.DataFrame.from_dict({"genes":genes})
    dfg.to_csv("TABLE_OF_GENES.tsv",sep="\t")

end_time = datetime.now()
print('Duration probaility calculation: {}'.format(end_time - start_time))
