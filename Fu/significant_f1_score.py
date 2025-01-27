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

def sig_medianF1(p):
    cvs = ["CV_1","CV_2","CV_3","CV_4","CV_5"]
    F1 = {}
    for cv in cvs:
        f = open(p+"/"+cv+"/F1_scoreReport.txt","r")
        for line in f:
            if line.startswith("Cell type"):
                l = line.strip("\n").split(":")
                value = float(l[-1].replace(" ",""))
                cell_type = l[1].split(" score")[0][1:]
                if cell_type not in F1:
                    F1[cell_type]=[value,]
                else:
                    F1[cell_type]+=[value,]
        f.close()
    mean_F1 = {}
    for c in F1:
        mean_F1[c]=np.mean(F1[c])
    median_F1 = np.median(list(mean_F1.values()))

def inter_sig_medianF1(path_dir):
    cvs = os.listdir(path_dir)
    if "custom" in cvs:
        cvs.remove("custom")
    if "AnnolistsBuilder_results" in cvs:
        cvs.remove("AnnolistsBuilder_results")
    for cv in cvs:
        f = open(path_dir+"/"+cv+"/F1_scoreReport.txt")
        F1 = {}
        for line in f:
            if line.startswith("Cell type"):
                l = line.strip("\n").split(":")
                value = float(l[-1].replace(" ",""))
                cell_type = l[1].split(" score")[0][1:]
                if value == 0.0:
                    continue
                if cell_type not in F1:
                    F1[cell_type]=[value,]
                else:
                    F1[cell_type]+=[value,]
        f.close()
        median_F1 = np.median(list(F1.values()))
        print(path_dir)
        print(cv)
        print(median_F1)


if __name__ == "__main__":
    path_dir = sys.argv[1]

    if "INTRA" in path_dir:
        run = sig_medianF1(path_dir)
    else:
        run = inter_sig_medianF1(path_dir)

