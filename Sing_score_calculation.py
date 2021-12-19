# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 09:59:28 2021

@author: Srinath
"""
#Script used for computing Epi score, Mes scores, HIF1 score and AMPK score 

import os 
import pandas
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from singscore.singscore import score
os.chdir('<path>') #enter path to directory containing the input gene expression files 
sigs_list = ['../Gene signatures/AMPK_geneset.txt', '../Gene signatures/HIF1_geneset.txt']  
num = len(sigs_list)
filelist = os.listdir()    

for i in range(len(filelist)): 
    data = pandas.read_csv(open(filelist[i], 'r'), header = 'infer', sep='\t')  #the gene expression file
    data = data.set_index(keys='Gene')  #Name of the first column in gene expression matrix set to 'Gene'. 
    GSE = (filelist[i].split('_'))[0]   #List containing GSEIDs
    for j in range(num): 
        sigs = pandas.read_csv(open(sigs_list[j], 'r'), header = 'infer') #Input gene set file. Contains only one column with the header 'Gene'
        sigs = list(sigs['Gene']) 
        scored_data_single = score(up_gene = sigs, sample=data, norm_method='theoretical', full_data=True) 
        fh = open('Data files/' + GSE + '_sing_scores_' + sigs_list[j], 'w') 
        fh.write(scored_data_single.to_csv(sep = '\t', header = True))   
        fh.close()