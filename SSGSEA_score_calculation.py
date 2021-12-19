# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 15:19:08 2020

@author: Srinath
"""
#Script used for computing Hallmark EMT score, Hallmark glycolysis and OXPHOS scores

import os
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt 

os.chdir(<path>) #enter path to directory containing the input gene expression files
FileList = os.listdir() 
GSEID = [element.split('_')[0] for element in FileList]
num = len(GSEID) 

for i in range(num): 
    filename = FileList[i]
    GMT_file = '../Gene signatures/Glycolysis_signature.gmt'
    ss = gp.ssgsea(data= filename,
         gene_sets= GMT_file,
         outdir= '../../Output/GSEA_scores/' + GSEID[i], 
         sample_norm_method='rank', 
         no_plot = True, processes=4, format='png')
    