#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 12:03:20 2022

@author: katarina
***
Modified by: Mauricio Suarez-Duran
"""

today = '0627'

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import json
import codecs
from scipy import stats
from cch_class import CCH, COMBINE_RESULTS, READ_ADST

DIR_OUT = "results/"
#
# Reading ADST files as input
#
fileList = sys.argv[2:]
if not len(fileList):
    print("\n=============================")
    print("Please add path to ADST files\n")
    exit()
'''
q_and_cq_histos = READ_ADST()
#
# Reading ADST files via PYROOT/adst.py 
# and putting histos into dictionaries 
#
cq_dict, q_dict = q_and_cq_histos.load_data(fileList)
dictCQ, dictQ = q_and_cq_histos.load_data(fileList)

with open("jsonFiles/cq_histos_dict_"+sys.argv[1]+".json", "w") as fp:
    json.dump(dictCQ, fp, indent = 4, sort_keys = True)

with open("jsonFiles/q_histos_dict_"+sys.argv[1]+".json", "w") as fp:
    json.dump(dictQ, fp, indent = 4, sort_keys = True)
#
exit(1)
'''
# Opening JSON file
# 
cq_dict = {}
q_dict = {} 
# Reading for February
#
'''
daysFeb = ['09', '10']
for i in daysFeb:
    file_cq_dict = open('jsonFiles/cq_histos_dict_ADST2022_02_'+i+'.json')
    file_q_dict = open('jsonFiles/q_histos_dict_ADST2022_02_'+i+'.json')
    #
    # returns JSON object as a dictionary
    tmp_cq_dict = json.load(file_cq_dict)
    tmp_q_dict = json.load(file_q_dict)
    #
    # Updating dictionaries
    cq_dict.update(tmp_cq_dict)
    q_dict.update(tmp_q_dict)
print("MSD, data from February charged")
'''
# Reading for July
#
'''
for i in range(18, 32):
    file_cq_dict = open('jsonFiles/cq_histos_dict_ADST2022_07_'+str(i)+'.json')
    file_q_dict = open('jsonFiles/q_histos_dict_ADST2022_07_'+str(i)+'.json')
    # returns JSON object as a dictionary
    #
    tmp_cq_dict = json.load(file_cq_dict)
    tmp_q_dict = json.load(file_q_dict)
    # Updating dictionaries
    #
    cq_dict.update(tmp_cq_dict)
    q_dict.update(tmp_q_dict)
print("MSD, data from July charged")
'''
# Reading for August
#

for i in range(1, 17):
    if i < 10:
        file_cq_dict = open('jsonFiles/cq_histos_dict_ADST2022_08_0'+str(i)+'.json')
        file_q_dict = open('jsonFiles/q_histos_dict_ADST2022_08_0'+str(i)+'.json')
    else:
        file_cq_dict = open('jsonFiles/cq_histos_dict_ADST2022_08_'+str(i)+'.json')
        file_q_dict = open('jsonFiles/q_histos_dict_ADST2022_08_'+str(i)+'.json')
    # returns JSON object as a dictionary
    #
    tmp_cq_dict = json.load(file_cq_dict)
    tmp_q_dict = json.load(file_q_dict)
    # Updating dictionaries
    #
    cq_dict.update(tmp_cq_dict)
    q_dict.update(tmp_q_dict)
print("MSD, data from August charged")

# Reading for October
#

file_cq_dict = open('jsonFiles/cq_histos_dict_ADST2022_10_07.json')
file_q_dict = open('jsonFiles/q_histos_dict_ADST2022_10_07.json')
# returns JSON object as a dictionary
#
tmp_cq_dict = json.load(file_cq_dict)
tmp_q_dict = json.load(file_q_dict)
# Updating dictionaries
#
cq_dict.update(tmp_cq_dict)
q_dict.update(tmp_q_dict)
print("MSD, data from October charged")

#
# Creating one CCH object per pair Q-CQ histo
analysis_obj = [ CCH() for i in range( len(cq_dict.keys()) ) ]
cch_list = []
for obj, key_histos in zip(analysis_obj, cq_dict.keys()):
    obj.load_data( np.array(cq_dict[key_histos]['hist']), \
                   np.array(q_dict[key_histos]['hist']), \
                   cq_dict[key_histos]['cqpk'], \
                   cq_dict[key_histos]['cqpkErr'], \
                   q_dict[key_histos]['qpk'], \
                   q_dict[key_histos]['qpkErr'], \
                   cq_dict[key_histos]['fromHistogram'], \
                   cq_dict[key_histos]['signal'], \
                   q_dict[key_histos]['ldfR'], \
                   q_dict[key_histos]['spDist'], \
                   q_dict[key_histos]['Energy'], \
                   q_dict[key_histos]['EnergyErr'], \
                   q_dict[key_histos]['Zenith'], key_histos )
    cch_list.append( key_histos )

print("MSD, CCH objects created")
#
# Fitting histograms
for obj, cch_name in zip(analysis_obj, cch_list):
    obj.fit_all()
print("MSD, histograms fitted")
#
# MSD, comparing with Offline fitting
'''
success_fit = 0
for obj, cch_name in zip(analysis_obj, cch_list):
    success_fit += obj.comparison_offline()
'''
#
# Calculating deltas
tmp = 1
for obj, cch_name in zip(analysis_obj, cch_list):
    try:
        obj.calculate_all_deltas()
    except:
        tmp = 2
        #print('Error delta calculation: '\
                #+cch_name.split()[0]+"_"+cch_name.split()[1])
print("MSD, Deltas calculated")
#
exit()
