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

q_and_cq_histos = READ_ADST()
#
# Reading ADST files via PYROOT/adst.py 
# and putting histos into dictionaries 
#
cq_dict, q_dict = q_and_cq_histos.load_data(fileList)
dictCQ, dictQ = q_and_cq_histos.load_data(fileList)
exit(1)
with open("jsonFiles2/cq_histos_dict_"+sys.argv[1]+".json", "w") as fp:
    json.dump(dictCQ, fp, indent = 4, sort_keys = True)

with open("jsonFiles2/q_histos_dict_"+sys.argv[1]+".json", "w") as fp:
    json.dump(dictQ, fp, indent = 4, sort_keys = True)

exit()

# Opening JSON file
# 
cq_dict = {}
q_dict = {} 
# Reading for February
#

daysFeb = ['09', '10']
for i in daysFeb:
    file_cq_dict = open('jsonFiles2/cq_histos_dict_ADST2022_02_'+i+'.json')
    file_q_dict = open('jsonFiles2/q_histos_dict_ADST2022_02_'+i+'.json')
    # returns JSON object as a dictionary
    #
    tmp_cq_dict = json.load(file_cq_dict)
    tmp_q_dict = json.load(file_q_dict)
    # Updating dictionaries
    #
    for i in tmp_cq_dict.keys():
        print(i.split()[1])
    cq_dict.update(tmp_cq_dict)
    q_dict.update(tmp_q_dict)
print("MSD, data from February charged")

# Reading for July
#
'''
for i in range(18, 32):
    file_cq_dict = open('jsonFiles2/cq_histos_dict_ADST2022_07_'+str(i)+'.json')
    file_q_dict = open('jsonFiles2/q_histos_dict_ADST2022_07_'+str(i)+'.json')
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
'''
for i in range(1, 17):
    if i < 10:
        file_cq_dict = open('jsonFiles2/cq_histos_dict_ADST2022_08_0'+str(i)+'.json')
        file_q_dict = open('jsonFiles2/q_histos_dict_ADST2022_08_0'+str(i)+'.json')
    else:
        file_cq_dict = open('jsonFiles2/cq_histos_dict_ADST2022_08_'+str(i)+'.json')
        file_q_dict = open('jsonFiles2/q_histos_dict_ADST2022_08_'+str(i)+'.json')
    # returns JSON object as a dictionary
    #
    tmp_cq_dict = json.load(file_cq_dict)
    tmp_q_dict = json.load(file_q_dict)
    # Updating dictionaries
    #
    cq_dict.update(tmp_cq_dict)
    q_dict.update(tmp_q_dict)
print("MSD, data from August charged")
'''
# Reading for October
#
'''
file_cq_dict = open('jsonFiles2/cq_histos_dict_ADST2022_10_07.json')
file_q_dict = open('jsonFiles2/q_histos_dict_ADST2022_10_07.json')
# returns JSON object as a dictionary
#
tmp_cq_dict = json.load(file_cq_dict)
tmp_q_dict = json.load(file_q_dict)
# Updating dictionaries
#
cq_dict.update(tmp_cq_dict)
q_dict.update(tmp_q_dict)
print("MSD, data from October charged")
'''
#
# Creating one CCH object per pair Q-CQ histo
analysis_obj = [ CCH() for i in range( len(cq_dict.keys()) ) ]
cch_list = []
for obj, key_histos in zip(analysis_obj, cq_dict.keys()):
    #print(key_histos)
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
for obj, cch_name in zip(analysis_obj, cch_list):
    try:
        obj.calculate_all_deltas()
    except:
        print('Error delta calculation: '\
                +cch_name.split()[0]+"_"+cch_name.split()[1])
print("MSD, Deltas calculated")
#
'''
how_many=[0, 0, 0]
for obj, cch_name in zip(analysis_obj, cch_list):
    sd = list(obj.uberdict.keys())[0]
    for pmt in range(1, 4):
        try:
            pvalcch = obj.uberdict[sd][pmt]['p_val_cch']
            pvaloch = obj.uberdict[sd][pmt]['p_va   l_och']
            pvalval = obj.uberdict[sd][pmt]['p_val_valley']
            Qerrcch = obj.uberdict[sd][pmt]['qpk_cch'][1]
            Qerroch = obj.uberdict[sd][pmt]['qpk_och'][1]
        except:
            print('Some pval(s) missing in: '\
                    +cch_name.split()[0]+"_"+cch_name.split()[1]\
                    +', '+str(pmt)) 
            continue
        how_many[pmt-1] += 1
        if Qerroch>50 or Qerrcch>50:
            try:
                obj.plot_fit(sd, pmt, save=[ciel])
            except:
                print('Error plotting fit in: '\
                        +cch_name.split()[0]+"_"+cch_name.split()[1]+', '\
                      +str(pmt))
'''
exit(1)
#
res = COMBINE_RESULTS(analysis_obj, DIR_OUT)
res.average_objects()
res.calculate_deltas_mean()
#res.plot_EventsPerStation()
print("MSD, Combined_results done")
'''
res.set_delta_mean({1: [0.025111184088437355,0.0007504875572834235],\
                    2: [0.0662085184325302,0.001261373355716716],\
                        3: [0.025111184088437355,0.0007504875572834235]})
'''
#
plt.rcParams.update({'font.size': 24})
res.plot_Delta_vs_vh(save=True, errs=True, exclude={1:[], 2:['871'], 3:[]}) #, exclude={1:['871'],2:['871'],3:['871']})
res.plot_delta_bands(save=True)
res.plot_delta_histo(save=True, fit_distrib=True, only_peak_fit=True)
res.plot_deltaMix_histo(save=True)
res.plot_logP(save=True,log_counts=False, log=False)
res.plot_signalsRelDiffPmt(save=True)
res.plot_signalsResolution(save=True)
print("MSD, plotting done")
#
#res.vh_check()
#res.deltas_check()
'''
#
# How many Deltas per station?
# make data matrix sd, pmt, time, QpkCCH, QpkCCh, errQpkOCH, QpkOCH err
matrixx = []
for obj, cch_name in zip(analysis_obj,cch_list):
    sd=list(obj.uberdict.keys())[0]
    for pmt in range(1,4):
        try:
            matrixx.append([sd, pmt, obj.uberdict[sd]['utc_times'][0],\
                *obj.uberdict[sd][pmt]['qpk_cch'], *obj.uberdict[sd][pmt]['qpk_och']])
        except:
            continue

matrixxx = []
for key in res.results.keys():
    for pmt in range(1,4):
        try:
            matrixxx.append([key, pmt, *res.results[key][pmt]['Delta'],\
                             *res.results[key][pmt]['v/h']])
        except:
            continue

matrixx=np.array(matrixx)
matrixxx=np.array(matrixxx)

# storing results  
np.save('fit_results.npy', matrixx)
np.save('delta_results.npy', matrixxx)
  
import pickle
with open('results_dict.pickle', 'wb') as handle:
    pickle.dump(res.results, handle, protocol=pickle.HIGHEST_PROTOCOL)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
ciel = PdfPages(DIR_OUT+'plots/Delta_band_highest_853.pdf')
sd, pmt = 853, 2
for i, obj in enumerate(analysis_obj):
    whichsd=list(obj.uberdict.keys())[0]
    # gpst=np.round(obj.uberdict[whichsd]['gps_times'][0])
    if sd==whichsd and obj.uberdict[sd][pmt]['Delta'][0]>0.2:
        obj.plot_fit(sd, pmt, save=[ciel])
        break
ciel.close()
'''
exit()
