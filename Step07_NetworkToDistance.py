#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Microsoft VS header
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
if os.name == 'nt' or platform == 'win32':
    print("Running on Windows")
    if 'ptvsd' in sys.modules:
        print("Running in Visual Studio")
        try:
            os.chdir(os.path.dirname(__file__))
            print('CurrentDir: ', os.getcwd())
        except:
            pass
#--------------------------------------------------#
    else:
        print("Running outside Visual Studio")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            pass
#--------------------------------------------------#
from rdkit import Chem
from rdkit.Chem import AllChem
#--------------------------------------------------#
import ast
import copy
import pickle
import scipy.io
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from random import shuffle
#--------------------------------------------------#
import seaborn as sns
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#--------------------------------------------------#
from scipy import stats
from matplotlib import pylab as pl
##############################################################################################################
##############################################################################################################
loading_folder = Path("MNX_data/")
saving_folder = Path("MNX_ECFP_savings/")
##############################################################################################################
##############################################################################################################
# Step07_paired_cmpds_list  : [ [paired_cmpds_set, distance], [ set([frozenset(), frozenset()]) , d ], 
#                               [ set([fs(), fs()]) , d ], 
#                               [ set([ fs , fs ]) , d ], ..., [], [] , ... ]
# Step07_all_pairs_list     : [ paired_cmpds_set, 
#                               set( [  frozenset( cmpd_x, ), frozenset(cmpd_x, cmpd_x, cmpd_x, ...)  ] ), 
#                               set( [fs(), fs()] ), 
#                               set( [ fs , fs ] ), ..., set([]), ...  ]

def update_paired_cmpds_list(distance, paired_cmpds_set, paired_cmpds_list, all_pairs_list):
    # Used in parse_one_pathway(one_pathway, paired_cmpds_list, all_pairs_list)
    if paired_cmpds_set not in all_pairs_list:
        all_pairs_list.append(paired_cmpds_set)
        paired_cmpds_list.append([paired_cmpds_set,distance])
    else: 
        for i in range(len(paired_cmpds_list)):
            if paired_cmpds_list[i][0]==paired_cmpds_set and paired_cmpds_list[i][1]>distance: # means we have a update in distance 
                #print ("update for shorter distance found")
                paired_cmpds_list[i][1]=distance
    return

##############################################################################################################
##############################################################################################################

def Step07_main(loading_folder, saving_folder): # Time-consuming
    #====================================================================================================#
    # Open the Step01 temp file and retrieve useful info.
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict_merged","rb")
    cmpd_mnxid_smiles_dict=pickle.load(pickle_in1)
    pickle_in1.close()
    
    #====================================================================================================#
    # Open the Step06 temp file and retrieve useful info.
    pickle_in1=open(saving_folder / "Step06_RXN_Network_N_step_dict","rb")
    RXN_Network_N_step_dict=pickle.load(pickle_in1)
    pickle_in1.close()

    #====================================================================================================#
    # Again, use brute-force here since it wont take too long for this amount of data
    print ("Start working on Cmpd Pairs & Distance", len(RXN_Network_N_step_dict))
    len_pwy = len(list(RXN_Network_N_step_dict.values())[0])
    print ("len_pwy (steps of expansion): ", len_pwy) # steps of expansion

    #====================================================================================================#
    paired_cmpds_list=[]
    all_pairs_list=[]
    for one_cmpd in tqdm(RXN_Network_N_step_dict.keys()): # one_cmpd here is MNX_id

        for i in range(len_pwy):
            for another_cmpd in RXN_Network_N_step_dict[one_cmpd][i]:
                distance = i+1
                paired_cmpds=[(cmpd_mnxid_smiles_dict[one_cmpd],),(cmpd_mnxid_smiles_dict[another_cmpd],)]
                paired_cmpds_set=set(frozenset(i) for i in paired_cmpds) # set?????
                update_paired_cmpds_list(distance, paired_cmpds_set, paired_cmpds_list, all_pairs_list)

    save_pickle_paired_cmpds_list = saving_folder / "Step07_paired_cmpds_list"
    save_pickle_all_pairs_list = saving_folder / "Step07_all_pairs_list"
    
    pickle_out1=open(save_pickle_paired_cmpds_list,"wb")
    pickle.dump(paired_cmpds_list, pickle_out1)
    pickle_out1.close()

    pickle_out2=open(save_pickle_all_pairs_list,"wb")
    pickle.dump(all_pairs_list, pickle_out2)
    pickle_out2.close()

    print ("Step07 Done!")

    return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
if __name__ == '__main__':
    Step07_main(loading_folder, saving_folder) # Time-consuming


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

def music(num): # Alarm when the computating is done!
    import platform
    try:
        import winsound
        beep = winsound.Beep
    except ImportError:
        def beep(f, d):
            s = 8000
            hp = int(s/f/2)
            b = chr(255)*hp+chr(0)*hp
            b *= int(d*f)
            a = file('/dev/audio', 'wb')
            a.write(b)
            a.close()    

    c = [(880, 1000),
         (587, 1000),
         (698, 500), 
         (880, 500), 
         (587, 1000),
         (698, 500),
         (880, 250),
         (1046, 250), 
         (988, 500),
         (784, 500),
         (699, 230),
         (784, 250), 
         (880, 500),
         (587, 500),
         (523, 250),
         (659, 250),
         (587, 750)]

    s = c * num
    for f, d in s:
        beep(f, d)
    return

##############################################################################################################
##############################################################################################################