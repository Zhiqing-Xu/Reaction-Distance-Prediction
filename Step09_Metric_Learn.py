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
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
#--------------------------------------------------#
import ast
import copy
import pickle
import scipy.io
import subprocess
import numpy as np
import pandas as pd
from numpy import *
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
#--------------------------------------------------#
from AP_RDKIT_FP import *
from Step07_NetworkToDistance import *
#--------------------------------------------------#
##############################################################################################################
##############################################################################################################
loading_folder = Path("MNX_data/")
saving_folder = Path("MNX_ECFP_savings/")
##############################################################################################################
##############################################################################################################
# all_cmpds     :  list( ["X","X",...]         ) # list
# all_ecfps     :  set ( ["ecfp", "ecfp", ...] ) # set
# all_pairs     :  [{{},{}}, {{},{}}, {{},{}},... ]
# all_info      :  [   [  { fr{}, fr{} }, d  ],   [  { fr{}, fr{} }, d  ],  [  { fr{}, fr{} }, d  ], ....  ]
##############################################################################################################
##############################################################################################################
# Args
# Select ECFP encodings
#-------------------      0        1        2        3          4         5        6     
ECFP_encodings_list = ["ECFP2", "ECFP4", "ECFP6", "JTVAE", "MorganFP", "ECFP8", "ECFPX"]
ECFP_encodings = ECFP_encodings_list[1]
ECFP_type = ECFP_encodings[-1] if ECFP_encodings in ["ECFP2", "ECFP4", "ECFP6"] else "6" # 2, 4, 6


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def list_smiles_to_ecfp_through_dict(smiles_list, all_cmpds_ecfps_dict):
    ecfp_list=[]
    for one_smiles in smiles_list:
        ecfp_list=ecfp_list+all_cmpds_ecfps_dict[one_smiles]
    return ecfp_list

#====================================================================================================#
def parse_one_pair_info(one_pair_info, all_ecfps, all_cmpds_ecfps_dict):
    dimension=len(all_ecfps)
    X1i=[0]*dimension
    X2i=[0]*dimension
    X1i_ecfp_list=list_smiles_to_ecfp_through_dict(list(list(one_pair_info[0])[0]),all_cmpds_ecfps_dict)
    X2i_ecfp_list=list_smiles_to_ecfp_through_dict(list(list(one_pair_info[0])[1]),all_cmpds_ecfps_dict)
    distance=one_pair_info[1]
    for one_ecfp in X1i_ecfp_list:
        X1i[all_ecfps.index(one_ecfp)]=X1i_ecfp_list.count(one_ecfp)
    for one_ecfp in X2i_ecfp_list:
        X2i[all_ecfps.index(one_ecfp)]=X2i_ecfp_list.count(one_ecfp)
    Yi=distance
    return (X1i,X2i,Yi)

#====================================================================================================#
def list_subtract(list_a,list_b):
    list_out=[]
    for i in range(len(list_a)):
        list_out.append(list_a[i]-list_b[i])
    return list_out

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


def Step09_main(loading_folder, saving_folder, ECFP_encodings):
    
    #====================================================================================================#
    pickle_in1=open(saving_folder  / "Step07_paired_cmpds_list","rb")
    paired_smiles_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(saving_folder  / "Step07_all_pairs_list","rb")
    all_pairs_list=pickle.load(pickle_in2)
    pickle_in2.close()
    #====================================================================================================#
    pickle_in1=open(saving_folder  / ("Step08_all_cmpds_"+ECFP_encodings),"rb")
    all_smiles=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(saving_folder  / ("Step08_all_ecfps_"+ECFP_encodings),"rb")
    all_ecfps=pickle.load(pickle_in2)
    pickle_in2.close()
    pickle_in3=open(saving_folder  / ("Step08_all_cmpds_ecfps_dict_"+ECFP_encodings),"rb")
    all_smiles_ecfps_dict=pickle.load(pickle_in3)
    pickle_in3.close()
    #====================================================================================================#
    for one_pair_info in paired_smiles_list:
        if len(one_pair_info[0])!=2:
            print (one_pair_info[0])
            print ("wtf?")
            paired_smiles_list.remove(one_pair_info)
    print ("screened!")
    #====================================================================================================#
    all_ecfps=list(all_ecfps)
    X_Diff=[]
    Y_Distance=[]

    for one_pair_info in tqdm(paired_smiles_list):
        (X1i, X2i, Yi)=parse_one_pair_info(one_pair_info,all_ecfps,all_smiles_ecfps_dict)
        X_Diff.append(list_subtract(X1i, X2i))
        Y_Distance.append(Yi)
    Step09_processed_data_dict = {"X_data": X_Diff, "y_data": Y_Distance}
    #====================================================================================================#

    pickle_out1=open(saving_folder / "Step09_processed_data_"+ ECFP_encodings,"wb")
    pickle.dump(Step09_processed_data_dict, pickle_out1)
    pickle_out1.close()


    print("Step09_main Done!")

    return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

if __name__ == '__main__':
    Step09_main(loading_folder, saving_folder, ECFP_encodings)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
