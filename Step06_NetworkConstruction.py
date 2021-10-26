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
import pickle
import scipy.io
import subprocess
import numpy as np
from tqdm import tqdm
from pathlib import Path
from random import shuffle
##############################################################################################################
##############################################################################################################
loading_folder = Path("MNX_data/")
saving_folder = Path("MNX_ECFP_savings/")
##############################################################################################################
##############################################################################################################

def expand_rxn_tree(one_MNXid, RXN_Network_N_step_dict, RXN_dict, max_len=10 ):
    # A brute-force algorithm is used here (INEFFICIENT). 
    # Computation does not take long for this easy problem
    # Efficient graph search algorithm can be used here if constructing a large network

    #print(one_MNXid)
    RXN_tree_N_step=[[] for i in range(max_len)]
    for i in range(max_len):
        if i==0:
            RXN_tree_N_step[i]=list(set(RXN_dict[one_MNXid]))
        else:
            cmpds_tb_expanded=RXN_tree_N_step[i-1]
            next_lv_cmpds=[]
            duplicates=set([one_MNXid,])
            for j in range(i):
                k=i-j-1
                duplicates=duplicates.union(set(RXN_tree_N_step[k]))
            for one_cmpd in cmpds_tb_expanded:
                next_lv_cmpds=[item for item in set(RXN_dict[one_cmpd]) if item not in duplicates]
            RXN_tree_N_step[i]=next_lv_cmpds
    #print(RXN_tree_N_step)
    RXN_Network_N_step_dict[one_MNXid]=RXN_tree_N_step
    return

##############################################################################################################
##############################################################################################################    
def Step06_main(loading_folder, saving_folder):
    #====================================================================================================#
    # 1. Open all pickles and retrieve compound/reaction info.
    pickle_in1=open(saving_folder / "Step03_screened_mnxid_rxn_list","rb")
    screened_mnxid_rxn_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step03_Used_MNXid_set","rb")
    Used_MNXid_set=pickle.load(pickle_in1)
    pickle_in1.close()

    #====================================================================================================#
    # 2. Format the reactions into a dictionary, 
    # Values are a list of MNXid's that are linked to the key with one step reaction
    RXN_dict=dict([])
    for one_MNXid in Used_MNXid_set:
        prod_list=[]
        for one_pair in screened_mnxid_rxn_list:
            if one_MNXid in one_pair:
                temp_pair=list(one_pair)
                temp_pair.remove(one_MNXid)
                prod_list.append(temp_pair[0])
        RXN_dict[one_MNXid]=prod_list

    #====================================================================================================#
    RXN_Network_N_step_dict=dict([])
    print("len(Used_MNXid_set): ", len(Used_MNXid_set))
    for one_MNXid in tqdm(Used_MNXid_set):
        expand_rxn_tree(one_MNXid, RXN_Network_N_step_dict, RXN_dict)

    #for i in tqdm(range(len(list(Used_MNXid_set)))):
    #    one_MNXid = list(Used_MNXid_set)[i]
    #    expand_rxn_tree(one_MNXid, RXN_Network_N_step_dict, RXN_dict)


    pickle_out1=open(saving_folder / "Step06_RXN_Network_N_step_dict","wb")
    pickle.dump(RXN_Network_N_step_dict, pickle_out1)
    pickle_out1.close()

    print ("Step06 done!")

    return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
if __name__ == '__main__':
    Step06_main(loading_folder, saving_folder)

