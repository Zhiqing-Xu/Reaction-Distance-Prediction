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
from tqdm import tqdm
from pathlib import Path
from copy import deepcopy
from random import shuffle
#--------------------------------------------------#
import statistics
#--------------------------------------------------#
from AP_RDKIT_FP import similarity_score
from AP_RDKIT_FP import generate_fingerprint
from AP_RDKIT_FP import similarity_metric_select
##############################################################################################################
##############################################################################################################
loading_folder = Path("MNX_data/")
saving_folder = Path("MNX_ECFP_savings/")
##############################################################################################################
##############################################################################################################
def remove_cofactors(one_rxn_info):
    # More cofactor id might be updated later
    cofactor_mnx_id=["MNXM01","MNXM1","MNXM2","MNXM3","MNXM4","MNXM5","MNXM6","MNXM7","MNXM8",\
                    "MNXM9","MNXM10","MNXM11","MNXM12","MNXM13","MNXM14","MNXM15","MNXM15"]
    for one_cofactor in cofactor_mnx_id:
        if one_cofactor in one_rxn_info[1]:
            one_rxn_info[0].remove(one_rxn_info[0][one_rxn_info[1].index(one_cofactor)])
            one_rxn_info[1].remove(one_cofactor)
        if one_cofactor in one_rxn_info[3]:
            one_rxn_info[2].remove(one_rxn_info[2][one_rxn_info[3].index(one_cofactor)])
            one_rxn_info[3].remove(one_cofactor)     
    return one_rxn_info
#====================================================================================================#
def remove_trivial_cmpds(one_rxn_info,cmpd_mnxid_smiles_dict):
    # Remove metal ions, small compounds, NA SMILES strings and empty strings from the list
    # The algorithm here is not perfected, probably need to take into account more conditions
    #print one_rxn_info
    rxn_info=[[],[],[],[]]
    for i in range(len(one_rxn_info[1])):
        one_cmpd=one_rxn_info[1][i]
        cmpd_smiles=cmpd_mnxid_smiles_dict[one_cmpd]
        if len(cmpd_smiles)>7 and cmpd_smiles.find(".")==-1:
            rxn_info[0].append(one_rxn_info[0][i])
            rxn_info[1].append(one_rxn_info[1][i])        

    for i in range(len(one_rxn_info[3])):
        one_cmpd=one_rxn_info[3][i]
        cmpd_smiles=cmpd_mnxid_smiles_dict[one_cmpd]
        if len(cmpd_smiles)>7 and cmpd_smiles.find(".")==-1:
            rxn_info[2].append(one_rxn_info[2][i])
            rxn_info[3].append(one_rxn_info[3][i])        
    #print rxn_info
    return rxn_info
##############################################################################################################
##############################################################################################################
def Step03_main(loading_folder, saving_folder):
    #====================================================================================================#
    #1. Inputs
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict_merged","rb")
    cmpd_mnxid_smiles_dict=pickle.load(pickle_in1)
    pickle_in1.close()

    pickle_in1=open(saving_folder / "Step02_mnxid_rxn_list","rb")
    mnxid_rxn_list=pickle.load(pickle_in1)
    pickle_in1.close()
    #====================================================================================================#
    #2. 1st Screening
    screened_mnxid_rxn_list_1=[]
    for one_rxn_info in mnxid_rxn_list:
        one_rxn_info=remove_cofactors(one_rxn_info)
        one_rxn_info=remove_trivial_cmpds(one_rxn_info,cmpd_mnxid_smiles_dict)

        if len(one_rxn_info[0])!=0 and len(one_rxn_info[2])!=0 \
            and len(one_rxn_info[0]) == len(one_rxn_info[1]) \
            and len(one_rxn_info[2]) == len(one_rxn_info[3]):
            screened_mnxid_rxn_list_1.append(one_rxn_info)
    print ("Done 1st Screening")
    del(mnxid_rxn_list)
    mnxid_rxn_list=copy.copy(screened_mnxid_rxn_list_1)
    #====================================================================================================#
    #3. 2nd Screening, get reactions in the form (A -> B) 
    # since we only reconstructs linear pathways here.
    screened_mnxid_rxn_list_2=[]

    parameter_1="ECFP"
    for one_rxn_info in mnxid_rxn_list:
        #print one_rxn_info
        if len(one_rxn_info[1])==1 and len(one_rxn_info[3])==1 \
            and one_rxn_info[1][0]!=one_rxn_info[3][0] \
            and set([one_rxn_info[1][0],one_rxn_info[3][0]]) not in screened_mnxid_rxn_list_2: # Consider simple reactions here (A->B)
            screened_mnxid_rxn_list_2.append(set([one_rxn_info[1][0],one_rxn_info[3][0]]))


    print ("Done 2nd Screening")
    print ("Size of screened MNX reaction list: ", len(screened_mnxid_rxn_list_2))
    #-------------------- (4) --------------------#
    # Get all SMILES needed for further analysis
    # No need to convert all SMILES (through VAE) to latent vectors
    Used_MNXid_set=set([])
    for i in screened_mnxid_rxn_list_2:
        Used_MNXid_set=Used_MNXid_set.union(i)
    #print (Used_MNXid_set)

    Used_SMILES_set=set([])
    for one_mnxid in Used_MNXid_set:
        Used_SMILES_set.add(cmpd_mnxid_smiles_dict[one_mnxid])
    #print (Used_SMILES_set)

    pickle_out1=open(saving_folder / "Step03_screened_mnxid_rxn_list","wb")
    pickle.dump(screened_mnxid_rxn_list_2, pickle_out1)
    pickle_out1.close()

    pickle_out1=open(saving_folder / "Step03_Used_MNXid_set","wb")
    pickle.dump(Used_MNXid_set, pickle_out1)
    pickle_out1.close()

    pickle_out1=open(saving_folder / "Step03_Used_SMILES_set","wb")
    pickle.dump(Used_SMILES_set, pickle_out1)
    pickle_out1.close()

    print ("Step03 Done!")

    return



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def Step03_main_extended(loading_folder, saving_folder):
    #====================================================================================================#
    #1. Inputs
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict_merged","rb")
    cmpd_mnxid_smiles_dict=pickle.load(pickle_in1)
    pickle_in1.close()

    pickle_in1=open(saving_folder / "Step02_mnxid_rxn_list","rb")
    mnxid_rxn_list=pickle.load(pickle_in1)
    pickle_in1.close()
    #====================================================================================================#
    #2. 1st Screening
    screened_mnxid_rxn_list_1=[]
    for one_rxn_info in mnxid_rxn_list:
        one_rxn_info=remove_cofactors(one_rxn_info)
        one_rxn_info=remove_trivial_cmpds(one_rxn_info,cmpd_mnxid_smiles_dict)

        if len(one_rxn_info[0])!=0 and len(one_rxn_info[2])!=0 \
            and len(one_rxn_info[0]) == len(one_rxn_info[1]) \
            and len(one_rxn_info[2]) == len(one_rxn_info[3]):
            screened_mnxid_rxn_list_1.append(one_rxn_info)
    print ("Done 1st Screening")
    del(mnxid_rxn_list)
    mnxid_rxn_list=copy.copy(screened_mnxid_rxn_list_1)
    #====================================================================================================#
    #3. 2nd Screening, get reactions in the form (A -> B) 
    # since we only reconstructs linear pathways here.
    screened_mnxid_rxn_list_2=[]

    count_x=0
    count_y=0
    count_u=0
    count_t=0

    pair_AB=[]
    pair_multi_rct=[]
    pair_multi_prt=[]
    pair_multi_r_p=[]


    parameter_1="ECFP"
    for one_rxn_info in mnxid_rxn_list:
        #print one_rxn_info
        if len(one_rxn_info[1])==1 and len(one_rxn_info[3])==1 \
            and one_rxn_info[1][0]!=one_rxn_info[3][0] \
            and set([one_rxn_info[1][0],one_rxn_info[3][0]]) not in screened_mnxid_rxn_list_2: # Consider simple reactions here (A->B)
            screened_mnxid_rxn_list_2.append(set([one_rxn_info[1][0],one_rxn_info[3][0]]))

            count_x+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_AB.append(         ssscore )

        elif len(one_rxn_info[1])==1 and len(one_rxn_info[3])>1:
            count_y+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_multi_rct.append(  ssscore )

        elif len(one_rxn_info[1])>1 and len(one_rxn_info[3])==1:
            count_u+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_multi_prt.append(  ssscore )

        elif len(one_rxn_info[1])>1 and len(one_rxn_info[3])>1:
            count_t+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_multi_r_p.append(  ssscore )

        else:
            pass

    print (count_x,count_y,count_u,count_t)

    print (pair_AB)
    print (pair_multi_rct)
    print (pair_multi_prt)
    print (pair_multi_r_p)

    print (statistics.mean(pair_AB))
    print (statistics.mean(pair_multi_rct))
    print (statistics.mean(pair_multi_prt))
    print (statistics.mean(pair_multi_r_p))

    print ("Done 2nd Screening")
    print (len(screened_mnxid_rxn_list_2))
    #====================================================================================================#
    #4. Get all SMILES needed for further analysis
    # No need to convert all SMILES (through VAE) to latent vectors
    Used_MNXid_set=set([])
    for i in screened_mnxid_rxn_list_2:
        Used_MNXid_set=Used_MNXid_set.union(i)
    print (Used_MNXid_set)

    Used_SMILES_set=set([])
    for one_mnxid in Used_MNXid_set:
        Used_SMILES_set.add(cmpd_mnxid_smiles_dict[one_mnxid])
    print (Used_SMILES_set)

    pickle_out1=open(saving_folder / "Step03_screened_mnxid_rxn_list","wb")
    pickle.dump(screened_mnxid_rxn_list_2, pickle_out1)
    pickle_out1.close()

    pickle_out1=open(saving_folder / "Step03_Used_MNXid_set","wb")
    pickle.dump(Used_MNXid_set, pickle_out1)
    pickle_out1.close()

    pickle_out1=open(saving_folder / "Step03_Used_SMILES_set","wb")
    pickle.dump(Used_SMILES_set, pickle_out1)
    pickle_out1.close()

    print ("Step03 Done!")

    return 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
if __name__ == '__main__':
    Step03_main(loading_folder, saving_folder)


