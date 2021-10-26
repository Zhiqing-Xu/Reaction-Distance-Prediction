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
    print("running on windows")
    if 'ptvsd' in sys.modules:
        print("running in visual studio")
        try:
            os.chdir(os.path.dirname(__file__))
            print('currentdir: ', os.getcwd())
        except:
            pass
#--------------------------------------------------#
    else:
        print("running outside visual studio")
        try:
            if not 'workbookdir' in globals():
                workbookdir = os.getcwd()
                print('workbookdir: ' + workbookdir)
                os.chdir(workbookdir)
        except:
            pass
#--------------------------------------------------#
if os.name != 'nt' and platform != 'win32':
    print("Not Running on Windows")
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
import gzip
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


def similarity_metric_select(fp_a,fp_b,parameter_1,parameter=2):
    if (parameter_1=="top"):
        similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
    elif (parameter_1=="MACCS"):
        similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
    elif (parameter_1=="atom_pairs"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    elif (parameter_1=="vec_pairs"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    elif (parameter_1=="torsions"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    elif (parameter_1=="FCFP"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    else: # ECFP
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    return similarity

#============================================================================================================================#
def generate_fingerprint(smiles_a, parameter_1, parameter_2=2):
    try: 
        cmpd_a=Chem.MolFromSmiles(str(smiles_a))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
    except Exception:
        #print ("Rdkit ERROR: generate fingerprint, ", smiles_a)
        cmpd_a=Chem.MolFromSmiles(str('O'))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
    return fp_a

#============================================================================================================================#
def similarity_score(smiles_a, smiles_b, parameter_1="ECFP", parameter_2=2): # Return the similarity of two compounds
    try:
        # parameter_1 is similarity metric selected
        cmpd_a=Chem.MolFromSmiles(str(smiles_a))
        cmpd_b=Chem.MolFromSmiles(str(smiles_b))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
            fp_b=FingerprintMols.FingerprintMol(cmpd_b)  
            similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
            fp_b=MACCSkeys.GenMACCSKeys(cmpd_b)
            similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
            fp_b=Pairs.GetAtomPairFingerprint(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
            fp_b=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
            fp_b=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
            fp_b=AllChem.GetMorganFingerprint(cmpd_b,parameter_2,useFeatures=True)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
            fp_b=AllChem.GetMorganFingerprint(cmpd_b,parameter_2)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    except Exception:
        if smiles_a.find("CoA")==-1 and smiles_b.find("CoA")==-1:
            similarity=0
            #print ("Rdkit ERROR: similarity score, ", smiles_a, smiles_b)
        else:
            similarity=1
    return similarity
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

def maxsimscore(list_a,list_b,fptype,parameter_2=2):
    score_list=[]
    for smiles_a in list_a:
        for smiles_b in list_b:
            score_list.append(similarity_score(smiles_a, smiles_b, fptype, parameter_2=2))
    return max(score_list)

#============================================================================================================================#
def replace_n(str1, n, str2):
    letters = (
    str2 if i == n else char
        for i, char in enumerate(str1)
    )
    return ''.join(letters)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# ECFP from CDK java file
def CDK_ECFP(smiles_str,ecfp_type,iteration_number):    
    # Use java file CDKImpl class to get ECFP from cmd line
    query_str1='java -cp .;cdk-2.2.jar CDKImpl ' + smiles_str + ' ' + ecfp_type + ' ' + str(iteration_number)
    query_result = subprocess.check_output(query_str1, shell=True)
    query_result = query_result.decode("gb2312")
    query_result=query_result.replace('[','')
    query_result=query_result.replace(']','')
    query_result=query_result.replace(' ','')
    query_result=query_result.replace('\n','')
    query_result=query_result.replace('\r','')
    if query_result!="":
        if query_result[-1]==',':
            query_result=query_result[0:-1]
        list_of_ecfp=query_result.split(",")
    else:
        list_of_ecfp=[]
    return list_of_ecfp 
#====================================================================================================#
def get_full_ecfp(smiles_str,ecfp_type,iteration_number):   
    # ECFP4 + itr2 or ECFP2 + itr1
    full_ecfp_list=[]
    for i in range(iteration_number+1):
        full_ecfp_list=full_ecfp_list+CDK_ECFP(smiles_str,ecfp_type,i)
    return full_ecfp_list
#====================================================================================================#
def generate_all_ECFPs(list_smiles,ecfp_type="ECFP6",iteration_number=3):
# return a list of ECFPs of all depth for a list of compounds (UNIQUE!!!)
    all_ecfps=set([])
    for smiles_a in list_smiles:
        discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
        #print(smiles_a)
        all_ecfps=all_ecfps.union(set(discriptors))
    return all_ecfps
#====================================================================================================#
def update_all_ECFPs(new_list_smiles, list_smiles, all_ecfps,ecfp_type="ECFP6", iteration_number=3):
    for smiles_a in new_list_smiles:
        if smiles_a not in list_smiles:
            discriptors = get_full_ecfp(smiles_a, ecfp_type, iteration_number)
            #print(smiles_a)
            all_ecfps = all_ecfps.union(set(discriptors))
    return all_ecfps
#====================================================================================================#
def generate_all_smiles_ecfps_dict(list_smiles, ecfp_type="ECFP6", iteration_number=3):
    all_smiles_ecfps_dict=dict([])
    for smiles_a in list_smiles:
        #print(smiles_a)
        all_smiles_ecfps_dict[smiles_a]=get_full_ecfp(smiles_a,ecfp_type,iteration_number)
    return all_smiles_ecfps_dict
#====================================================================================================#
def update_all_smiles_ecfps_dict(new_list_smiles, list_smiles, all_smiles_ecfps_dict, ecfp_type="ECFP6", iteration_number=3):
    for smiles_a in new_list_smiles:
        if smiles_a not in list_smiles:
            #print(smiles_a)
            all_smiles_ecfps_dict[smiles_a] = get_full_ecfp(smiles_a, ecfp_type, iteration_number)
    return all_smiles_ecfps_dict
#====================================================================================================#
def generate_all_smiles_ecfps_list_dict(list_smiles, ecfp_type="ECFP6", iteration_number=3):
    all_ecfps=set([])
    all_smiles_ecfps_dict=dict([])
    for smiles_a in tqdm(list_smiles):
        discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
        #print(smiles_a)
        all_smiles_ecfps_dict[smiles_a]=discriptors
        all_ecfps=all_ecfps.union(set(discriptors))
    return list(all_ecfps),all_smiles_ecfps_dict ##### ????? why list here?
#====================================================================================================#
def update_all_ecfps_and_all_smiles_ecfps_dict(new_list_smiles, list_smiles, all_ecfps, all_smiles_ecfps_dict, ecfp_type="ECFP6", iteration_number=3):
    #count_x=0
    for smiles_a in new_list_smiles:
        if smiles_a not in list_smiles:
            discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
            #print count_x
            #count_x+=1
            all_smiles_ecfps_dict[smiles_a]=discriptors
            all_ecfps=all_ecfps.union(set(discriptors))
    return all_ecfps, all_smiles_ecfps_dict ##### ????? why NOT list here?
###################################################################################################################
###################################################################################################################
# Encode Substrates.
def list_smiles_to_ecfp_through_dict(smiles_list, all_smiles_ecfps_dict):
    ecfp_list=[]
    for one_smiles in smiles_list:
        ecfp_list=ecfp_list + all_smiles_ecfps_dict[one_smiles]
    return ecfp_list
#====================================================================================================#
def smiles_to_ECFP_vec( smiles_x, all_ecfps, all_smiles_ecfps_dict):
    dimension=len(all_ecfps)
    Xi=[0]*dimension
    Xi_ecfp_list=list_smiles_to_ecfp_through_dict( [smiles_x, ] ,all_smiles_ecfps_dict)
    for one_ecfp in Xi_ecfp_list:
        Xi[all_ecfps.index(one_ecfp)]=Xi_ecfp_list.count(one_ecfp)
    return np.array(Xi)
#====================================================================================================#
def Get_ECFPs_encoding(X_subs_representations, all_ecfps, all_smiles_ecfps_dict):
    X_subs_encodings=[]
    for one_smiles in X_subs_representations:
        one_subs_encoding = smiles_to_ECFP_vec(one_smiles, all_ecfps, all_smiles_ecfps_dict) # substrate_encoding
        X_subs_encodings.append(one_subs_encoding)
    return X_subs_encodings
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#








