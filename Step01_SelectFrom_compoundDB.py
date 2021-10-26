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
#########################################################################################################
#########################################################################################################
loading_folder = Path("MNX_data/")
saving_folder = Path("MNX_ECFP_savings/")
#########################################################################################################
#########################################################################################################
def readlines_chem_prop(file, mark="", skip_header=False):
    count_x=0
    cmpd_mnxid_smiles_dict=dict([])
    # skip lines not containing information
    # There are 388 lines of header (trivial information about the DB)
    if skip_header==True:
        for i in range(388): 
            line = file.readline()

    # Process 10^6 at a time in case there are problematic lines in the data file.
    # It is easy to combine the dictionaries afterwards.
    for i in tqdm(range(100000)): 
        count_x+=1
        line = file.readline()
        if line != "" and line != "\n" :
            cmpd_info_list=line.split("\t")
            #print (cmpd_info_list[0], count_x+1)
            cmpd_mnxid_smiles_dict[cmpd_info_list[0]]=cmpd_info_list[6]
        else:
            return cmpd_mnxid_smiles_dict
    return cmpd_mnxid_smiles_dict
#########################################################################################################
#########################################################################################################
def Step01_main_1(loading_folder, saving_folder):
    # Process 10^6 at a time in case there are problematic lines in the data file.
    # Write the compound info into 7 dictionaries, each info of with up to 10^6 compounds
    # It is easy to combine the dictionaries with a few lines of code.
    file_MNX_cmpd_address = loading_folder / "chem_prop.tsv"
    file_MNX_cmpd = open(file_MNX_cmpd_address)
    cmpd_mnxid_smiles_dict1=readlines_chem_prop(file_MNX_cmpd,"",True)
    cmpd_mnxid_smiles_dict2=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict3=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict4=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict5=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict6=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict7=readlines_chem_prop(file_MNX_cmpd,"",False)

    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict1","wb")
    pickle.dump(cmpd_mnxid_smiles_dict1, pickle_out1)
    pickle_out1.close()
    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict2","wb")
    pickle.dump(cmpd_mnxid_smiles_dict2, pickle_out1)
    pickle_out1.close()
    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict3","wb")
    pickle.dump(cmpd_mnxid_smiles_dict3, pickle_out1)
    pickle_out1.close()
    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict4","wb")
    pickle.dump(cmpd_mnxid_smiles_dict4, pickle_out1)
    pickle_out1.close()
    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict5","wb")
    pickle.dump(cmpd_mnxid_smiles_dict5, pickle_out1)
    pickle_out1.close()
    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict6","wb")
    pickle.dump(cmpd_mnxid_smiles_dict6, pickle_out1)
    pickle_out1.close()
    pickle_out1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict7","wb")
    pickle.dump(cmpd_mnxid_smiles_dict7, pickle_out1)
    pickle_out1.close()

#########################################################################################################
#########################################################################################################
def Step01_main_2(loading_folder, saving_folder):
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict1","rb")
    cmpd_mnxid_smiles_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict2","rb")
    cmpd_mnxid_smiles_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict3","rb")
    cmpd_mnxid_smiles_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict4","rb")
    cmpd_mnxid_smiles_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict5","rb")
    cmpd_mnxid_smiles_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict6","rb")
    cmpd_mnxid_smiles_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open(saving_folder / "Step01_cmpd_mnxid_smiles_dict7","rb")
    cmpd_mnxid_smiles_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    cmpd_mnxid_smiles_dict=cmpd_mnxid_smiles_dict1.copy()
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict2)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict3)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict4)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict5)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict6)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict7)

    # Add some keys to the dict to avoid errors
    cmpd_mnxid_smiles_dict["BIOMASS"]=""
    cmpd_mnxid_smiles_dict["MNXM0"]=""
    print ("done merge")

    pickle_out1 = open(saving_folder / "Step01_cmpd_mnxid_smiles_dict_merged","wb")
    pickle.dump(cmpd_mnxid_smiles_dict, pickle_out1)
    pickle_out1.close()

    print("Step01 Done!")

#########################################################################################################
#########################################################################################################
if __name__ == '__main__':
    Step01_main_1(loading_folder, saving_folder)
    Step01_main_2(loading_folder, saving_folder)
