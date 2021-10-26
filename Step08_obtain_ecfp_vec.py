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
from sklearn.metrics import confusion_matrix, classification_report
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
##############################################################################################################
##############################################################################################################
loading_folder = Path("MNX_data/")
saving_folder = Path("MNX_ECFP_savings/")
##############################################################################################################
##############################################################################################################
# all_smiles     :  list( ["X","X",...]         ) # list
# all_ecfps     :  set ( ["ecfp", "ecfp", ...] ) # set
# all_pairs     :  [{{},{}}, {{},{}}, {{},{}},... ]
# all_info      :  [   [  { fr{}, fr{} }, d  ],   [  { fr{}, fr{} }, d  ],  [  { fr{}, fr{} }, d  ], ....  ]
##############################################################################################################
##############################################################################################################
# Args
# Select ECFP encodings
#-------------------      0        1        2        3          4         5        6     
ECFP_encodings_list = ["ECFP2", "ECFP4", "ECFP6", "JTVAE", "MorganFP", "ECFP8", "ECFPX"]
ECFP_encodings = ECFP_encodings_list[0]
ECFP_type = ECFP_encodings[-1] if ECFP_encodings in ["ECFP2", "ECFP4", "ECFP6"] else "6" # 2, 4, 6


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def Step08_open_Step07_pickles(Step07_paired_cmpds_list_file_address, Step07_all_pairs_list_file_address):
    # Step07_paired_cmpds_list  : [ [paired_cmpds_set, distance], [ set([frozenset(), frozenset()]) , d ], 
    #                               [ set([fs(), fs()]) , d ], 
    #                               [ set([ fs , fs ]) , d ], ..., [], [] , ... ]
    # Step07_all_pairs_list     : [ paired_cmpds_set, 
    #                               set( [  frozenset( cmpd_x, ), frozenset(cmpd_x, cmpd_x, cmpd_x, ...)  ] ), 
    #                               set( [fs(), fs()] ), 
    #                               set( [ fs , fs ] ), ..., set([]), ...  ]
 
    print ("loading data!")
    pickle_in1=open(Step07_paired_cmpds_list_file_address,"rb")
    paired_smiles_list=pickle.load(pickle_in1)
    pickle_in1.close()

    pickle_in2=open(Step07_all_pairs_list_file_address,"rb")
    all_pairs_list=pickle.load(pickle_in2)
    pickle_in2.close()
    print ("pickle data loaded!")

    return paired_smiles_list, all_pairs_list

##############################################################################################################
##############################################################################################################
def get_all_smiles_from_all_pairs_list(all_pairs_list):
    # all_pairs_list     : [ paired_cmpds_set, 
    #                        set( [  frozenset( cmpd_x, ), frozenset(cmpd_x, cmpd_x, cmpd_x, ...)  ] ), 
    #                        set( [fs(), fs()] ), 
    #                        set( [ fs , fs ] ), ..., set([]), ...  ]

    all_smiles=[]
    for one_pair_set in all_pairs_list:
        for one_frozen_set in one_pair_set:
            for one_cmpd_x in one_frozen_set:
                if one_cmpd_x not in all_smiles:
                    all_smiles.append(one_cmpd_x)
    return all_smiles


##############################################################################################################
##############################################################################################################
def Initialize_all_smiles_all_ecfps(loading_folder, saving_folder, ECFP_encodings):

    #====================================================================================================#
    Step07_paired_cmpds_list_file_address = saving_folder  / "Step07_paired_cmpds_list" # "_0" files are old files
    Step07_all_pairs_list_file_address = saving_folder / "Step07_all_pairs_list" # "_0" files are old files
    paired_smiles_list, all_pairs_list = Step08_open_Step07_pickles(Step07_paired_cmpds_list_file_address, Step07_all_pairs_list_file_address)
    
    all_smiles = get_all_smiles_from_all_pairs_list(all_pairs_list)
    print (len(all_smiles))
    print (len(all_pairs_list))
    (all_ecfps, all_smiles_ecfps_dict) = generate_all_smiles_ecfps_list_dict(all_smiles, 
                                                                             ecfp_type = ECFP_encodings, 
                                                                             iteration_number = round(int(ECFP_encodings[-1])/2.) )

    #====================================================================================================#
    pickle_out1=open(saving_folder  / ("Step08_all_cmpds_" + ECFP_encodings),"wb")
    pickle.dump(all_smiles, pickle_out1)
    pickle_out1.close()
    pickle_out2=open(saving_folder  / ("Step08_all_ecfps_" + ECFP_encodings),"wb")
    pickle.dump(all_ecfps, pickle_out2)
    pickle_out2.close()
    pickle_out3=open(saving_folder  / ("Step08_all_cmpds_ecfps_dict_" + ECFP_encodings),"wb")
    pickle.dump(all_smiles_ecfps_dict, pickle_out3)
    pickle_out3.close()

    print ("Done Initialize_all_smiles_all_ecfps! ")

    return

##############################################################################################################
##############################################################################################################
def Step08_show_similarities_distn(loading_folder, saving_folder, ECFP_encodings):
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
    print(len(paired_smiles_list))
    shuffle(paired_smiles_list)
    #====================================================================================================#
    dis_list=[]
    sim_list=[]
    sim_inferred_list=[]
    dis_inferred_list=[]
    count_x=0
    for one_pair_info in tqdm(paired_smiles_list): #[   [ { fr{}, fr{} },d ],   [ { fr{}, fr{} },d ],  [{{},{}},d], ....  ]
        count_x+=1
        #print count_x
        if count_x>=10002:
            break
        if len(list(list(one_pair_info[0])))==2:
            score=maxsimscore(list(list(one_pair_info[0])[0]),list(list(one_pair_info[0])[1]),"MACCS")
            
            distance = one_pair_info[1]
            dis_list.append(distance)

            sim_inferred = 1/(distance+1.)
            sim_inferred_list.append(sim_inferred)

            similarity=score
            sim_list.append(similarity)

            dis_inferred = round(1./(similarity+0.1) - 1)
            dis_inferred_list.append(dis_inferred)

    sim_dist_1=[]
    sim_dist_2=[]
    sim_dist_3=[]
    sim_dist_4=[]
    sim_dist_5=[]
    sim_dist_6=[]
    sim_dist_7=[]
    sim_dist_8=[]
    sim_dist_9=[]
    sim_dist_10=[]
    #====================================================================================================#
    for i in tqdm(range(5000)):
        distance = dis_list[i]
        prediction = sim_list[i]
        if distance==1:
            sim_dist_1.append(prediction)
        if distance==2:
            sim_dist_2.append(prediction)
        if distance==3:
            sim_dist_3.append(prediction)
        if distance==4:
            sim_dist_4.append(prediction)
        if distance==5:
            sim_dist_5.append(prediction)
        if distance==6:
            sim_dist_6.append(prediction)
        if distance==7:
            sim_dist_7.append(prediction)
        if distance==8:
            sim_dist_8.append(prediction)
        if distance==9:
            sim_dist_9.append(prediction)
        if distance==10:
            sim_dist_10.append(prediction)
    #====================================================================================================#
    x=[]
    x.append(np.array(sim_dist_1))
    x.append(np.array(sim_dist_2))
    x.append(np.array(sim_dist_3))
    x.append(np.array(sim_dist_4))
    x.append(np.array(sim_dist_5))
    x.append(np.array(sim_dist_6))
    x.append(np.array(sim_dist_7))
    x.append(np.array(sim_dist_8))  
    x.append(np.array(sim_dist_9)) 
    x.append(np.array(sim_dist_10)) 

    #====================================================================================================#
    col_list=["red","orange","goldenrod","darkgreen","cyan","darkblue","purple","grey", "saddlebrown", "pink"]
    #====================================================================================================#
    def plot_type_1():
        plt.figure()
        for i in [0,1,2,3,4,5,6,7,8,9]:
            plt.subplot(10, 1, i+1)
            #sns.kdeplot(x[i], bw = 0.01 , color="darkred")
            #sns.distplot(x[i], hist = 1, bins=np.arange(0.5,10.5,0.333), rug=False, kde=True, kde_kws={'bw':0.1},hist_kws=dict(ec="k"), color=col_list[i])
            sns.distplot(x[i], hist = 1, bins=40, rug=False, kde=True, kde_kws={'bw':0.05},hist_kws=dict(ec="k"), color=col_list[i])
            plt.xlim((0,1))
        plt.show()
        return
    #====================================================================================================#
    def plot_type_2():
        plt.figure()
        for i in [0,1,2,3,4,5,6,7,8,9]:
            plt.subplot(10, 1, i+1)
            #sns.kdeplot(x[i], bw = 0.01 , color="darkred")
            #sns.distplot(x[i], hist = 1, bins=np.arange(0.5,10.5,0.333), rug=False, kde=True, kde_kws={'bw':0.1},hist_kws=dict(ec="k"), color=col_list[i])
            sns.distplot(x[i], hist = 1, bins=40, rug=False, kde=True, kde_kws={'bw':0.02},hist_kws=dict(ec="k"), color=col_list[i])
            plt.xlim((0,1))
        plt.show()
        return
    #====================================================================================================#
    plot_type_1()
    plot_type_2()
    #====================================================================================================#
    print (corrcoef(sim_list,sim_inferred_list)[1,0])
    print (corrcoef(dis_list,dis_inferred_list)[1,0])
    #====================================================================================================#
    cm=confusion_matrix(dis_list,dis_inferred_list)
    confusion_matrix_df = pd.DataFrame(cm)

    plt.figure()
    sns.heatmap(confusion_matrix_df, cmap="Reds", center=550, annot=True, fmt="d")
    #sns.heatmap(confusion_matrix_df, center=250, annot=True, fmt="d")
    plt.show()

##############################################################################################################
##############################################################################################################
def test_difference_bwt_old_and_new(loading_folder, saving_folder):

    Step07_paired_cmpds_list_file_address = saving_folder  / "Step07_paired_cmpds_list_0" # "_0" files are old files
    Step07_all_pairs_list_file_address = saving_folder / "Step07_all_pairs_list_0" # "_0" files are old files
    open_lists = Step08_open_Step07_pickles(Step07_paired_cmpds_list_file_address ,Step07_all_pairs_list_file_address)
    paired_smiles_list_0=open_lists[0]
    all_pairs_list_0=open_lists[1]

    Step07_paired_cmpds_list_file_address = saving_folder  / "Step07_paired_cmpds_list" # "_0" files are old files
    Step07_all_pairs_list_file_address = saving_folder / "Step07_all_pairs_list" # "_0" files are old files
    open_lists = Step08_open_Step07_pickles(Step07_paired_cmpds_list_file_address, Step07_all_pairs_list_file_address)
    paired_smiles_list=open_lists[0]
    all_pairs_list=open_lists[1]

    print(len(paired_smiles_list_0))
    print(len(paired_smiles_list))
    print(len(all_pairs_list_0))
    print(len(all_pairs_list))

    count_x=0
    for one_pair in tqdm(all_pairs_list):
        if one_pair in all_pairs_list_0:
            count_x+=1
    print(count_x)


    return


##############################################################################################################
##############################################################################################################
def Step08_main(loading_folder, saving_folder, ECFP_encodings):
    Initialize_all_smiles_all_ecfps(loading_folder, saving_folder, ECFP_encodings)
    print("Step08 Done!")
    return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
if __name__ == '__main__':

    #test_difference_bwt_old_and_new()

    #Initialize_all_smiles_all_ecfps(loading_folder, saving_folder, ECFP_encodings)

    Step08_show_similarities_distn(loading_folder, saving_folder, ECFP_encodings)






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Useless funcs 

def update_all_smiles_all_ecfps(): # Outdated and Deprecated DO NOT USE !!!!!
    print ("start update_all_smiles_all_ecfps ...")
    open_pickle_paired_smiles_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    open_pickle_all_pairs_list="../zz_metric_learn_savings/zz0_all_pairs_list_000.pickle"
    open_lists = Step08_open_Step07_pickles(open_pickle_paired_smiles_list,open_pickle_all_pairs_list)
    paired_smiles_list=open_lists[0]
    all_pairs_list=open_lists[1]
    new_all_smiles = get_all_smiles(paired_smiles_list, all_pairs_list)
    print ("total compounds number: ", len(new_all_smiles))
    #====================================================================================================#
    pickle_in1=open("../zz_metric_learn_savings/zz1_all_cmpds_000","rb")
    all_smiles=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open("../zz_metric_learn_savings/zz1_all_ecfps_000","rb")
    all_ecfps=pickle.load(pickle_in2)
    pickle_in2.close()
    pickle_in3=open("../zz_metric_learn_savings/zz1_all_cmpds_ecfps_dict_000","rb")
    all_smiles_ecfps_dict=pickle.load(pickle_in3)
    pickle_in3.close()
    #====================================================================================================#
    #all_ecfps=update_all_ECFPs(new_all_smiles,all_smiles,all_ecfps,ecfp_type="ECFP2",iteration_number=1)
    #all_smiles_ecfps_dict=update_all_smiles_ecfps_dict(new_all_smiles,all_smiles,all_smiles_ecfps_dict,ecfp_type="ECFP2",iteration_number=1)
    (all_ecfps,all_smiles_ecfps_dict)=update_all_ecfps_and_all_smiles_ecfps_dict(new_all_smiles,all_smiles,all_ecfps,all_smiles_ecfps_dict,ecfp_type="ECFP2",iteration_number=1)
    #====================================================================================================#
    pickle_out1=open("../zz_metric_learn_savings/zz1_all_cmpds_000","wb")
    pickle.dump(new_all_smiles, pickle_out1)
    pickle_out1.close()
    pickle_out2=open("../zz_metric_learn_savings/zz1_all_ecfps_000","wb")
    pickle.dump(all_ecfps, pickle_out2)
    pickle_out2.close()
    pickle_out3=open("../zz_metric_learn_savings/zz1_all_cmpds_ecfps_dict_000","wb")
    pickle.dump(all_smiles_ecfps_dict, pickle_out3)
    pickle_out3.close()

    print ("DONE UPDATE update_all_smiles_all_ecfps")


##############################################################################################################
##############################################################################################################