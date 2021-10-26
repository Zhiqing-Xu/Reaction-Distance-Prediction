# ==================================================================================== #
#!/usr/bin/python

# ==================================================================================== #
# header
import os, os.path
from sys import platform
if os.name == 'nt' or platform == 'win32':
    os.chdir(os.path.dirname(__file__))

import numpy as np
import json
import pickle
from model import MoleculeVAE
from utils import encode_smiles, decode_latent_molecule, interpolate
from math import *

def euclidean_distance(x,y):
    return sqrt(sum(pow(a-b,2) for a, b in zip(x, y)))

def smiles_to_latent_vec(one_smiles, model, charset):
    def modify_smiles(one_smiles):
        # Remove Radicals *, since it is not trained in the VAE model
        # Delete in the following order so that left chars are correct SMILES strings
        one_smiles=one_smiles.replace("([*])","")
        one_smiles=one_smiles.replace("(*)","")
        one_smiles=one_smiles.replace("[*]","")
        one_smiles=one_smiles.replace("(-*)","")
        one_smiles=one_smiles.replace("-*","")
        one_smiles=one_smiles.replace("\*","")
        one_smiles=one_smiles.replace("/*","")
        one_smiles=one_smiles.replace("*","")
        one_smiles=one_smiles.replace(":","")
        return one_smiles
    latent_vec = encode_smiles(modify_smiles(one_smiles), model, charset)
    return list(latent_vec[0])

def latent_vec_to_smiles(one_latent_vec, model, charset, latent_dim=292):
    smiles = decode_latent_molecule(np.array([one_latent_vec]), model, charset, latent_dim)
    return smiles

def smiles_list_to_latent_vec(smiles_list, model, charset):
    smiles_latent_dict=dict([])
    count_x=0
    for one_smiles in smiles_list:
        if len(one_smiles)<=120:
            count_x+=1
            print(count_x)
            smiles_latent_dict[one_smiles]=smiles_to_latent_vec(one_smiles, model, charset)
    return smiles_latent_dict



def main():
    #-------------------- (0) --------------------#
    # Number of dimensions trained the VAE
    # ref: http://chembl.blogspot.com/2017/07/using-autoencoders-for-molecule.html
    # trained_model 0.99 validation accuracy
    # trained with 80% of ALL chembl molecules, validated on the other 20.
    latent_dim = 292
    trained_model = 'chembl_23_model.h5'
    charset_file = 'charset.json'

    with open('charset.json', 'r') as outfile:
        charset = json.load(outfile)
    model = MoleculeVAE()
    model.load(charset, trained_model, latent_rep_size = latent_dim)



    #-------------------- (1) --------------------#
    pickle_in1=open("./Step03_Used_SMILES_set","rb")
    Used_SMILES_set=pickle.load(pickle_in1, encoding='latin1')
    pickle_in1.close()
    smiles_VAEVEC_dict = smiles_list_to_latent_vec(list(Used_SMILES_set), model, charset)
    print ("Done Conversion: smiles_list_to_latent_vec")
    
    pickle_out1=open("./Step04_smiles_VAEVEC_dict","wb")
    pickle.dump(smiles_VAEVEC_dict,pickle_out1)
    pickle_out1.close()
    
    #temp_file=open("./Step04_smiles_VAEVEC_dict_not_pickle","w")
    #temp_file.write(str(smiles_VAEVEC_dict))

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':

    main()