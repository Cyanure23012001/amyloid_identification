import pandas as pd
import numpy as np
from Bio import Align
from Bio.Align import PairwiseAligner
import pdb_numpy.DSSP as DSSP
import pdb_numpy
import dtale
import pickle
from pdb_numpy.format import mmcif, pdb
from pdb_numpy import Coor, Model
import pdb_numpy.DSSP as DSSP
import freesasa
from collections import Counter
import time
import requests
import concurrent.futures
import gzip
import os
from tqdm import tqdm
import copy
import logging

# Custom modules
import align # Alignement avec Needleman et Wunsch du DSSP
import calculations # Calcul de pleins d'informations
import coordinates_editor # Calcul des coordonnées 3D
import hydro_calc # Calcul de l'hydropathie et de l'hydrophobicité
import download # Télécharge les structures

def process_struct(struct, is_amyloid):
	#try:
    try:
        local_filename = download.download_pdb(struct)
        new_row = calculations.PDB(local_filename)
    except:
        new_row = None
        
    if new_row != None:
        amyloid = {'is_amyloid': is_amyloid}
        new_row2 = {**new_row,**amyloid}
        #print(new_row2)
        #print(f"{local_filename} succeed")
        return new_row2
    else:
        print(f"{local_filename} failed")
        return None

df4_list = []

df4 = pd.DataFrame()
# Chargement des fichiers de données
with open('pickle_models/CPAD.pickle', 'rb') as f:
    CPAD_list = pickle.load(f)

print(CPAD_list)

from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError

with ProcessPool(max_workers=os.cpu_count()) as pool:
    futures = [pool.schedule(process_struct, args=(struct, True), timeout=20) for struct in CPAD_list]

    for future in futures:
        new_row2 = None
        try:
            new_row2 = future.result()
        except TimeoutError:
            print(f"Killed worker due to timeout: {future}")
        except ProcessExpired as error:
            print(f"Process expired: {error}")
        if new_row2 is not None:
            df4_list.append(new_row2)


"""

for x in tqdm(CPAD_list):
    
    local_filename = download.download_pdb(x)
    new_row = calculations.PDB(local_filename)
    df4_list.append(new_row)
"""

df4_list = [i for i in df4_list if i is not None]


df4 = pd.DataFrame.from_dict(df4_list)
 
df4 = df4[[col for col in df4 if not (col.startswith('DF') or col.startswith('AA'))]+[col for col in df4 if (col.startswith('DF') or col.startswith('AA'))]]


with open('CPAD_calculations','wb') as f:
    pickle.dump(df4,f)
"""
with open('tempresults2', 'rb') as f:
    tempresults2 = pickle.load(f)


frames = [tempresults2, df4]
tempresults3 = pd.concat(frames)

with open('tempresults3','wb') as f:
    pickle.dump(tempresults3,f)
"""
print(df4)
import dtale
d = dtale.show(df4,subprocess=False)
print(d._url)
