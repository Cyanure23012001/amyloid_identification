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

    local_filename = download.download_pdb(struct)
    new_row = calculations.PDB(local_filename)
    if new_row != None:
        amyloid = {'is_amyloid': is_amyloid}
        new_row2 = {**new_row,**amyloid}
        #print(new_row2)
        #print(f"{local_filename} succeed")
        return new_row2
    else:
        print(f"{local_filename} failed")
        return None
	#except:
	#return None
"""

def process_struct(struct, is_amyloid):

	local_filename = download.download_pdb(struct)
	new_row = calculations.PDB(local_filename)
	amyloid = {'is_amyloid': is_amyloid}
	new_row2 = {**new_row,**amyloid}
	print(new_row2)
	return new_row2
"""


start = time.time()
df4_list = []

df4 = pd.DataFrame()
# Chargement des fichiers de données
with open('mehdi.pickle', 'rb') as f:
    mehdi_list = pickle.load(f)

with open('TESTSET_500.pickle', 'rb') as f:
    test_set = pickle.load(f)

with open('tempresults', 'rb') as f:
    tempresults = pickle.load(f)


print(len(df4))
tempresults = pd.DataFrame()
df4 = copy.deepcopy(tempresults)

#print(len(tempresults))

# On regarde ce qui n'a pas encore été calculé
if len(tempresults) != 0:
    list_pdb_calcules_mehdi = tempresults[tempresults['is_amyloid'] == True]['PDB_ID'].values.tolist()
    mehdi_list = [x for x in mehdi_list if x not in list_pdb_calcules_mehdi]
    print(len(mehdi_list))

# Permet d'obtenir tout les pdb qui n'ont pas été calculés pour le moment
    list_pdb_calcules = tempresults[tempresults['is_amyloid'] == False]['PDB_ID'].values.tolist()
    test_set = [x for x in test_set if x not in list_pdb_calcules]
    print(len(test_set))
else:
    print(len(mehdi_list))
    print(len(test_set))


#print(process_struct('6SJX',True))

# Multiprocess execution
"""
with concurrent.futures.ProcessPoolExecutor() as executor:

    # Process the `mehdi_list` items
    results = [executor.submit(process_struct, struct, True) for struct in mehdi_list]
    for future in concurrent.futures.as_completed(results):
        new_row2 = None
        try:
            new_row2 = future.result(timeout=60) # Timeout set to 1 minute
        except concurrent.futures.TimeoutError:
            future.cancel() # Kill the worker if the timeout is reached
            print(f"Killed worker due to timeout: {future}")
        if new_row2 is not None:
            df4_list.append(new_row2)

    # Process the `test_set` items
    results = [executor.submit(process_struct, struct, False) for struct in test_set]
    for future in concurrent.futures.as_completed(results):
        new_row2 = None
        try:
            new_row2 = future.result(timeout=60) # Timeout set to 1 minute
        except concurrent.futures.TimeoutError:
            future.cancel() # Kill the worker if the timeout is reached
            print(f"Killed worker due to timeout: {future}")
        if new_row2 is not None:
            df4_list.append(new_row2)
"""

"""
for x in tqdm(mehdi_list):
    try:
        df4_list.append(process_struct(x,True))
    except:
        pass

for x in tqdm(test_set):
    try:
        df4_list.append(process_struct(x,False))
    except:
        pass

#print(df4_list)

# Suppression des NoneTypes
df4_list = [i for i in df4_list if i is not None]
"""
"""
newusers = dict()
for ud in users:
    newusers[ud.pop('id')] = ud
print newusers
"""

"""
df4 = pd.DataFrame.from_dict(df4_list)
dtale.show(df4)
with open('tempresults2','wb') as f:
    pickle.dump(df4,f)
"""

df4_list = []
with open('results.pickle', 'rb') as f:
    res = pickle.load(f)

for x in tqdm(res):
    try:
        local_filename = download.download_pdb(x)
        new_row = calculations.PDB(local_filename)
        df4_list.append(new_row)
    except:
        pass
df4_list = [i for i in df4_list if i is not None]
df5 = pd.DataFrame.from_dict(df4_list)
with open('calculations_1','wb') as f:
    pickle.dump(df5,f)

end = time.time()
print("FIN")
print(end - start)
print(d._url)


