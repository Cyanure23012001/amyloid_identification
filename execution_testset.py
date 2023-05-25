import pandas as pd
import numpy as np
import dtale
import pickle
import freesasa
import concurrent.futures
import os
from tqdm import tqdm
import copy

from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError

# Custom modules
import calculations # Calcul de pleins d'informations
import download # Télécharge les structures

# Fonction executant le téléchargement et le calcul des descripteurs

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
# Chargement des fichiers de données, ici les resultats de la requete
with open('pickle_models/results.pickle', 'rb') as f:
    res = pickle.load(f)


print(len(res))
import time
time.sleep(5)

with ProcessPool(max_workers=os.cpu_count()) as pool:
    futures = [pool.schedule(process_struct, args=(struct, True), timeout=60) for struct in res]

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
    
 

#

# On enlève les valeurs où il y a des None
results_list = [i for i in df4_list if i is not None]


results = pd.DataFrame.from_dict(df4_list)
 
# On deplace à droite les colonnes prenant beaucoup de place
results = results[[col for col in results if not (col.startswith('DF') or col.startswith('AA'))]+[col for col in results if (col.startswith('DF') or col.startswith('AA'))]]


# On sauvegarde au format pickle le résultat

with open('pickle_models/calculations_results_request','wb') as f:
    pickle.dump(results,f)

