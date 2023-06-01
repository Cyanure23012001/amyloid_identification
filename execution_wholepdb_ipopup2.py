import pandas as pd
import numpy as np
import dtale
import pickle
import freesasa
import concurrent.futures
import os
from tqdm import tqdm
import copy
import glob
import re
import time
import gzip
import shutil
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError

# Custom modules
import calculations # Calcul de pleins d'informations
import download # Télécharge les structures

# Fonction executant le téléchargement et le calcul des descripteurs


def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        shutil.copyfileobj(s_file, d_file, block_size)
"""
with open("list_files.txt", "w") as f:
    n1 = '\n'
    f.write(f"{n1.join(glob.glob('/amyloid_identification/pdb_db/*'))}")
    f.write(f"{len(glob.glob('/amyloid_identification/pdb_db/*'))}")
"""

def process_struct(struct, is_amyloid):
	#try:
    try:
	#/amyloid_identification/pdb_db/pdb3rrr.ent.gz
        first_pdb_path = f"/pdb_db/all/pdb/pdb{struct.lower()}.ent.gz"
        second_pdb = os.readlink(first_pdb_path)
        third_pdb = os.path.abspath(second_pdb)
        pdb_path = f"/pdb_db{third_pdb}"
        #print(third_pdb)
        #print(pdb_path)
        if os.path.exists(pdb_path):
            print("EXISTE")
            with open("logging/file_exist.txt", "a") as f:
            	f.write(f"{pdb_path}\n")
	
        local_filename = f"{os.getcwd()}/download_folder/{struct}.pdb"
        print(local_filename)
        gunzip_shutil(pdb_path, local_filename)

        #local_filename = download.download_pdb(struct)
        new_row = calculations.PDB(local_filename)
    except Exception as e:
        with open("logging/error.txt", "a") as f:
            f.write(f"{pdb_path} | {pdb_path[-8:-4]} {str(e)}\n")
        new_row = None
        
    if new_row != None:
        amyloid = {'is_amyloid': is_amyloid}
        new_row2 = {**new_row,**amyloid}

    
        # delete local_filename, test if pdb exist in download_folder_assembly, in modified structures and also deleted them
        try:
            if os.path.exists(local_filename):
                os.remove(local_filename)
            if os.path.exists(f"download_folder_assembly/{struct}.pdb"):
                os.remove(f"download_folder_assembly/{struct}.pdb.gz")
                os.remove(f"download_folder_assembly/{struct}.pdb")
            if os.path.exists(f"modified_structures/{struct}.pdb"):
                os.remove(f"modified_structures/{struct}.pdb")
        except:
            pass

        
        #print(new_row2)
        #print(f"{local_filename} succeed")
        return new_row2

    else:
        try:
            if os.path.exists(local_filename):
                os.remove(local_filename)
            if os.path.exists(f"download_folder_assembly/{struct}.pdb"):
                os.remove(f"download_folder_assembly/{struct}.pdb.gz")
                os.remove(f"download_folder_assembly/{struct}.pdb")
            if os.path.exists(f"modified_structures/{struct}.pdb"):
                os.remove(f"modified_structures/{struct}.pdb")
        except:
            pass

        print(f"{local_filename} failed")
        return None

df4_list = []

df4 = pd.DataFrame()
# Chargement des fichiers de données, ici les resultats de la requete
with open('pickle_models/found_pdbs.pickle', 'rb') as f:
    res = pickle.load(f)
print(f"Res size: {len(res)}")
with open("nb_file.log", "a") as f:
    f.write(f"{len(res)}\n")



already_computed = glob.glob('pickle_models/calculations_results_temp_*')
if len(already_computed) > 1:
    numb = 0
    for x in already_computed:
        curr_number = re.findall(r'\d+',x)
        if int(curr_number[0]) > numb:
            numb = int(curr_number[0])

    
    # On ouvre les fichiers pickle et on soustrait les valeurs de la colonne 'PDB_ID' à la liste res
    # Permet de skip les valeurs deja calculées
    with open(f'pickle_models/calculations_results_temp_{numb}.pickle','rb') as f:
        temp_df = pickle.load(f)
        

    res = list(set(res) - set(temp_df['PDB_ID'].tolist()))
else:
    temp_df = pd.DataFrame()

print(f"Res size: {len(res)}")
with open("logging/nb_file.log", "a") as f:
    f.write(f"{len(res)}\n")

with ProcessPool(max_workers=os.cpu_count()) as pool:
    futures = [pool.schedule(process_struct, args=(struct, True), timeout=60) for struct in res]
    with tqdm(total=len(futures)) as pbar:
        for i, future in enumerate(futures):
            new_row2 = None
            try:
                new_row2 = future.result()
            except TimeoutError:
                print(f"Killed worker due to timeout: {future}")
            except ProcessExpired as error:
                print(f"Process expired: {error}")
            if new_row2 is not None:
                df4_list.append(new_row2)
            
            # Systeme de backup des résultats
            if i % 100 == 0:
                try:
                    print("SAVING_results")
                    results_list = [i for i in df4_list if i is not None]
                    df4 = pd.DataFrame.from_dict(df4_list)

                    # on y sauvegarde les résultats obtenus précédemment
                    results = pd.concat([temp_df, df4], ignore_index=True)
                
                    #results = results[[col for col in results if not (col.startswith('DF') or col.startswith('AA'))]+[col for col in results if (col.startswith('DF') or col.startswith('AA'))]]
                    with open(f'pickle_models/calculations_results_temp_{i}.pickle','wb') as f:
                        pickle.dump(results,f)
                except Exception as e:
                    with open("logging/error_2.txt", "a") as f:
                        f.write(f"{i} | {str(e)}\n")

            pbar.update(1)


 

#

# On enlève les valeurs où il y a des None
results_list = [i for i in df4_list if i is not None]


results = pd.DataFrame.from_dict(df4_list)
 
# On deplace à droite les colonnes prenant beaucoup de place
results = results[[col for col in results if not (col.startswith('DF') or col.startswith('AA'))]+[col for col in results if (col.startswith('DF') or col.startswith('AA'))]]


# On sauvegarde au format pickle le résultat

with open('pickle_models/calculations_results_requesttest','wb') as f:
    pickle.dump(results,f)

