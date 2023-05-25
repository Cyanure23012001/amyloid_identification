import align
import coordinates_editor
import hydro_calc
import pandas as pd
from pdb_numpy.format import mmcif, pdb
from pdb_numpy import Coor, Model
import pdb_numpy.DSSP as DSSP
import numpy as np
import freesasa
from collections import Counter
import api_informations
import local_DSSP
import time
import pdb_numpy
import download_assembly

def PDB(pdb_path):
    #print(pdb_path)
    #on récupere les données sur la pdb
    api_information = api_informations.get_info(pdb_path[-8:-4])
    #print(api_information)

    #Si cristallo: on reconstruit, sinon on recupère l'assembly
    if api_information['experimental_method'] == "X-RAY DIFFRACTION":
        coords = Coor(pdb_path)
        test = coords.select_atoms("protein and not altloc B C D E")

        sequenceAA = test.get_aa_seq()
        temp = {val: key for key, val in sequenceAA.items()}
        res = {val:key for key, val in temp.items()}    
        AAseq = ''.join(list(res.values())).replace('-','')

        test = coordinates_editor.cooredit(test)
        path_struct_complete = f'modified_structures/{pdb_path[-8:-4]}.pdb'
        test.write(path_struct_complete)
    
    else:
        pdb_path = download_assembly.download_pdb(pdb_path[-8:-4])
        coords = Coor(pdb_path)
        test = coords.select_atoms("protein and not altloc B C D E")

        # Generation de la sequence d'AA
        sequenceAA = test.get_aa_seq()
        temp = {val: key for key, val in sequenceAA.items()}
        res = {val:key for key, val in temp.items()}    
        AAseq = ''.join(list(res.values())).replace('-','')

    
    # Pour diminuer les temps de calculs/ remplacé par un timeout maintenant
    """
    if len(AAseq) > 500:
        #print(AAseq)
        print(f"CANCELLED {pdb_path}, {len(AAseq)} AA")
        return None
    """

    print(pdb_path[-8:-4])
    api_information = api_informations.get_info(pdb_path[-8:-4])
    #print(api_information)


    #print(f"CHAINS = {str(test.chain)}")
    #coor.write_pdb('7ZH7_new.pdb')
    #print(test.transformation)

    # Calcul de la matrice de liaisons hydrogènes
    #try:
    try:
        HBOND = DSSP.compute_Hbond_matrix(test).view()
    except Exception as e:
        with open("error.txt", "a") as f:
            f.write(f"{pdb_path[-8:-4]} {str(e)}\n")
        return None

    
    #result_dssp = DSSP.compute_DSSP(test)
    result_dssp = local_DSSP.DSSP(HBOND,test)
    
#print(len(HBOND))


    df = pd.DataFrame(HBOND)
    # Get number of cells in df dataframe
    number_max_of_hbonds = df.size




# Calcul des noms de residus correspondants à chaque chaine

    result = pd.unique(["%d_%s" % t for t in zip(test.uniq_resid, test.chain)])
    #print(result)


#print(result)
#print(df)

    df.columns = result
    df.index = result

# Matrice Oui/non sur liaisons hydrogènes mais pour les résidus
    coord = np.where(df)
    coordinates = [(x,y) for x, y in zip(coord[0], coord[1])]
#print(coordinates)
    names = pd.unique(test.chain)
    df2 = pd.DataFrame(0, index=names, columns=names)

    for coord in coordinates:
        i, j = coord
    #print(f"Coords: {coord}, Index: {df.index[i]}, Column: {df.columns[j]}")
        row = df.index[i][-1:]
        column = df.columns[j][-1:]
    #print(f"Row: {row}, Column = {column}")
        df2[column][row] += 1

    #d = dtale.show(df,notebook=True)
    #d = dtale.show(df2)
    #print(d._main_url)

    unique_combinations = []

    # Determine le nombre de fois ou il y a des liaisons (matrice qui n'est pas symmetrique)
    for i in names:
        for j in names:
            pair = [i,j]
            if i != j and (pair[::-1] not in unique_combinations) :
                unique_combinations.append(pair)

#print(unique_combinations)
    total = {}
    numberofh = 0
    
    # On récupère toutes les liaisons différentes trouvées
    for combinations in unique_combinations:
        total[f"{combinations[0]}_{combinations[1]}"] = int(df2[combinations[0]][combinations[1]])+int(df2[combinations[1]][combinations[0]])
        numberofh = numberofh + int(df2[combinations[0]][combinations[1]])+int(df2[combinations[1]][combinations[0]])


    numbertotalh = numberofh/len(pd.unique(test.uniq_resid))*100


    numberHchains = numberofh/(len(pd.unique(test.uniq_resid)) * len(pd.unique(test.chain)))

    # Recuperation de la nouvelle sequence
    sequenceAA_full = test.get_aa_seq()
    temp = {val: key for key, val in sequenceAA_full.items()}
    res = {val:key for key, val in temp.items()}    
    AAseq_full = ''.join(list(res.values())).replace('-','')
    #except:
     #   return None
    #Ecriture du fichier PDB
    #pathpdb = f"DOWN/{pdb_code}.pdb"
    #test.write(pathpdb)

    #Calcul avec FreeSasa du SASA

    # faire sur monomère et structure reconstruite

    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    #print(result.totalArea())
    area_classes = freesasa.classifyResults(result, structure)

    if api_information['experimental_method'] == "X-RAY DIFFRACTION":
        structure_total = freesasa.Structure(path_struct_complete)
        result_full = freesasa.calc(structure_total)
        result_total = result_full.totalArea()
    else:
        result_total = 0


    # Calcul de la frequence de chaque lettre avec DSSP
    df3=pd.DataFrame.from_dict(result_dssp[0],orient="index", columns=["DSSP"])
    conc = ""
    for item in result_dssp[0]:
        conc += result_dssp[0][item]
    conc = conc.replace(" ", "")
    #print(conc)
    DSSP_letters = ["G","H","I","T","E","B","S"]
    letter_freq = Counter(conc)


    for lettre in DSSP_letters:
        if lettre not in letter_freq:
            letter_freq[lettre] = 0
        else:
            letter_freq[lettre] = letter_freq[lettre]/len(conc)

    df4 = pd.DataFrame()

    #Alignement du DSSP

    if api_information['experimental_method'] == "X-RAY DIFFRACTION":
        total, alignement_total = align.alignement(result_dssp)
    else:
        total, alignement_total = align.alignement(result_dssp)
        #alignement_total = 0
    try:
        total_hydrophobicity,total_hydropathy = hydro_calc.hydro_calculations(result)
    except Exception as e:
        with open("error.txt", "a") as f:
            f.write(f"{pdb_path[-8:-4]} {str(e)}\n")
        return None

    descriptors = {'PDB_ID':pdb_path[-8:-4], 'nb H_bond / nres':numbertotalh,'nb H_bond total':numberofh, 'H_bond / nb_res * nb_chains':numberHchains, 'Max_hbonds_number':number_max_of_hbonds, 'Number of chains':len(names), 'SASA_unedited':result.totalArea(), 'SASA_full':result_total, 'SASA_theo':result.totalArea()*len(names), 'Alignement_total':total, 'Alignement_moyen': alignement_total, 'Hydrophobicity':total_hydrophobicity, 'Hydropathy':total_hydropathy,'AAseq':AAseq, 'AAseq_full':AAseq_full, 'DF_DSSP':df3,'DF_H':df2}
    new_row = {**descriptors,**api_information,**letter_freq}
    
    #print(new_row)
    return new_row
