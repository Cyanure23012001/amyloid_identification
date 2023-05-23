from Bio import Align
from Bio.Align import PairwiseAligner
import numpy as np
#@profile
def alignement(result_dssp):
    aligner = PairwiseAligner()
    count = []
    longueur = []
    chains_letters = []

    for chain in result_dssp[0]:
        chains_letters.append(chain)
        #print(result_dssp[0][chain])
        # E beta ladder
        # B beta bridge
        # S Bend
        count.append(result_dssp[0][chain].count('E'))
        longueur.append(len(result_dssp[0][chain]))
    #print(len(chains_letters))

    alignements_results = []
    #print(chains_letters)
    for value in range(1,len(chains_letters)):
        alignement = aligner.align(result_dssp[0][chains_letters[0]],result_dssp[0][chains_letters[value]])
        #print(alignement.score)
        bit_score = (alignement.score / len(result_dssp[0][chains_letters[0]])*100)
        #print(bit_score)
        alignements_results.append(bit_score)

    countnp = np.array(count)
    longueurnp = np.array(longueur)


    resultats_alignements_np = np.array(alignements_results)


    total = (np.sum(countnp)/np.sum(longueurnp))*100

# On enlève les valeurs qui ne sont pas numériques ?
    #print(resultats_alignements_np)
    #resultats_filtered = resultats_alignements_np[~np.isnan(resultats_alignements_np).any(axis=1)]
    #print(resultats_filtered)
    alignement_moyen = np.mean(resultats_alignements_np)

    #print(total)
    #print(alignement_moyen)

    return total, alignement_moyen
