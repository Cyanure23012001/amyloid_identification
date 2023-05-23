import pdb_numpy
import numpy as np
from pdb_numpy import geom
from scipy.spatial import distance_matrix

def angle_vec(vec_a, vec_b):
    unit_vec_a = vec_a / np.linalg.norm(vec_a)
    unit_vec_b = vec_b / np.linalg.norm(vec_b)

    dot_product = np.dot(unit_vec_a, unit_vec_b)

    angle = np.arccos(dot_product)

    return angle


def compute_bend(CA_sel):
    n_res = len(CA_sel.uniq_resid)
    bend = np.array([False for _ in range(n_res)])

    for i in range(2, n_res - 2):
        CA_i = CA_sel.xyz[i]
        CA_i_minus_2 = CA_sel.xyz[i - 2]
        CA_i_plus_2 = CA_sel.xyz[i + 2]

        vec_i_1 = CA_i - CA_i_minus_2
        vec_i_2 = CA_i_plus_2 - CA_i

        vec_i_1 /= np.linalg.norm(vec_i_1)
        vec_i_2 /= np.linalg.norm(vec_i_2)

        # print(i, np.degrees(geom.angle_vec(vec_i_1, vec_i_2)))
        if np.degrees(angle_vec(vec_i_1, vec_i_2)) > 70.0:
            bend[i] = True

    return bend

#@codon.jit

def DSSP(Hbond_mat,coor):   
    #from python import numpy as np 
    model = coor.models[0]
    cutoff = 8

    CA_sel = coor.select_atoms("protein and name CA and not altloc B C D E")
    unique_residues = np.unique(CA_sel.uniq_resid)
    chain_array = CA_sel.chain
    n_res = len(unique_residues)

    max_dist = 0

    SS_list = []
    SS_seq = np.array([" " for i in range(n_res)])
    H_seq = np.array([False for i in range(n_res)])
    G_seq = np.array([False for i in range(n_res)])
    I_seq = np.array([False for i in range(n_res)])
    E_seq = np.array([False for i in range(n_res)])
    S_seq = compute_bend(CA_sel.models[0])

    # N-turn
    for i in range(n_res - 3):
        if i < n_res - 4 and Hbond_mat[i, i + 4]:
            H_seq[i] = True
        if Hbond_mat[i, i + 3]:
            G_seq[i] = True
        if i < n_res - 5 and Hbond_mat[i, i + 5]:
            I_seq[i] = True

    # Beta sheet
    # PART TO ACCELERATE
    for i in range(1, n_res - 1):
        for j in range(i + 3, n_res - 1):
            if (
                (Hbond_mat[i - 1, j] and Hbond_mat[j, i + 1])
                or (Hbond_mat[j - 1, i] and Hbond_mat[i, j + 1])
            ) or (
                (Hbond_mat[i, j] and Hbond_mat[j, i])
                or (Hbond_mat[i - 1, j + 1] and Hbond_mat[j - 1, i + 1])
            ):
                E_seq[i] = True
                E_seq[j] = True

    # Assign secondary structure sequence (order follows the list above)
    # Bend
    for i in range(n_res):
        if S_seq[i]:
            SS_seq[i] = "S"

    for i in range(n_res - 1):
        if I_seq[i]:
            SS_seq[i + 1 : i + 5] = "T"
        if H_seq[i]:
            SS_seq[i + 1 : i + 4] = "T"
        if G_seq[i]:
            SS_seq[i + 1 : i + 3] = "T"

    # A minimal helix is defined by two consecutive n-turns

    for i in range(1, n_res - 3):
        if G_seq[i] and G_seq[i - 1]:
            SS_seq[i : i + 3] = "G"

    for i in range(1, n_res - 1):
        if E_seq[i] and E_seq[i - 1]:
            SS_seq[i - 1 : i + 1] = "E"
        elif E_seq[i] and not (E_seq[i + 1] and E_seq[i - 1]):
            SS_seq[i] = "B"

    for i in range(1, n_res - 4):
        # A minimal helix is defined by two consecutive n-turns
        if H_seq[i] and H_seq[i - 1]:
            SS_seq[i : i + 4] = "H"
    #  In 2012, DSSP was rewritten so that the assignment of π helices
    #  was given preference over α helices, resulting in better detection of π helices.
    for i in range(1, n_res - 5):
        if I_seq[i] and I_seq[i - 1]:
            SS_seq[i : i + 5] = "I"

    # Two overlapping minimal helices offset by two or three residues are joined
    # into one helix:
    for i in range(1, n_res - 2):
        if SS_seq[i - 1] == "H" and (SS_seq[i + 1] == "H" or SS_seq[i + 2] == "H"):
            SS_seq[i] = "H"

    for i in range(1, n_res - 1):
        if SS_seq[i - 1] in ["E", "B"] and (SS_seq[i + 1] in ["E", "B"]):
            SS_seq[i - 1 : i + 2] = "E"

    seq_dict = {}

    for SS, chain in zip(SS_seq, chain_array):
        if chain not in seq_dict:
            seq_dict[chain] = SS
        else:
            seq_dict[chain] += SS

    SS_list.append(seq_dict)
    return SS_list
