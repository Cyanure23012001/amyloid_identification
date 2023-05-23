def hydro_calculations(freesasa_result):
    # Eisenberg scale 
    # source: https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
    hydrophobicity_scale = {
        'ALA': 0.62,
        'CYS': 0.29,
        'ASP': -0.9,
        'GLU': -0.74,
        'PHE': 1.19,
        'GLY': 0.48,
        'HIS': -0.40,
        'ILE': 1.38,
        'LYS': -1.50,
        'LEU': 1.06,
        'MET': 0.64,
        'ASN': -0.78,
        'PRO': 0.12,
        'GLN': -0.85,
        'ARG': -2.53,
        'SER': -0.18,
        'THR': -0.05,
        'VAL': 1.08,
        'TRP': 0.81,
        'TYR': 0.26
    }

    #source: DOI:10.2174/1874464810902030171
    hydropathy_scale = {
        'ALA': 1.8,
        'CYS': 2.5,
        'ASP': -3.5,
        'GLU': -3.5,
        'PHE': 2.8,
        'GLY': -0.4,
        'HIS': -3.2,
        'ILE': 4.5,
        'LYS': -3.9,
        'LEU': 3.8,
        'MET': 1.9,
        'ASN': -3.5,
        'PRO': -1.6,
        'GLN': -3.5,
        'ARG': -4.5,
        'SER': -0.8,
        'THR': -0.7,
        'VAL': 4.2,
        'TRP': -0.9,
        'TYR': -1.3   
    }

    # Calculate SASA for entire structure
    structure_sasa = freesasa_result
    test = structure_sasa.residueAreas()
    dicthydrophobicity = {}
    dicthydropathy = {}
    for x in test:
        for y in test[x]:
            #print(test[x][y].residueType)
            #print(hydrophobicity_scale[test[x][y].residueType])
            # Calcul de l'hydrophobicité relative à chaque résidu, en fonction du SASA
            # doi:10.1016/s0021-9673(03)00182-1. ISSN 0021-9673. PMID 12877193, source wikipedia hydrophobicity scales
            #print(float(hydrophobicity_scale[test[x][y].residueType])*test[x][y].total)
            dicthydrophobicity[f'{x}_{y}'] = float(hydrophobicity_scale[test[x][y].residueType])*test[x][y].total
            dicthydropathy[f'{x}_{y}'] = hydropathy_scale[test[x][y].residueType]

    total_hydrophobicity = sum(dicthydrophobicity.values())
    total_hydropathy = sum(dicthydropathy.values())
    return total_hydrophobicity,total_hydropathy
