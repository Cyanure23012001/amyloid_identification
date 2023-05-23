# Module permetttant d'avoir des informations sur la structure actuelle, avec l'api de la RCSB
import requests
import re

def get_info(pdb_id):
    # URL pour récupérer les informations de base de la structure
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    # Récupération des données
    response = requests.get(url)
    data = response.json()
    # Initialisation du dictionnaire d'informations
    informations = {}
    # Récupération de la méthode expérimentale utilisée
    informations['experimental_method'] = data['exptl'][0]['method']
    # Récupération de l'ID de l'assemblage
    assembly_id = data['rcsb_entry_container_identifiers']['assembly_ids'][0]

    list_pdb_fasta = []
   
    url_fasta = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url_fasta)

    list_lines = response._content.decode().split('\n')
    for line in range(len(list_lines)):
        if list_lines[line].startswith('>'):
            fasta_dict = {}
            list_line = list_lines[line].split('|')
            pdb_id_local = list_line[0][1:]

            pattern = r'[A-Z]'
            # Find all matches in the input string
            matches = re.findall(pattern, list_line[1])

            fasta_dict['pdb_id_local'] = pdb_id_local
            fasta_dict['chains'] = matches[1:]
            fasta_dict['seq'] = list_lines[line+1]
            list_pdb_fasta.append(fasta_dict)
            
    # URL pour récupérer les informations sur l'assemblage
    url_assembly = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}"
    # Récupération des données
    response = requests.get(url_assembly)
    data = response.json()
    # Récupération de l'état oligomérique
    informations['oligomeric_state1'] = data['pdbx_struct_assembly']['oligomeric_details']
    informations['oligomeric_state2'] = data['rcsb_struct_symmetry'][0]['oligomeric_state']
    informations['atom_count'] = data['rcsb_assembly_info']['atom_count']
    informations['H_atom_count'] = data['rcsb_assembly_info']['hydrogen_atom_count']
    informations['symmetry'] = data['rcsb_struct_symmetry'][0]['type']
    informations['fasta'] = list_pdb_fasta

    return informations

print(get_info('7ZH7'))


    