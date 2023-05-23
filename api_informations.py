# Module permetttant d'avoir des informations sur la structure actuelle, avec l'api de la RCSB
import requests

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

    return informations

print(get_info('4HHB'))


    