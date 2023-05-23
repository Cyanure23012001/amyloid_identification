import os
import requests
def download_pdb(pdb_id):
    # Define the URL to download the PDB assembly file
    local_filename = f"download_folder/{pdb_id}.pdb"
    if not os.path.exists(local_filename):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

        # Use the requests module to get the content of the URL
        response = requests.get(url)


        # Save the content of the response to the local filename
        with open(local_filename, "wb") as f:
            f.write(response.content)



    #print(f"Downloaded {local_filename}")
    return local_filename

