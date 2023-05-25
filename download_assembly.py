import os
import requests
def download_pdb(pdb_id):
    # Define the URL to download the PDB assembly file
    local_filename = f"download_folder_assembly/{pdb_id}.pdb.gz"
    if not os.path.exists(local_filename):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"

        # Use the requests module to get the content of the URL
        response = requests.get(url)

        # Save the content of the response to the local filename
        with open(local_filename, "wb") as f:
            f.write(response.content)
        
        import gzip
        import shutil

        def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
            with gzip.open(source_filepath, 'rb') as s_file, \
                    open(dest_filepath, 'wb') as d_file:
                shutil.copyfileobj(s_file, d_file, block_size)

        gunzip_shutil(local_filename, local_filename[:-3])

        




    #print(f"Downloaded {local_filename}")
    return local_filename[:-3]

