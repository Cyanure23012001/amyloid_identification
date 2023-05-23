import download
import calculations
struct = '4E0O'

#@profile
def exec(struct):
    local_filename = download.download_pdb(struct)
    new_row = calculations.PDB(local_filename)
    return new_row

print(exec(struct))

