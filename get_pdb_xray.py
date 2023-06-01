from pypdb import *
found_pdbs = Query("X-RAY DIFFRACTION", query_type='ExpTypeQuery').search()
print(found_pdbs[:10])
 
import pickle
with open('found_pdbs.pickle', 'wb') as f:
    pickle.dump(found_pdbs, f)



