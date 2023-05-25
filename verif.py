 
import urllib.request
import tkinter as tk
import pandas
import pickle
with open('pickle_models/calculations_full', 'rb') as f:
    resultdf = pickle.load(f)
 
# open the pdb_save_list pickle, get the first element for each list in the list, and substract theses elements from the pdb_list var
try:
    with open('pdb_save_list.pickle', 'rb') as f:
        pdb_save_list = pickle.load(f)
    pdb_list = resultdf['PDB_ID'].values.tolist()
    pdb_list = [x for x in pdb_list if x not in [i[0] for i in pdb_save_list]]
except FileNotFoundError:
    pdb_list = resultdf['PDB_ID'].values.tolist()





from PIL import Image, ImageTk
if len(pdb_save_list) == 0:
	pdb_save_new = []
else:
	pdb_save_new = pdb_save_list

quit_var = 0
for pdb_code in pdb_list:
    if quit_var == 1:
        break
    url = f"https://cdn.rcsb.org/images/structures/{pdb_code.lower()}_assembly-1.jpeg"
    try:
        urllib.request.urlretrieve(url, f"{pdb_code}.jpeg")

        # show the image

        def save_pdb_code(answer):
            pdb_save_new.append([pdb_code, answer])
            root.destroy()

        root = tk.Tk()
        root.title("Save PDB Code")
        root.geometry("500x600")
        img = Image.open(f"{pdb_code}.jpeg")
        img = img.resize((500, 500), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(img)
        panel = tk.Label(root, image=img)
        panel = tk.Label(root, image=img)
        panel.pack(side="top", expand="no")
        yes_button = tk.Button(root, text="Yes", command=lambda: save_pdb_code(True))
        yes_button.pack(side=tk.LEFT, padx=10)

        no_button = tk.Button(root, text="No", command=lambda: save_pdb_code(False))
        no_button.pack(side=tk.RIGHT, padx=10)
        def break_boucle():
            with open('pdb_save_list.pickle', 'wb') as f:
                pickle.dump(pdb_save_new, f)    
            print("quit")
            global quit_var 
            quit_var = 1
        exit_button = tk.Button(root, text="Exit", command=lambda: break_boucle())
        exit_button.pack(side=tk.RIGHT, padx=50)


        
    

        root.mainloop()

        """
        user_input = input("Image downloaded successfully. Do you want to save this PDB code? (y/n): ")
        if user_input.lower() == 'y':
            pdb_save_list.append([pdb_code,user_input])
        if user_input.lower() == 'n':
            pdb_save_list.append([pdb_code,user_input])
        """
    except:
        print("Error: Image could not be downloaded.")

    


print("PDB codes with downloaded images:", pdb_save_list)
with open('pdb_save_list.pickle', 'wb') as f:
    pickle.dump(pdb_save_new, f)  


 
resultdf.loc[resultdf['PDB_ID'] == '7MU1', 'is_amyloid'] = True
 
resultdf.loc[resultdf['experimental_method'] != 'X-RAY DIFFRACTION', 'SASA_full'] = resultdf['SASA_unedited']




