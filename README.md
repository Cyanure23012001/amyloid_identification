
# Amyloid Identification

This folder contains scripts to aid in the identification of Beta-Amyloid structures using simple descriptors. The calculations include:

- Number of hydrogen bonds between chains, corrected with the number of chains

## Execution

To ensure everything works, execute `pip install -r requirements.txt` and `python setup.py`. The `execution*.py` files are the main scripts used to launch calculations on a specific set of PDB_IDs. They contain almost all the same code, with the most recent one being `execution_testset.py`, which calculates the descriptors for 3300 structures. The results are saved in a `.pickle` format in the `pickle_models` folder. Note that executing one of these scripts will replace the results.

To train and test the model, Jupyter notebooks are used. `Model_training` generates training and test sets from the data of `execution2.py`, while `Model_test` tests the model on a larger set (the results of `execution_testset.py`).

`Correct_amyloid_values` edits some values from CPAD with manual correction to select only beta amyloids. A little program, `verif.py`, is used to do the verification part.

In the `pickle_models` folder:

- `amypro.pickle`, `CPAD.pickle`, and `mehdi.pickle` contain PDB IDs and are used to train the model.
- `results.pickle` contains the results of the request on the PDB, the 4000-ish PDB IDs.
- `TESTSET*` are sets from the PISCES server with various PDB IDs, trying to represent the diversity of possible structures.
- `tempresults*` are the calculations of descriptors on the training sets, and `calculations*` on the test set.

`align.py` allows the alignment of DSSP sequences to see if the different chains are structurally close. It uses a simple Needleman and Wunsch algorithm.

`api_informations.py` downloads various information on the structure, sequence, and method used for resolution from the RCSB.

`download*.py` allows the download of either the PDB structure or the assembly directly when the structure is not from X-ray crystallography.

`Coordinates_editor.py` uses the `pdb_numpy` library to reconstruct the structure of X-ray crystallography-issued PDBs in 3D.

`hydro_calc.py` uses simple hydropathy and hydrophobicity scales to do the calculations.

`local_DSSP.py` is a modified version of `pdb_numpy` that allows skipping one calculation of the hydrogen bond map.

`calculations.py` does all the job of synchronizing one by one all the calculations.

`Single_pdb_analysis` allows getting a profile for the execution of a structure, the function that takes time, useful information.


