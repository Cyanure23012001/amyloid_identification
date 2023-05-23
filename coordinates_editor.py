def cooredit(coor):
    try:
        local_rmsd = []
        ref = coor.select_atoms("protein")
        ref_seq_dict = ref.get_aa_seq()

        if len(ref_seq_dict) <= 6:
            ref = ref.add_symmetry()
            ref.compute_chains_CA(Ca_cutoff=4.2)
            ref_seq_dict = ref.get_aa_seq()
            ref = ref.apply_transformation()
            ref.compute_chains_CA(Ca_cutoff=4.2)
            ref_seq_dict = ref.get_aa_seq()
            ref = ref.remove_overlap_chain()
            ref_seq_dict = ref.get_aa_seq()
            ref.compute_chains_CA(Ca_cutoff=4.2)

        if len(ref_seq_dict) <= 6:
            ref = ref.copy_box(x=2, y=2, z=2)
            ref.compute_chains_CA(Ca_cutoff=4.2)
            ref = ref.remove_overlap_chain()
            ref_seq_dict = ref.get_aa_seq()
            ref.compute_chains_CA(Ca_cutoff=4.2)
        print("using REF")
        return ref
    except:
        print("USING COOR")
        return coor
