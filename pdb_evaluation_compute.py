from Bio.PDB import PDBParser, PDBIO
from Bio.PDB import Chain
import numpy as np
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data, get_atom_data
from Bio.PDB import PDBParser, Superimposer

# List of standard amino acids
standard_amino_acids = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
]

# Define the is_hetero method to check if a chain is hetero (non-standard amino acid chain)
def is_hetero(chain: Chain) -> bool:
    for residue in chain:
        # Check if the residue name is a standard amino acid
        if residue.get_resname() in standard_amino_acids:
            return False  # Not a hetero chain if any standard amino acid exists
    return True  # hetero chain if no standard amino acids found


# TM-align: performs better than RMSD alignment for short distances, worse for long distances
def align_structures_tm(experimental_complex, predicted_complex):
    # TM alignment based on Cα atoms of the receptor
    # Get the first chain from each structure
    experimental_chain = next(experimental_complex.get_chains())
    predicted_chain = next(predicted_complex.get_chains())

    # Get residue data (Cα coordinates and sequences) for the chains
    coords1, seq1 = get_residue_data(experimental_chain)
    coords2, seq2 = get_residue_data(predicted_chain)

    # Perform alignment using TM-align algorithm
    result = tm_align(coords1, coords2, seq1, seq2)
    
    # Get coordinates of all atoms
    all_predicted_coords = np.array([atom.get_coord() for atom in predicted_complex.get_atoms()])
    
    # Apply rotation and translation
    rotated_coords = np.dot(all_predicted_coords, result.u) + result.t
    
    # Update atom coordinates in the predicted_complex
    for atom, new_coord in zip(predicted_complex.get_atoms(), rotated_coords):
        atom.coord = new_coord
    
    # Print rotation matrix and translation vector
    print("Rotation matrix:")
    print(result.u)
    print("Translation vector:")
    print(result.t)
    
    # Get TM-score and RMSD (Cα only)
    tm = result.tm_norm_chain2
    rms = result.rmsd
    print(f"The tmalign-CA-TM-score between the two structures is: {tm}")
    print(f"The tmalign-CA-RMSD between the two structures is: {rms} Å")

    # Create PDBIO object to write output file
    io = PDBIO()
    io.set_structure(predicted_complex)
    # Save the aligned structure to a new PDB file
    io.save("aligned_predicted_complex.pdb")
    return predicted_complex


# RMSD alignment: suitable for all cases
def align_structures_rmsd(experimental_complex, predicted_complex):
    # Create Superimposer object
    sup = Superimposer()
    
    # Get all Cα atoms from experimental and predicted structures
    experimental_atoms_ca = [atom for atom in experimental_complex.get_atoms() if atom.get_name() == 'CA']
    predicted_atoms_ca = [atom for atom in predicted_complex.get_atoms() if atom.get_name() == 'CA']
    predicted_atoms = list(predicted_complex.get_atoms())
    
    # Set atom groups for alignment
    sup.set_atoms(experimental_atoms_ca, predicted_atoms_ca)  # (fixed, moving)
    
    # Perform alignment and update coordinates
    sup.apply(predicted_atoms)
    
    # Get rotation matrix and translation vector
    rot_matrix, tran_vector = sup.rotran
    
    # Print rotation matrix and translation vector
    print("Rotation matrix:")
    print(rot_matrix)
    print("Translation vector:")
    print(tran_vector)

    # Get RMSD value (Cα only)
    rmsd = sup.rms
    print(f"The rmsdalign-CA-RMSD between the two structures is: {rmsd} Å")
    
    # Create PDBIO object to write output file
    io = PDBIO()
    io.set_structure(predicted_complex)
    # Save the aligned structure to a new PDB file
    io.save("aligned_predicted_complex.pdb")
    return predicted_complex


def create_graph(atom_list, atom_ids):
    import networkx as nx

    G = nx.Graph()

    for i, atom_i in enumerate(atom_list):
        cr_i = COVALENT_RADIUS[atom_ids[i]]
        for j, atom_j in enumerate(atom_list):
            cr_j = COVALENT_RADIUS[atom_ids[j]]
            distance = np.linalg.norm(atom_i - atom_j)
            threshold = (cr_i + cr_j + BOND_TOLERANCE) if i != j else 1
            if distance < threshold:  # Adjust threshold as needed
                G.add_edge(i, j)

    return G


def calculate_pae_rmsd(experimental_complex, predicted_complex):
    # Create mappings to search atoms by residue name, atom name, and residue index
    predicted_atom_map = {}
    experimental_atom_map = {}
    predicted_ligand_atom_map = {}
    experimental_ligand_atom_map = {}    
    
    # Define ligand chains
    predicted_ligand_chains = [chain for chain in predicted_complex.get_chains() if is_hetero(chain)]
    # Define experimental structure ligand chains
    experimental_ligand_chains = [chain for chain in experimental_complex.get_chains() if is_hetero(chain)]
    
    # Build atom mappings
    for chain in predicted_complex.get_chains():
        residue_id = 0  # Avoid issues with different starting indices in structures
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                predicted_atom_map[key] = atom
            residue_id += 1
    
    # Mapping for experimental structure
    for chain in experimental_complex.get_chains():
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                experimental_atom_map[key] = atom
            residue_id += 1
    
    # Build atom mapping for predicted structure ligands
    for chain in predicted_ligand_chains:
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                predicted_ligand_atom_map[key] = atom
            residue_id += 1
    
    # Build atom mapping for experimental structure ligands
    for chain in experimental_ligand_chains:
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                experimental_ligand_atom_map[key] = atom
            residue_id += 1
    
    # Ensure atom counts match between predicted and experimental structures
    predicted_keys = set(predicted_atom_map.keys()).intersection(set(experimental_atom_map.keys()))
    
    # Get list of matching ligand atom keys in predicted structure
    predicted_ligand_keys = set(predicted_ligand_atom_map.keys()).intersection(set(experimental_ligand_atom_map.keys()))
    
    # Initialize lists for PAE and RMSD
    pae_list = []
    rmsd_list = []
    # Lists for ligand PAE and RMSD
    ligand_pae_list = []
    ligand_rmsd_list = []
    
    # Calculate ligand PAE and RMSD (ensure all four key components match)
    for key in predicted_ligand_keys:
        if key in experimental_ligand_atom_map:
            pred_atom = predicted_ligand_atom_map[key]
            exp_atom = experimental_ligand_atom_map[key]
            distance = np.linalg.norm(pred_atom.get_coord() - exp_atom.get_coord())
            ligand_rmsd_list.append(distance**2)
            ligand_pae_list.append(distance)
    
    # Calculate overall PAE and RMSD
    for key in predicted_keys:
        if key in experimental_atom_map:
            pred_atom = predicted_atom_map[key]
            exp_atom = experimental_atom_map[key]
            distance = np.linalg.norm(pred_atom.get_coord() - exp_atom.get_coord())
            rmsd_list.append(distance**2)
            pae_list.append(distance)
    
    # Calculate average ligand PAE
    if ligand_pae_list:
        ligand_average_pae = np.mean(ligand_pae_list)
    else:
        ligand_average_pae = None
    
    # Calculate average ligand RMSD
    if ligand_rmsd_list:
        ligand_average_rmsd = np.sqrt(np.mean(ligand_rmsd_list))
    else:
        ligand_average_rmsd = None
    
    # Calculate average PAE
    if pae_list:  # Ensure the list is not empty
        average_pae = np.mean(pae_list)
    else:
        average_pae = None  # Return None if no matching atoms
    
    # Calculate average RMSD
    if rmsd_list:  # Ensure the list is not empty
        average_rmsd = np.sqrt(np.mean(rmsd_list))
    else:
        average_rmsd = None  # Return None if no matching atoms
    
    return ligand_average_pae, average_pae, ligand_average_rmsd, average_rmsd


def calculate_tm_score(experimental_complex, predicted_complex):
    # Create mappings to search atoms by residue name, atom name, and residue index
    predicted_atom_map = {}
    experimental_atom_map = {}
    predicted_ligand_atom_map = {}
    experimental_ligand_atom_map = {}
    
    # Build atom mappings
    for chain in predicted_complex.get_chains():
        residue_id = 0  # Avoid issues with different starting indices in structures
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                predicted_atom_map[key] = atom
            residue_id += 1
    
    # Mapping for experimental structure
    for chain in experimental_complex.get_chains():
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                experimental_atom_map[key] = atom
            residue_id += 1
    
    # Define ligand chains in predicted structure
    predicted_ligand_chains = [chain for chain in predicted_complex.get_chains() if is_hetero(chain)]
    # Define ligand chains in experimental structure
    experimental_ligand_chains = [chain for chain in experimental_complex.get_chains() if is_hetero(chain)]

    # Build atom mapping for predicted structure ligands
    for chain in predicted_ligand_chains:
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                predicted_ligand_atom_map[key] = atom
            residue_id += 1
    
    # Build atom mapping for experimental structure ligands
    for chain in experimental_ligand_chains:
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                experimental_ligand_atom_map[key] = atom
            residue_id += 1
    
    # Get list of matching atom keys in predicted structure
    predicted_keys = set(predicted_atom_map.keys()).intersection(set(experimental_atom_map.keys()))
    print(f"Number of matching atoms: {len(list(predicted_keys))}")
    
    # Get list of matching ligand atom keys in predicted structure
    predicted_ligand_keys = set(predicted_ligand_atom_map.keys()).intersection(set(experimental_ligand_atom_map.keys()))
    print(f"Number of matching ligand atoms: {len(list(predicted_ligand_keys))}")
    
    # Initialize lists for TM-score calculation
    tm_list = []
    # Initialize list for ligand TM-score
    tm_ligand_list = []
    
    # Calculate overall TM-score
    for key in predicted_keys:
        if key in experimental_atom_map:
            pred_atom = predicted_atom_map[key]
            exp_atom = experimental_atom_map[key]
            distance = np.linalg.norm(pred_atom.get_coord() - exp_atom.get_coord())
            tm_list.append(distance)
    
    # Calculate number of atoms in experimental structure
    L_N = len(list(experimental_complex.get_atoms()))
    d0 = max(4.5, min(8, 1.24 * (L_N - 15) ** (1/3) - 1.8))

    tm_score = 1 / L_N * sum(1 / (1 + (d / d0) ** 2) for d in tm_list)
    
    # Calculate ligand TM-score
    for key in predicted_ligand_keys:
        if key in experimental_ligand_atom_map:
            pred_atom = predicted_ligand_atom_map[key]
            exp_atom = experimental_ligand_atom_map[key]
            distance = np.linalg.norm(pred_atom.get_coord() - exp_atom.get_coord())
            tm_ligand_list.append(distance)
    
    # Calculate number of ligand atoms in experimental structure
    lig_L_N = sum(len(list(chain.get_atoms())) for chain in experimental_ligand_chains)
    lig_d0 = max(4.5, min(8, 0.6 * (lig_L_N - 0.5) ** (1/2) - 2.5))
    
    # Calculate ligand TM-score
    ligand_tm_score = 1 / lig_L_N * sum(1 / (1 + (d / lig_d0) ** 2) for d in tm_ligand_list)

    return tm_score, ligand_tm_score


def calculate_lddt(experimental_complex, predicted_complex, r=5, tolerances=[0.5, 1, 2, 4]):
    # Create mappings
    experimental_atom_map = {}
    predicted_atom_map = {}
    predicted_ligand_atom_map = {}
    experimental_ligand_atom_map = {}
    experimental_protein_atom_map = {}
    predicted_protein_atom_map = {}
    
    # Build atom mapping for experimental structure
    for chain in experimental_complex.get_chains():
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                experimental_atom_map[key] = atom
                if is_hetero(chain):
                    experimental_ligand_atom_map[key] = atom
                else:
                    experimental_protein_atom_map[key] = atom
            residue_id += 1
    
    # Build atom mapping for predicted structure
    for chain in predicted_complex.get_chains():
        residue_id = 0
        for residue in chain:
            for atom in residue:
                key = (chain.get_id(), residue_id, residue.get_resname(), atom.get_name())
                predicted_atom_map[key] = atom
                if is_hetero(chain):
                    predicted_ligand_atom_map[key] = atom
                else:
                    predicted_protein_atom_map[key] = atom
            residue_id += 1

    # Select atom pairs for all atoms
    selected_pairs = []
    exp_keys = list(experimental_atom_map.keys())
    for i in range(len(exp_keys)):
        for j in range(i + 1, len(exp_keys)):
            key1 = exp_keys[i]
            key2 = exp_keys[j]
            # Ensure residues are different and sequence separation is at least r
            if abs(key1[1] - key2[1]) >= r:
                atom1 = experimental_atom_map[key1]
                atom2 = experimental_atom_map[key2]
                distance_exp = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
                if distance_exp < 15:
                    selected_pairs.append((key1, key2))
                    
    # Select protein-ligand atom pairs
    protein_ligand_pairs = []
    for exp_lig_key, exp_lig_atom in experimental_ligand_atom_map.items():
        for exp_prot_key, exp_prot_atom in experimental_protein_atom_map.items():
            distance_exp = np.linalg.norm(exp_lig_atom.get_coord() - exp_prot_atom.get_coord())
            if distance_exp < 15:
                protein_ligand_pairs.append((exp_lig_key, exp_prot_key))
    
    # Calculate distance differences for ligands
    ligand_distance_diffs = []
    for (exp_lig_key, exp_prot_key) in protein_ligand_pairs:
        exp_lig_atom = experimental_ligand_atom_map[exp_lig_key]
        exp_prot_atom = experimental_protein_atom_map[exp_prot_key]
        pred_lig_atom = predicted_ligand_atom_map.get(exp_lig_key)
        pred_prot_atom = predicted_protein_atom_map.get(exp_prot_key)
        
        if pred_lig_atom is None or pred_prot_atom is None:
            continue
        
        distance_exp = np.linalg.norm(exp_lig_atom.get_coord() - exp_prot_atom.get_coord())
        distance_pred = np.linalg.norm(pred_lig_atom.get_coord() - pred_prot_atom.get_coord())
        distance_diff = abs(distance_exp - distance_pred)
        ligand_distance_diffs.append(distance_diff)
    
    # Apply tolerances and calculate ligand LDDT
    ligand_lddt_scores = []
    for tolerance in tolerances:
        within_tolerance = sum(1 for diff in ligand_distance_diffs if diff < tolerance)
        lddt_score = within_tolerance / len(ligand_distance_diffs) if ligand_distance_diffs else 0
        ligand_lddt_scores.append(lddt_score)
    
    ligand_lddt_average = np.mean(ligand_lddt_scores) if ligand_lddt_scores else 0
    
    # Calculate distance differences for overall structure
    distance_diffs = []
    for (exp_key1, exp_key2) in selected_pairs:
        exp_atom1 = experimental_atom_map[exp_key1]
        exp_atom2 = experimental_atom_map[exp_key2]
        pred_atom1 = predicted_atom_map.get(exp_key1)
        pred_atom2 = predicted_atom_map.get(exp_key2)
        
        if pred_atom1 is None or pred_atom2 is None:
            continue  # Skip missing atoms
        
        distance_exp = np.linalg.norm(exp_atom1.get_coord() - exp_atom2.get_coord())
        distance_pred = np.linalg.norm(pred_atom1.get_coord() - pred_atom2.get_coord())
        distance_diff = abs(distance_exp - distance_pred)
        distance_diffs.append(distance_diff)
    
    # Apply tolerances and calculate LDDT scores
    lddt_scores = []
    for tolerance in tolerances:
        within_tolerance = sum(1 for diff in distance_diffs if diff < tolerance)
        lddt_score = within_tolerance / len(distance_diffs) if distance_diffs else 0
        lddt_scores.append(lddt_score)
    
    # Calculate average LDDT
    lddt_average = np.mean(lddt_scores) if lddt_scores else 0
    
    return lddt_average, ligand_lddt_average


if __name__ == '__main__':
    # Load PDB structures and calculate PAE, RMSD, TM-score, LDDT
    experimental_complex_pdb = "./input_files/from_rcsb/5fui_complex_real.pdb"  # Path to experimental complex PDB file
    predicted_complex_pdb = "./input_files/from_rcsb/5fui_complex_tank.pdb"  # Path to predicted complex PDB file
    
    # Create PDB parser
    parser = PDBParser()
    
    # Load experimental complex structure
    experimental_complex = parser.get_structure("ExperimentalComplex", experimental_complex_pdb)
    
    # Load predicted complex structure
    predicted_complex = parser.get_structure("PredictedComplex", predicted_complex_pdb)
    
    # Align structures using RMSD and TM-align (based on Cα atoms)
    predicted_complex = align_structures_rmsd(experimental_complex, predicted_complex)
    predicted_complex = align_structures_tm(experimental_complex, predicted_complex)
    
    # Calculate PAE and RMSD
    ligand_pae, pae, ligand_rmsd, rmsd = calculate_pae_rmsd(experimental_complex, predicted_complex)
    
    # Calculate TM-score
    tm_score, tm_ligand_score = calculate_tm_score(experimental_complex, predicted_complex)
    
    # Calculate LDDT
    lddt, ligand_lddt = calculate_lddt(experimental_complex, predicted_complex)
    
    # Print results
    print(f"The overall LDDT for ligands between predicted and experimental complex structures is: {ligand_lddt}")
    print(f"The overall LDDT between predicted and experimental complex structures is: {lddt}")
    print(f"The overall TM-score for ligands between predicted and experimental complex structures is: {tm_ligand_score}")
    print(f"The overall TM-score between predicted and experimental complex structures is: {tm_score}")
    
    if ligand_pae is not None:
        print(f"The overall PAE for ligands between predicted and experimental complex structures is: {ligand_pae} Å")
    else:
        print("No matching ligand atoms found for PAE calculation.")
    
    if pae is not None:
        print(f"The overall PAE between predicted and experimental complex structures is: {pae} Å")
    else:
        print("No matching atoms found for PAE calculation.")
    
    if ligand_rmsd is not None:
        print(f"The overall RMSD for ligands between predicted and experimental complex structures is: {ligand_rmsd} Å")
    else:
        print("No matching ligand atoms found for RMSD calculation.")
    
    if rmsd is not None:
        print(f"The overall RMSD between predicted and experimental complex structures is: {rmsd} Å")
    else:
        print("No matching atoms found for RMSD calculation.")
