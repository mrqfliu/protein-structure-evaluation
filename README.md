# Protein Structure Evaluation Toolkit  
A Python toolkit for evaluating structural similarity between predicted and experimental protein complexes, with support for ligand-binding site analysis.  


## Overview  
This toolkit provides a comprehensive workflow to compare predicted protein structures (e.g., from AlphaFold, RoseTTAFold) with experimental structures (from PDB). It includes structural alignment algorithms and multiple evaluation metrics to assess global and local structural accuracy, with special focus on ligand-protein interactions.  


## Key Features  
- **Structural Alignment**: Supports both RMSD-based and TM-align alignment (Cα atom-focused).  
- **Evaluation Metrics**:  
  - RMSD (Root Mean Square Deviation) for global and ligand-specific accuracy.  
  - TM-score (Template Modeling score) for assessing fold similarity.  
  - LDDT (Local Distance Difference Test) for local structural quality.  
  - PAE (Predicted Aligned Error) for residue-level distance errors.  
- **Ligand-Specific Analysis**: Separately evaluates ligand-binding site accuracy by distinguishing hetero chains (ligands) from protein chains.  


## Dependencies  
- Python 3.8+  
- Required libraries:  
  ```bash
  biopython>=1.81  # For PDB parsing
  tmtools>=0.5.0   # For TM-align algorithm
  numpy>=1.24.0    # For numerical calculations
  networkx>=3.1    # For graph-based structure analysis (optional)
  ```  


## Installation  
1. Create the environment
   ```bash
   conda create -n protein-eval python=3.8 -y
   conda activate protein-eval
   ```
   
2. Clone the repository:  
   ```bash
   git clone https://github.com/your-username/protein-structure-evaluation.git
   cd protein-structure-evaluation
   ```  

3. Install dependencies:  
   ```bash
   pip install -r requirements.txt
   ```  

   *(If `tmtools` is unavailable via PyPI, install from source: `pip install git+https://github.com/your-username/tmtools.git`)*  


## Usage  
### Input Files  
- Experimental structure: PDB file of the experimentally determined protein complex (e.g., from RCSB PDB).  
- Predicted structure: PDB file of the computationally predicted complex.  


### Basic Workflow  
1. Prepare input PDB files and place them in `./input_files/from_rcsb/` (or modify paths in the script).  

2. Run the evaluation script:  
   ```bash
   python evaluate_structures.py
   ```  


### Customization  
Modify the following parameters in `evaluate_structures.py` to fit your needs:  
- `experimental_complex_pdb`: Path to the experimental PDB file.  
- `predicted_complex_pdb`: Path to the predicted PDB file.  
- Alignment method: Toggle between `align_structures_rmsd()` and `align_structures_tm()`.  
- LDDT tolerances: Adjust `tolerances=[0.5, 1, 2, 4]` in `calculate_lddt()` for local quality assessment.  


## Output  
The toolkit generates:  
1. An aligned PDB file (`aligned_predicted_complex.pdb`) for visual inspection (e.g., in PyMOL).  
2. Console output with key metrics:  
   - Global and ligand-specific RMSD, TM-score, LDDT, and PAE.  
   - Rotation matrices and translation vectors from alignment steps.  


## Example Output  
```  
Rotation matrix:
[[ 0.998  0.021 -0.054]
 [-0.020  0.999  0.012]
 [ 0.054 -0.010  0.998]]
Translation vector:
[1.234, -0.567, 2.345]
The tmalign-CA-TM-score between the two structures is: 0.892
The rmsdalign-CA-RMSD between the two structures is: 1.23 Å

Number of matching atoms: 234
Number of matching ligand atoms: 12

The overall LDDT between predicted and experimental complex structures is: 0.87
The overall LDDT for ligands between predicted and experimental complex structures is: 0.76
The overall TM-score between predicted and experimental complex structures is: 0.89
The overall TM-score for ligands between predicted and experimental complex structures is: 0.72
The overall RMSD between predicted and experimental complex structures is: 1.56 Å
The overall RMSD for ligands between predicted and experimental complex structures is: 0.89 Å
```  


## Notes  
- **PDB Format**: Ensure input PDB files include standard residue names (e.g., "ALA", "GLU") and hetero chains (ligands) are labeled correctly (e.g., "ATP", "HEM").  
- **Ligand Detection**: Hetero chains (ligands) are automatically identified by checking for non-standard amino acids.  
- **Alignment Order**: The script first runs RMSD alignment, then refines with TM-align (modify the order in `__main__` if needed).  


## License  
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.  


## Citation  
If you use this toolkit in your research, please cite:  
*(Add your citation here if published)*  


---  
Contact: [your-email@example.com] for issues or feature requests.
