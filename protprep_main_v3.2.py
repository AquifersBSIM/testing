from protprep_packages.pdb_processor import PDBProcessor
from protprep_packages import sbatch_manager
from protprep_packages.pdb_combiner import PDBCombiner
import pymol
from pymol import cmd
from pathlib import Path
import glob
import re
import os
import pyfiglet

"""
The main script that starts the machine and deep learning
Created by: BSIM
Contact via email: bsim@swinburne.edu.my OR 105572795@students.swinburne.edu.my
"""

def extract_single_chain_and_clean(pdb_file, chain_id, pdb_id, cofactors=None, cosubstrates=None, metal_ions=None):
    # Ensure pdb_file is a Path object.
    pdb_file = Path(pdb_file)
    
    # Convert cofactors, cosubstrates, and metal_ions to lists if a single string is provided.
    if cofactors is not None and isinstance(cofactors, str):
        cofactors = [cofactors]
    if cosubstrates is not None and isinstance(cosubstrates, str):
        cosubstrates = [cosubstrates]
    if metal_ions is not None and isinstance(metal_ions, str):
        metal_ions = [metal_ions]
    
    # Open the original PDB file and read all lines.
    with pdb_file.open('r') as infile:
        lines = infile.readlines()

    filtered_lines = []
    for line in lines:
        # Process ATOM records: include only if they belong to the desired chain.
        if line.startswith("ATOM"):
            if len(line) > 21 and line[21].strip() == chain_id:
                filtered_lines.append(line)
        # Process HETATM records: include only if the residue is in one of the allowed lists
        # and belongs to the desired chain.
        elif line.startswith("HETATM"):
            # Ensure the HETATM record belongs to the desired chain.
            if len(line) > 21 and line[21].strip() != chain_id:
                continue
            resname = line[17:20].strip()  # Residue name is typically in columns 18-20.
            if ((cofactors is not None and resname in cofactors) or
                (cosubstrates is not None and resname in cosubstrates) or
                (metal_ions is not None and resname in metal_ions)):
                filtered_lines.append(line)
        # Preserve all other lines (headers, remarks, etc.)
        else:
            filtered_lines.append(line)
    
    # Create a new file name for the cleaned file.
    new_pdb_file = pdb_file.with_name("rec_" + pdb_file.stem + f"_{chain_id}_clean{pdb_file.suffix}")
    
    # Write the filtered lines to the new PDB file.
    with new_pdb_file.open('w') as outfile:
        outfile.writelines(filtered_lines)
    
    return new_pdb_file
    
def process_pdb():
    processor = PDBProcessor()

    # Get user input and download PDB files/ligands
    pdb_file_dictionary = processor.get_basic_info()
    processor.logger.info("PDB entries collected: %s", list(pdb_file_dictionary.keys()))

    # Process each entry individually and store results in processor.processed_results
    processor.process_all_entries()  # This should loop internally over all entries

    # Now iterate over each processed entry
    for pdb_id, result in processor.processed_results.items():
        pdb_file = pdb_file_dictionary[pdb_id]['pdb_file']
        ligand_to_keep = pdb_file_dictionary[pdb_id]['ligand']
        chain_id = result['chain_id']
        cofactors = result['cofactors']
        cosubstrates = result['cosubstrates']
        metal_ions = result['metal_ions']

        processor.logger.info("Processing PDB: %s", pdb_id)
        processor.logger.info("File downloaded to: %s", pdb_file)

        # Extract and clean the ligand for this pdb entry
        processor.logger.info(f"Extracting the ligand {ligand_to_keep} now")
        sbatch_manager.create_and_run_sbatch_script_extract_and_clean_specific_ligands(pdb_id, pdb_file, ligand_to_keep, chain_id)
        #ligand_pdb = extract_and_clean_specific_ligands(pdb_file, pdb_id, ligand_to_keep, chain_id)
        #processor.logger.info("Ligand file saved to: %s", ligand_pdb)

        # Extract a single chain and clean based on cofactors, cosubstrates, and metal ions
        new_pdb_file = extract_single_chain_and_clean(pdb_file, chain_id, pdb_id, cofactors, cosubstrates, metal_ions)
        
        processor.logger.info(f"Cleaning and adding hydrogens to the receptor protein now")
        # Run the sbatch script for this pdb entry
        sbatch_manager.create_and_run_sbatch_script_add_h(pdb_id, new_pdb_file)

        # FANG ZHE HAW PART !

def combine_pdb():
    pdb_combiner = PDBCombiner()
    pdb_combiner.combine_pdb()

def print_menu():

    width = 85  # Width of the box
    width_1 = 100
    # Create border characters
    top_border = "┏" + "━" * width + "┓"
    middle_border = "┣" + "━" * width + "┫"
    bottom_border = "┗" + "━" * width + "┛"
    
    # Print the options box
    print(top_border)

    # Choose which logo to print (original or fancy)
    
    print("┃" + " ____            _        ____".center(width-15) + "┃".rjust(width-69))
    print("┃" + "|  _ \ _ __ ___ | |_     |  _ \ _ __ ___ _ __".center(width) + "┃")
    print("┃" + "| |_) | '__/ _ \| __|____| |_) | '__/ _ \ '_ \\".center(width) + "┃")
    print("┃" + "  |  __/| | | (_) | ||_____|  __/| | |  __/ |_) |".center(width) + "┃".rjust(width-85))
    print("┃" + "|_|   |_|  \___/ \__|    |_|   |_|  \___| .__/".center(width) + "┃")
    print("┃" + "                                        |_|".center(width-2) + "┃".rjust(width-82))

    print("┃" + "Protein Preparation Platform".center(width) + "┃")
    print("┃" + "Version 3.2".center(width) + "┃")

    print(middle_border)
    print(f"┃ {'Modes:'}{''.ljust(width - 7)}┃")
    
    # List of options
    options = [
        "1. Download and clean receptor protein and its corresponding inhibitory ligands",
        "2. Combine a cleaned receptor protein with a ligand"
    ]
    
    # Print options
    for option in options:
        print(f"┃ {option}{''.ljust(width - len(option) - 1)}┃")
    
    print(bottom_border)

def main():

    print_menu()

    mode_selection = input(f"Enter the mode to use: ")

    if mode_selection == "1":
        process_pdb()
    elif mode_selection == "2":
        combine_pdb()

if __name__ == "__main__":
    main()


