import pymol
from pymol import cmd
from pathlib import Path
import argparse 

def extract_and_clean_specific_ligands(pdb_file, pdb_id, ligand_to_keep, chain_id):
    # Launch PyMOL (using -cq flags for quiet command-line execution)
    pymol.finish_launching(['pymol', '-cq'])
    
    # Load the PDB file into an object named by pdb_id
    cmd.load(pdb_file, pdb_id)
    
    # Remove solvent molecules to clean the structure
    cmd.remove("solvent")
    
    # Select the ligand to keep on the specified chain.
    # "organic" ensures that we are working with small molecules.
    # "chain" filters by the provided chain identifier.
    # "resn" is used to filter by the residue name of the ligand.
    selection_str = f"organic and chain {chain_id} and resn {ligand_to_keep}"
    cmd.select("ligand_keep", selection_str)
    
    # Add hydrogens to the selected ligand.
    cmd.h_add("all")
    
    # Save the selected ligand to a new PDB file named "lig_{pdb_id}.pdb"
    output_file = f"lig_{pdb_id}.pdb"
    cmd.save(output_file, "ligand_keep")
    
    # Clean up and quit PyMOL (useful in script mode)
    cmd.delete("all")
    cmd.quit()

    return output_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create and run an sbatch script using a given pdb_id and pdb_file."
    )
    parser.add_argument('--pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('--pdb_id', type=str, help='PDB identifier')
    parser.add_argument('--ligand_to_keep', type=str, help='Name of the ligand to keep')
    parser.add_argument('--chain_id', type=str, help='Chain identifier')

    args = parser.parse_args()

    extract_and_clean_specific_ligands(args.pdb_file, args.pdb_id, args.ligand_to_keep, args.chain_id)