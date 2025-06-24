from pathlib import Path
import pymol
from pymol import cmd
import argparse

def add_hydrogens_water_and_clean_with_pymol(pdb_file, water_radius=5.0):
    pdb_file = Path(pdb_file)
    
    # Launch PyMOL in quiet/headless mode.
    pymol.finish_launching(['pymol', '-qc'])
    
    # Load the structure into PyMOL.
    obj_name = "structure"
    cmd.load(str(pdb_file), obj_name)
    
    # First remove existing solvent
    cmd.remove("solvent")
    
    # Add hydrogens based on standard chemical rules
    cmd.h_add(obj_name)
    
    # Add water molecules around the protein
    # The "water_radius" parameter determines how far from the protein to add water
    cmd.solvent(obj_name, radius=water_radius)
    
    # Save the modified structure back to the same file
    cmd.save(str(pdb_file), obj_name)
    
    print(f"Hydrogens and water added and saved to {pdb_file}")
    
    # Quit PyMOL
    cmd.quit()
    
    return pdb_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add hydrogens and water to a PDB file using PyMOL."
    )

    parser.add_argument(
        "--pdb_file",
        required=True,
        help="Path to the PDB file."
    )
    
    parser.add_argument(
        "--water_radius",
        type=float,
        default=5.0,
        help="Radius around the protein to add water molecules (default: 5.0 Å)."
    )
    
    args = parser.parse_args()
    add_hydrogens_water_and_clean_with_pymol(args.pdb_file, args.water_radius)