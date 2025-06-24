from pathlib import Path
import pymol
from pymol import cmd
import argparse

def extract_ligand_metadata(ligand_file):
    keywords = ["MODEL", "REMARK  Name =", "REMARK VINA RESULT:", 
                "REMARK INTER + INTRA:", "REMARK INTER:", "REMARK INTRA:", "REMARK UNBOUND:"]
    metadata = []
    model1_started = False
    try:
        with open(ligand_file, "r") as f:
            for line in f:
                line_stripped = line.strip()
                # Check for model header
                if line_stripped.startswith("REMARK MODEL"):
                    # If we haven't started yet and it's MODEL 1, mark as started and add the converted line.
                    if not model1_started and line_stripped == "REMARK MODEL 1":
                        metadata.append(line_stripped)
                        model1_started = True
                        continue
                    # If a new model header is found and we already collected MODEL 1, break out.
                    if model1_started and line_stripped != "REMARK MODEL 1":
                        break

                # Only process lines if we're in the MODEL 1 block.
                if model1_started:
                    for key in keywords:
                        if line_stripped.startswith(key) and key != "MODEL":
                            metadata.append(line_stripped)
                            break
    except Exception as e:
        print(f"Error reading ligand file {ligand_file}: {e}")
    
    return metadata

def insert_metadata_into_pdb_content(pdb_file, metadata_lines):
    try:
        with open(pdb_file, "r") as f:
            original_content = f.read()
        with open(pdb_file, "w") as f:
            for line in metadata_lines:
                f.write(line + "\n")
            f.write(original_content)
        print(f"Inserted metadata into {pdb_file}")
    except Exception as e:
        print(f"Error inserting metadata into pdb file {pdb_file}: {e}")

def add_hydrogens_and_clean_with_pymol(pdb_file):
    pdb_file = Path(pdb_file)
    
    # Launch PyMOL in quiet/headless mode.
    pymol.finish_launching(['pymol', '-qc'])
    
    # Load the structure into PyMOL.
    obj_name = "structure"
    cmd.load(str(pdb_file), obj_name)
    cmd.remove("solvent")
    
    # Add hydrogens based on standard chemical rules.
    cmd.h_add(obj_name)
    
    # Remove nonpolar hydrogens.
    cmd.remove(f"{obj_name} & hydro & not nbr. (don.|acc.)")
    
    # Save the modified structure back to the same file.
    cmd.save(str(pdb_file), obj_name)
    
    print(f"Hydrogens added and nonpolar hydrogens removed. Saved to {pdb_file}")
    
    # Quit PyMOL.
    cmd.quit()
    
    return pdb_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add hydrogens and remove nonpolar hydrogens from a given PDB file."
    )

    parser.add_argument(
        "--pdb_file",
        required=True,
        help="Path to the PDB file."
    )
    
    args = parser.parse_args()

    # Define the metadata keywords to search for.
    metadata_keywords = ["MODEL", "REMARK  Name =", "REMARK VINA RESULT:",
                         "REMARK INTER + INTRA:", "REMARK INTER:", "REMARK INTRA:", "REMARK UNBOUND:"]
    metadata_present = False
    # Check if the PDB file already contains any metadata lines.
    try:
        with open(args.pdb_file, "r") as file:
            for line in file:
                if any(line.strip().startswith(key) for key in metadata_keywords):
                    print(f"Metadata present !")
                    metadata_present = True
                    break
    except Exception as e:
        print(f"Error checking metadata in {args.pdb_file}: {e}")

    # If metadata is present, extract it.
    if metadata_present:
        metadata = extract_ligand_metadata(args.pdb_file)
    else:
        metadata = []

    # Process the file with PyMOL.
    pdb_file = add_hydrogens_and_clean_with_pymol(args.pdb_file)

    # If we extracted metadata, insert it at the top of the new file.
    if metadata:
        insert_metadata_into_pdb_content(pdb_file, metadata)
