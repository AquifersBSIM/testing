import pymol
from pymol import cmd
from openbabel import openbabel
import os

class PymolCombiner:
    def __init__(self, receptor_file, ligand_files):
        self.receptor_file = receptor_file
        self.ligand_files = ligand_files
        pymol.finish_launching(['pymol', '-cq'])
    
    def convert_pdbqt_to_pdb(self, pdbqt_file):
        """
        Convert PDBQT files to PDB format using OpenBabel
        
        Parameters:
            pdbqt_file (str): Path to the PDBQT file
        
        Returns:
            str: Path to the converted PDB file
        """
        pdb_file = os.path.splitext(pdbqt_file)[0] + ".pdb"
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdbqt", "pdb")
        mol = openbabel.OBMol()
        print(f"\nConverting {pdbqt_file} to PDB format...")
        obConversion.ReadFile(mol, pdbqt_file)
        obConversion.WriteFile(mol, pdb_file)
        print(f"Conversion complete: {pdb_file}\n")
        return pdb_file
    
    def load_structures(self):
        # Handle receptor files
        self.receptor_pdbs = []
        if isinstance(self.receptor_file, list):
            for receptor in self.receptor_file:
                if isinstance(receptor, str) and receptor.endswith('.pdbqt'):
                    self.receptor_pdbs.append(self.convert_pdbqt_to_pdb(receptor))
                else:
                    self.receptor_pdbs.append(receptor)
        else:
            if isinstance(self.receptor_file, str) and self.receptor_file.endswith('.pdbqt'):
                self.receptor_pdbs = [self.convert_pdbqt_to_pdb(self.receptor_file)]
            else:
                self.receptor_pdbs = [self.receptor_file]
        
        self.receptor_names = []
        for i, receptor_pdb in enumerate(self.receptor_pdbs, start=1):
            receptor_name = "receptor" if i == 1 else f"receptor_{i}"
            self.receptor_names.append(receptor_name)
            print(f"Loading receptor protein: {receptor_pdb} as {receptor_name}")
            cmd.load(receptor_pdb, receptor_name)
        
        # Process ligand files
        self.converted_ligands = []
        self.ligand_base_names = []
        if not isinstance(self.ligand_files, list):
            self.ligand_files = [self.ligand_files]
            
        for ligand in self.ligand_files:
            if isinstance(ligand, str) and ligand.endswith('.pdbqt'):
                converted_ligand = self.convert_pdbqt_to_pdb(ligand)
                self.converted_ligands.append(converted_ligand)
            else:
                self.converted_ligands.append(ligand)
            base_name = os.path.splitext(os.path.basename(ligand))[0]
            self.ligand_base_names.append(base_name)
        
        self.ligand_object_names = []
        for i, ligand in enumerate(self.converted_ligands, start=1):
            object_name = f"ligand_{i}"
            self.ligand_object_names.append(object_name)
            print(f"Loading ligand file: {ligand} as {object_name}")
            cmd.load(ligand, object_name)
    
    def combine_structures(self, output_prefix=None, individual_complexes=True):
        """
        Combines loaded receptor(s) and ligand(s) into one or more PDB files.
        If individual_complexes is True, each ligand will have its own output file combined with the receptor(s).
        Otherwise, all structures are combined into one file.
        """
        all_objects = cmd.get_object_list()
        if not all_objects:
            print("Error: No objects loaded to combine\n")
            return
        output_files = []
        if individual_complexes:
            for i, ligand_name in enumerate(self.ligand_object_names):
                if output_prefix:
                    output_file = f"{output_prefix}_{self.ligand_base_names[i]}.pdb"
                else:
                    output_file = f"complex_{self.ligand_base_names[i]}.pdb"
                complex_objects = self.receptor_names + [ligand_name]
                cmd.select("current_complex", " or ".join(complex_objects))
                complex_obj_name = f"complex_{i+1}"
                cmd.create(complex_obj_name, "current_complex")
                cmd.save(output_file, complex_obj_name)
                print(f"Complex saved as {output_file}\n")
                output_files.append(output_file)
        else:
            if output_prefix is None:
                if len(self.ligand_base_names) == 1:
                    output_file = f"complex_{self.ligand_base_names[0]}.pdb"
                else:
                    output_file = f"complex_{'_and_'.join(self.ligand_base_names[:2])}"
                    if len(self.ligand_base_names) > 2:
                        output_file += f"_plus_{len(self.ligand_base_names)-2}_more"
                    output_file += ".pdb"
            else:
                output_file = f"{output_prefix}.pdb"
            cmd.select("all_objects", " or ".join(all_objects))
            cmd.create("combined", "all_objects")
            cmd.save(output_file, "combined")
            print(f"Combined structure saved as {output_file}\n")
            output_files.append(output_file)
        return output_files
    
    def extract_ligand_metadata(self, ligand_file):
        keywords = ["MODEL", "REMARK  Name =", "REMARK VINA RESULT:", 
                    "REMARK INTER + INTRA:", "REMARK INTER:", "REMARK INTRA:", "REMARK UNBOUND:"]
        metadata = []
        model1_started = False
        try:
            with open(ligand_file, "r") as f:
                for line in f:
                    line_stripped = line.strip()
                    # Check for model header
                    if line_stripped.startswith("MODEL"):
                        # If we haven't started yet and it's MODEL 1, mark as started and add the converted line.
                        if not model1_started and line_stripped == "MODEL 1":
                            metadata.append("REMARK " + line_stripped)
                            model1_started = True
                            continue
                        # If a new model header is found and we already collected MODEL 1, break out.
                        if model1_started and line_stripped != "MODEL 1":
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

    def insert_metadata_into_pdb_content(self, pdb_file, metadata_lines):
        """
        Inserts metadata lines at the beginning of the PDB file.
        
        Parameters:
            pdb_file (str): Path to the PDB file to modify.
            metadata_lines (list): List of metadata lines to insert.
        """
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
    
    def combine_structures_with_metadata(self, output_prefix=None, individual_complexes=True):
        """
        Combines the structures as in combine_structures(), then extracts metadata from the
        original ligand file(s) and inserts those lines into the beginning of the output PDB file(s).
        
        In 'individual_complexes' mode, each output file will receive metadata from its corresponding ligand file.
        In combined mode, metadata from all ligand files are concatenated (duplicates removed) and inserted.
        
        Returns:
            list: List of output PDB file paths.
        """
        output_files = self.combine_structures(output_prefix=output_prefix, individual_complexes=individual_complexes)
        if not output_files:
            print("No output files to process for metadata insertion.")
            return None
        
        if individual_complexes:
            # For each complex file, use the corresponding ligand file for metadata.
            for i, output_file in enumerate(output_files):
                ligand_file = self.ligand_files[i]
                metadata = self.extract_ligand_metadata(ligand_file)
                self.insert_metadata_into_pdb_content(output_file, metadata)
        else:
            # For a single combined file, merge metadata from all ligand files.
            combined_metadata = []
            for ligand_file in self.ligand_files:
                metadata = self.extract_ligand_metadata(ligand_file)
                combined_metadata.extend(metadata)
            # Remove any duplicate metadata lines while preserving order
            unique_metadata = list(dict.fromkeys(combined_metadata))
            self.insert_metadata_into_pdb_content(output_files[0], unique_metadata)
        return output_files
