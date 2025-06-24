import os
import glob
import re
from .combinator_pdb import PymolCombiner
from protprep_packages import sbatch_manager

class PDBCombiner:
    def __init__(self):
        # Initialization can be expanded as needed
        pass

    def combine_pdb(self):
        # Get receptor files and ligand files
        receptor_protein = glob.glob("rec_*_clean.pdbqt")
        pdb_id = None
        if receptor_protein:
            match = re.search(r"rec_(.+?)_clean\.pdbqt", receptor_protein[0])
            if match:
                pdb_id = match.group(1)

        pdbqt_directory = input(
            "Please input the path to your CSV files (e.g., /path/to/file). "
            "Press Enter to use the current directory:\n"
        ).strip()

        # Use the current directory if the user presses Enter
        if not pdbqt_directory:
            pdbqt_directory = os.getcwd()

        # Get all .pdbqt files and filter out receptor files
        ligand_files = [
            f for f in glob.glob(os.path.join(pdbqt_directory, "*.pdbqt"))
            if not (os.path.basename(f).startswith("rec_") and os.path.basename(f).endswith("_clean.pdbqt"))
        ]
        
        if not ligand_files:
            print("No ligand .pdbqt files found.")
            return

        # Determine the maximum filename length for border formatting
        max_length = max([len(f) for f in ligand_files], default=0)
        width = max(max_length + 10, 30)  # Ensure a minimum width

        # Create border characters
        top_border = "┏" + "━" * width + "┓"
        middle_border = "┣" + "━" * width + "┫"
        bottom_border = "┗" + "━" * width + "┛"

        # Print header with border
        print(top_border)
        print("┃ " + "Ligand(s) found in structure:".ljust(width - 2) + " ┃")
        print(middle_border)

        # Display the ligand options
        for i, ligand in enumerate(ligand_files, start=1):
            print("┃ " + f"{i}. {os.path.basename(ligand)}".ljust(width - 2) + " ┃")

        print(bottom_border)
        
        # Ask user to select a ligand or choose all
        choice = input("Select a ligand by number, or type 'all' to use all ligands: ").strip()

        if choice.lower() == 'all':
            selected_ligands = ligand_files
        else:
            try:
                index = int(choice)
                if 1 <= index <= len(ligand_files):
                    selected_ligands = [ligand_files[index - 1]]
                else:
                    print("Invalid number selection. Exiting.")
                    return
            except ValueError:
                print("Invalid input. Please enter a number or 'all'. Exiting.")
                return

        print("Selected ligands:", selected_ligands)
        
        # Instantiate the PymolCombiner (assumed to be defined elsewhere)
        combiner = PymolCombiner(receptor_protein, selected_ligands)
        combiner.load_structures()
        # You can switch between combine_structures() and combine_structures_with_metadata() as needed
        output_file = combiner.combine_structures_with_metadata(individual_complexes=True)

        # Ask for post processing options
        post_processing_input = input(
            "Now would you like to add polar hydrogens and/or water? Enter (yes/no)\n"
        ).strip().lower()

        if post_processing_input == "yes":
            print("\nOptions:")
            print("1. Add polar hydrogens only")
            #print("2. Add waters only")
            #print("3. Add polar hydrogens and water\n")
            option_input = input("Enter your options:\n").strip()

            if option_input == "1":
                if isinstance(output_file, list):
                    for file in output_file:
                        sbatch_manager.create_and_run_sbatch_script_add_h(pdb_id, file)
                else:
                    sbatch_manager.create_and_run_sbatch_script_add_h(pdb_id, output_file)

            elif option_input == "2":
                if isinstance(output_file, list):
                    for file in output_file:
                        sbatch_manager.create_and_run_sbatch_script_add_water(pdb_id, file)
                else:
                    sbatch_manager.create_and_run_sbatch_script_add_water(pdb_id, output_file)

            elif option_input == "3":
                if isinstance(output_file, list):
                    for file in output_file:
                        sbatch_manager.create_and_run_sbatch_script_add_h_and_water(pdb_id, file)
                else:
                    sbatch_manager.create_and_run_sbatch_script_add_h_and_water(pdb_id, output_file)

