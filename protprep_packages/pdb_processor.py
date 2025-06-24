import subprocess
import os
import requests
from pathlib import Path
from datetime import datetime
import logging
from .ligand_parser import RCSBLigandParser

class PDBProcessor:
    def __init__(self, export_dir='protprep_logging'):
        self.pdb_id = None
        self.ligand_name = None
        self.chain_id = None
        self.cofactor_names = []
        self.cosubstrate_names = []
        self.metal_ion_names = []

        # Setup export directory
        self.export_dir = Path(export_dir)
        self._setup_export_dir()
        
        # Setup logging
        self._setup_logging()
        
        # Log initialization
        #self.logger.info(f"Initialized {self.__class__.__name__} with export directory: {self.export_dir}")

    def _setup_export_dir(self):
        """Create export directory if it doesn't exist"""
        try:
            self.export_dir.mkdir(parents=True, exist_ok=True)
            return True
        except Exception as e:
            print(f"Error creating export directory: {e}")
            return False

    def _setup_logging(self):
        """Configure logging with both file and console handlers"""
        self.logger = logging.getLogger(f"{self.__class__.__name__}")
        self.logger.setLevel(logging.INFO)
        
        # Clear existing handlers
        self.logger.handlers = []
        
        # File handler
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        fh = logging.FileHandler(self.export_dir / f'protprep_{timestamp}.log')
        fh.setLevel(logging.DEBUG)
        
        # Console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        
        # Formatter
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

    def wget_pdb_file(self):
        """
        Download PDB file using wget or requests as fallback.
        Returns the path to the downloaded file.
        """
        if not self.pdb_id:
            raise ValueError("PDB ID not set. Call get_basic_info() first.")

        url = f"https://files.rcsb.org/download/{self.pdb_id}.pdb"
        output_file = f"{self.pdb_id}.pdb"

        try:
            # First attempt: using wget
            result = subprocess.run(
                ["wget", url, "-O", output_file],
                capture_output=True,
                text=True,
                check=True
            )
            self.pdb_file_path = output_file
            self.logger.info(f"Successfully downloaded {output_file} using wget")
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Fallback: using requests if wget fails or is not available
            try:
                self.logger.info(f"Wget failed or not found, trying direct download...")
                response = requests.get(url)
                response.raise_for_status()
                
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                
                self.pdb_file_path = output_file
                self.logger.info(f"Successfully downloaded {output_file} using direct download")
                
            except requests.exceptions.RequestException as e:
                self.logger.info(f"Error downloading file: {e}")
                return None

        # Verify file exists and has content
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            return Path(output_file).absolute()
        else:
            self.logger.info(f"Download failed or file is empty")
            return None

    def get_basic_info(self):
        """Get one or more PDB and ligand information from user."""
        self.pdb_data = {}

        mode = input("Do you want to enter a single PDB file or multiple? (Enter 'single' or 'multiple'):\n").strip().lower()

        if mode == 'single':
            pdb_id = input("\nEnter the PDB ID in uppercase (e.g., 3ERK):\n")
            self.pdb_id = pdb_id
            pdb_file = self.wget_pdb_file()
            analyzer = RCSBLigandParser()
            analyzer.fetch_structure_page(pdb_id)
            analyzer.display_results()
            ligand = input("\nEnter the ligand of interest in uppercase (e.g., SB4):\n").upper()
            self.pdb_data[pdb_id] = {'pdb_file': pdb_file, 'ligand': ligand}
        elif mode == 'multiple':
            while True:
                pdb_id = input("\nEnter the PDB ID in uppercase (e.g., 3ERK) or type 'done' to finish:\n")
                if pdb_id.lower() == 'done':
                    break
                self.pdb_id = pdb_id
                pdb_file = self.wget_pdb_file()
                analyzer = RCSBLigandParser()
                analyzer.fetch_structure_page(pdb_id)
                analyzer.display_results()
                ligand = input("\nEnter the ligand of interest in uppercase (e.g., SB4):\n")
                self.pdb_data[pdb_id] = {'pdb_file': pdb_file, 'ligand': ligand}
        else:
            print("Invalid input. Please enter 'single' or 'multiple'.")
            return None

        return self.pdb_data 

    def process_chain_extraction(self):
        """Process chain extraction preferences."""
        extract_one_chain = input("\nDo you want to extract only one chain? (Y/N):\n")
        if extract_one_chain.upper() == "Y":
            self.chain_id = input("\nEnter the chain ID to extract (e.g., A):\n")

    def process_cofactors(self):
        """Process cofactor preferences."""
        extract_cofactor = input("\nDo you want to keep the cofactors? (Y/N):\n")
        if extract_cofactor.upper() == "Y":
            ask_more_than_one = input("Is the number of cofactors that you want to keep more than one? (Y/N):\n")
            
            if ask_more_than_one.upper() != "Y":
                cofactor = input("Enter the name of the cofactor you want to extract (e.g., HEM):\n")
                self.cofactor_names = [cofactor]
            else:
                number_of_cofactors = int(input("How many cofactors do you want to keep? (e.g., 3):\n"))
                for _ in range(number_of_cofactors):
                    cofactor = input("Enter the name of a cofactor:\n")
                    self.cofactor_names.append(cofactor)
                print("Cofactors to be kept:\n", self.cofactor_names)
                self.logger.info("Cofactors to be kept: ", self.cofactor_names)

    def process_cosubstrates(self):
        """Process cosubstrate preferences."""
        extract_cosubstrate = input("\nDo you want to keep the cosubstrate? (Y/N):\n")
        if extract_cosubstrate.upper() == "Y":
            ask_more_than_one = input("Is the number of cosubstrate that you want to keep more than one? (Y/N):\n")
            
            if ask_more_than_one.upper() != "Y":
                cosubstrate = input("Enter the name of the cosubstrate you want to extract (e.g., HEM):\n")
                self.cosubstrate_names = [cosubstrate]
            else:
                number_of_cosubstrates = int(input("How many cosubstrate do you want to keep? (e.g., 3):\n"))
                for _ in range(number_of_cosubstrates):
                    cosubstrate = input("Enter the name of a cosubstrate:\n")
                    self.cosubstrate_names.append(cosubstrate)
                print("Cosubstrate to be kept:\n", self.cosubstrate_names)
                self.logger.info("Cosubstrate to be kept: ", self.cosubstrate_names)

    def process_metal_ions(self):
        """Process metal ion preferences."""
        extract_metal_ion = input("\nDo you want to keep the metal ions? (Y/N):\n")
        if extract_metal_ion.upper() == "Y":
            ask_more_than_one = input("Is the number of metal ions that you want to keep more than one? (Y/N):\n")
            
            if ask_more_than_one.upper() != "Y":
                metal_ion = input("Enter the name of the metal ion that you want to keep (e.g., ZN):\n")
                self.metal_ion_names = [metal_ion]
            else:
                number_of_metal_ions = int(input("How many metal ions do you want to keep? (e.g., 3):\n"))
                for _ in range(number_of_metal_ions):
                    metal_ion = input("Enter the name of a metal ion:\n")
                    self.metal_ion_names.append(metal_ion)
                print("Metal ions to be kept:\n", self.metal_ion_names)
                self.logger.info("Metal ions to be kept: ", self.metal_ion_names)

    def process_all_entries(self):
        """Process all user inputs for each PDB and store results per PDB."""
        # Ensure processed_results is empty before starting
        self.processed_results = {}
        
        for pdb_id in self.pdb_data.keys():
            self.pdb_id = pdb_id

            # Create the message and fancy box borders
            message = f"Processing PDB entry: {pdb_id}"
            width = len(message) + 2
            top_border = "┏" + "━" * width + "┓"
            bottom_border = "┗" + "━" * width + "┛"

            # Print the boxed message
            print("\n" + top_border)
            print("┃ " + message + " ┃")
            print(bottom_border)
                        
            # Process chain extraction for this pdb
            self.process_chain_extraction()
            # Process cofactors for this pdb
            self.process_cofactors()
            # Process cosubstrates for this pdb
            self.process_cosubstrates()
            # Process metal ions for this pdb
            self.process_metal_ions()
            
            # Store the processed results for this pdb
            self.processed_results[pdb_id] = {
                'pdb_id': pdb_id,
                'ligand_name': self.pdb_data[pdb_id]['ligand'],
                'chain_id': self.chain_id,
                'cofactors': self.cofactor_names,
                'cosubstrates': self.cosubstrate_names,
                'metal_ions': self.metal_ion_names
            }
            
            # Reset instance attributes for the next pdb
            self.chain_id = None
            self.cofactor_names = []
            self.cosubstrate_names = []
            self.metal_ion_names = []

        self.logger.info("Processed results: %s", self.processed_results)


