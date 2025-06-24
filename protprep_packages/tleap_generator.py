class TleapInputGenerator:
    def __init__(self):
        """Initialize TleapInputGenerator."""
        self.base_sources = [
            "source leaprc.protein.ff14SB",
            "source leaprc.water.tip3p",
            "source leaprc.DNA.bsc1",
            "source leaprc.gaff"
        ]

    def generate_multiple_metals_input(self, pdb_id):
        """Generate tleap input for multiple metal ions configuration."""
        tleap_input = f"""\
{self._get_sources()}

# Load the protein and the ligand (cofactor)
protein = loadpdb rec_{pdb_id}.pdb
metal_ion = loadpdb {pdb_id}_all_metal_ions_single_chain.pdb

# Combine the protein and cofactor into a single structure
combined_structure = combine {{ protein metal_ion }}

# Check and save the combined structure
check combined_structure
savepdb combined_structure rec_{pdb_id}_clean.pdb
"""
        tleap_input_filename = f"{pdb_id}_tleap_2S_multiple_metals.in"
        self._write_tleap_file(tleap_input_filename, tleap_input)

    def generate_three_component_input(self, pdb_id, cofactor_name, metal_ion_name):
        """Generate tleap input for protein, cofactor, and metal ions."""
        tleap_input = f"""\
{self._get_sources()}

# Load the protein, cofactor, and metal_ion files
protein = loadpdb rec_{pdb_id}.pdb
cofactor = loadpdb {pdb_id}_single_chain_with_{cofactor_name}.pdb
metal_ion = loadpdb {pdb_id}_single_chain_with_{metal_ion_name}.pdb

# Combine the protein and cofactor into a single structure
combined_structure = combine {{ protein cofactor metal_ion}}

# Check and save the combined structure
check combined_structure
savepdb combined_structure rec_{pdb_id}_clean.pdb
"""
        tleap_input_filename = f"{pdb_id}_tleap_3S.in"
        self._write_tleap_file(tleap_input_filename, tleap_input)

    def generate_protein_only_input(self, pdb_id):
        """Generate tleap input for just protein."""
        tleap_input = f"""\
{self._get_sources()}

protein = loadpdb rec_{pdb_id}.pdb
check protein
savepdb protein rec_{pdb_id}_clean.pdb
"""
        tleap_input_filename = f"{pdb_id}_tleap.in"
        self._write_tleap_file(tleap_input_filename, tleap_input)

    def generate_protein_cofactor_input(self, pdb_id, cofactor_name):
        """Generate tleap input for protein and cofactor."""
        tleap_input = f"""\
{self._get_sources()}

# Load the protein and the ligand (cofactor)
protein = loadpdb rec_{pdb_id}.pdb
cofactor = loadpdb {pdb_id}_single_chain_with_{cofactor_name}.pdb

# Combine the protein and cofactor into a single structure
combined_structure = combine {{ protein cofactor }}

# Check and save the combined structure
check combined_structure
savepdb combined_structure rec_{pdb_id}_clean.pdb
"""
        tleap_input_filename = f"{pdb_id}_tleap_2S_PC.in"
        self._write_tleap_file(tleap_input_filename, tleap_input)

    def generate_protein_metal_input(self, pdb_id, metal_ion_name):
        """Generate tleap input for protein and metal ion."""
        tleap_input = f"""\
{self._get_sources()}

# Load the protein, cofactor, and metal_ion files
protein = loadpdb rec_{pdb_id}.pdb
metal_ion = loadpdb {pdb_id}_single_chain_with_{metal_ion_name}.pdb

# Combine the protein and cofactor into a single structure
combined_structure = combine {{ protein metal_ion}}

# Check and save the combined structure
check combined_structure
savepdb combined_structure rec_{pdb_id}_clean.pdb
"""
        tleap_input_filename = f"{pdb_id}_tleap_2S_PMI.in"
        self._write_tleap_file(tleap_input_filename, tleap_input)

    def generate_all_inputs(self, pdb_id, cofactor_name=None, metal_ion_name=None):
        """Generate all possible tleap input combinations based on available components."""
        self.generate_multiple_metals_input(pdb_id)
        if cofactor_name and metal_ion_name:
            self.generate_three_component_input(pdb_id, cofactor_name, metal_ion_name)
        self.generate_protein_only_input(pdb_id)
        if cofactor_name:
            self.generate_protein_cofactor_input(pdb_id, cofactor_name)
        if metal_ion_name:
            self.generate_protein_metal_input(pdb_id, metal_ion_name)

    def _get_sources(self):
        """Return the standard source commands as a string."""
        return "\n".join(self.base_sources)

    def _write_tleap_file(self, filename, content):
        """Write the tleap input content to a file."""
        with open(filename, "w") as tleap_file:
            tleap_file.write(content)