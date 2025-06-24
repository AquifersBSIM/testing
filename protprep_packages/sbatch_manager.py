import subprocess
import argparse 

def create_and_run_sbatch_script_add_h(pdb_id=None, pdb_file=None):

    sbatch_script = f"""\
#!/bin/bash
#
"""
    sbatch_script += f"#SBATCH --job-name=rec_protein_{pdb_id}_add_h\n"
    sbatch_script += f"#SBATCH --output=rec_protein_{pdb_id}_add_h.txt\n"
    sbatch_script += "#\n"
    sbatch_script += "#SBATCH --ntasks=1\n"
    sbatch_script += "#SBATCH --cpus-per-task=1\n"
    sbatch_script += "#SBATCH --ntasks-per-node=1\n"
    sbatch_script += "#SBATCH --time=00:10:00\n"
    sbatch_script += "#SBATCH --mem-per-cpu=2G\n"

    sbatch_script += f"ml conda\n"
    sbatch_script += f"conda activate /fred/oz241/BSIM/conda_meeko\n"
    sbatch_script += f"/usr/bin/time -v python protprep_packages/add_h.py --pdb_file {pdb_file}\n"

    # Write the sbatch script to a file
    script_filename = f"rec_protein_{pdb_id}_add_h.sh"
    with open(script_filename, "w") as sbatch_file:
        sbatch_file.write(sbatch_script)

    # Launch the sbatch script using subprocess
    subprocess.run(["sbatch", script_filename])

def create_and_run_sbatch_script_add_water(pdb_id=None, pdb_file=None):

    # Create the Add Water SBATCH script with dependency on Add Hydrogen job
    add_water_script = f"""\
#!/bin/bash
#
#SBATCH --job-name=rec_protein_{pdb_id}_add_water
#SBATCH --output=rec_protein_{pdb_id}_add_water.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2G

ml conda
conda activate /fred/oz241/BSIM/conda_meeko
/usr/bin/time -v python protprep_packages/add_water.py --pdb_file {pdb_file}
"""

    water_script_filename = f"rec_protein_{pdb_id}_add_water.sh"
    with open(water_script_filename, "w") as sbatch_file:
        sbatch_file.write(add_water_script)

    # Launch the sbatch script using subprocess
    subprocess.run(["sbatch", water_script_filename])

def create_and_run_sbatch_script_add_h_and_water(pdb_id=None, pdb_file=None):
    # Create the Add Hydrogen SBATCH script
    add_h_script = f"""\
#!/bin/bash
#
#SBATCH --job-name=rec_protein_{pdb_id}_add_h
#SBATCH --output=rec_protein_{pdb_id}_add_h.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2G

ml conda
conda activate /fred/oz241/BSIM/conda_meeko
/usr/bin/time -v python protprep_packages/add_h.py --pdb_file {pdb_file}
"""

    h_script_filename = f"rec_protein_{pdb_id}_add_h.sh"
    with open(h_script_filename, "w") as sbatch_file:
        sbatch_file.write(add_h_script)

    # Submit the Add Hydrogen job and capture the Job ID
    result = subprocess.run(["sbatch", h_script_filename], capture_output=True, text=True)
    job_id = None
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]  # Extract the Job ID from sbatch output
        print(f"Submitted add_h job with ID {job_id}")
    else:
        print(f"Failed to submit add_h job: {result.stderr}")
        return

    # Create the Add Water SBATCH script with dependency on Add Hydrogen job
    add_water_script = f"""\
#!/bin/bash
#
#SBATCH --job-name=rec_protein_{pdb_id}_add_water
#SBATCH --output=rec_protein_{pdb_id}_add_water.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2G
#SBATCH --dependency=afterok:{job_id}  # Dependency on add_h job

ml conda
conda activate /fred/oz241/BSIM/conda_meeko
/usr/bin/time -v python protprep_packages/add_water.py --pdb_file {pdb_file}
"""

    water_script_filename = f"rec_protein_{pdb_id}_add_water.sh"
    with open(water_script_filename, "w") as sbatch_file:
        sbatch_file.write(add_water_script)

    # Submit the Add Water job
    result = subprocess.run(["sbatch", water_script_filename], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Submitted add_water job with dependency on {job_id}")
    else:
        print(f"Failed to submit add_water job: {result.stderr}")

def create_and_run_sbatch_script_extract_and_clean_specific_ligands(pdb_id, pdb_file, ligand_to_keep, chain_id):

    sbatch_script = f"""\
#!/bin/bash
#
"""
    sbatch_script += f"#SBATCH --job-name=extracting_lig_{pdb_id}\n"
    sbatch_script += f"#SBATCH --output=extracting_lig_{pdb_id}.txt\n"
    sbatch_script += "#\n"
    sbatch_script += "#SBATCH --ntasks=1\n"
    sbatch_script += "#SBATCH --cpus-per-task=1\n"
    sbatch_script += "#SBATCH --ntasks-per-node=1\n"
    sbatch_script += "#SBATCH --time=00:10:00\n"
    sbatch_script += "#SBATCH --mem-per-cpu=2G\n"

    sbatch_script += f"ml conda\n"
    sbatch_script += f"conda activate /fred/oz241/BSIM/conda_meeko\n"
    sbatch_script += f"/usr/bin/time -v python protprep_packages/extract_and_clean_specific_ligands.py --pdb_file {pdb_file} --pdb_id {pdb_id} --ligand_to_keep {ligand_to_keep} --chain_id {chain_id}\n"

    # Write the sbatch script to a file
    script_filename = f"extracting_lig_{pdb_id}.sh"
    with open(script_filename, "w") as sbatch_file:
        sbatch_file.write(sbatch_script)

    # Launch the sbatch script using subprocess
    subprocess.run(["sbatch", script_filename])
