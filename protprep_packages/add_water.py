#!/usr/bin/env python3
import sys
import numpy as np
import random
import math
import argparse

def parse_pdb(filename):
    """
    Parse the PDB file to extract all atom lines and their coordinates.
    Returns a list of lines (for output later) and a numpy array of atom coordinates.
    """
    pdb_lines = []
    coords = []
    with open(filename, 'r') as f:
        for line in f:
            pdb_lines.append(line.rstrip())
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue
    return pdb_lines, np.array(coords)

def random_rotation_matrix():
    """
    Create a random 3D rotation matrix using Euler angles.
    """
    # random Euler angles (in radians)
    a, b, c = [random.uniform(0, 2*math.pi) for _ in range(3)]
    Rx = np.array([[1, 0, 0],
                   [0, math.cos(a), -math.sin(a)],
                   [0, math.sin(a), math.cos(a)]])
    Ry = np.array([[math.cos(b), 0, math.sin(b)],
                   [0, 1, 0],
                   [-math.sin(b), 0, math.cos(b)]])
    Rz = np.array([[math.cos(c), -math.sin(c), 0],
                   [math.sin(c), math.cos(c), 0],
                   [0, 0, 1]])
    return Rz @ Ry @ Rx

def create_water_molecule(origin):
    """
    Create a water molecule (O, H1, H2) with the oxygen at 'origin'.
    The water geometry is based on:
      O at (0,0,0)
      H1 at (0.9572, 0, 0)
      H2 at (-0.2399872, 0.927297, 0)
    A random rotation is applied to the molecule.
    Returns a list of tuples: (atom_name, x, y, z)
    """
    # Standard water geometry in Å
    water_template = {
        'O': np.array([0.0, 0.0, 0.0]),
        'H1': np.array([0.9572, 0.0, 0.0]),
        'H2': np.array([-0.2399872, 0.927297, 0.0])
    }
    R = random_rotation_matrix()
    water_atoms = []
    for atom, pos in water_template.items():
        pos_rot = R @ pos + origin
        water_atoms.append((atom, pos_rot[0], pos_rot[1], pos_rot[2]))
    return water_atoms

def write_pdb(filename, pdb_lines, water_mols, start_atom_number, start_residue_number):
    """
    Write a new PDB file that includes the original lines and the added waters.
    water_mols is a list where each element is a list of water atoms.
    """
    with open(filename, 'w') as f:
        # Write original PDB lines
        for line in pdb_lines:
            f.write(line + "\n")
        # Write water molecules
        atom_number = start_atom_number
        res_number = start_residue_number
        for water in water_mols:
            for atom in water:
                # Format: ATOM, atom number, atom name, residue name, chain, residue number, x,y,z, occupancy, B-factor, element
                record = "ATOM  {atom_number:5d} {atom_name:^4s} HOH A{res_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {element:>2s}"
                element = atom[0][0]  # first letter (e.g., O or H)
                f.write(record.format(atom_number=atom_number,
                                        atom_name=atom[0],
                                        res_number=res_number,
                                        x=atom[1],
                                        y=atom[2],
                                        z=atom[3],
                                        element=element) + "\n")
                atom_number += 1
            res_number += 1
        f.write("END\n")

def add_waters_to_structure(pdb_coords, margin=5.0, spacing=2.75, cutoff=2.2):
    """
    Create a grid of water molecules around the structure.
    - margin: extra distance (Å) added to the bounding box.
    - spacing: grid spacing in Å.
    - cutoff: distance threshold (Å); water oxygen is not placed if too close to any solute atom.
    Returns a list of water molecules.
    """
    if pdb_coords.size == 0:
        raise ValueError("No coordinates found in the PDB file.")

    min_coords = np.min(pdb_coords, axis=0) - margin
    max_coords = np.max(pdb_coords, axis=0) + margin

    water_mols = []
    # Generate grid points
    x_range = np.arange(min_coords[0], max_coords[0], spacing)
    y_range = np.arange(min_coords[1], max_coords[1], spacing)
    z_range = np.arange(min_coords[2], max_coords[2], spacing)

    for x in x_range:
        for y in y_range:
            for z in z_range:
                point = np.array([x, y, z])
                # Check if point is at least cutoff away from any solute atom.
                dists = np.linalg.norm(pdb_coords - point, axis=1)
                if np.min(dists) > cutoff:
                    water_mols.append(create_water_molecule(point))
    return water_mols

def main(input_pdb):

    output_pdb = "test.pdb"

    # Parse original PDB
    pdb_lines, pdb_coords = parse_pdb(input_pdb)
    print(f"Read {len(pdb_coords)} atoms from {input_pdb}")

    # Create water molecules
    water_mols = add_waters_to_structure(pdb_coords)
    print(f"Added {len(water_mols)} water molecules.")

    # Determine starting atom and residue numbers for water molecules
    last_atom_number = 0
    last_residue_number = 0
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom_num = int(line[6:11])
                res_num = int(line[22:26])
                last_atom_number = max(last_atom_number, atom_num)
                last_residue_number = max(last_residue_number, res_num)
            except ValueError:
                continue

    start_atom_number = last_atom_number + 1
    start_residue_number = last_residue_number + 1

    # Write output PDB file with added water molecules
    write_pdb(output_pdb, pdb_lines, water_mols, start_atom_number, start_residue_number)
    print(f"Output written to {output_pdb}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create and run an sbatch script using a given pdb_id and pdb_file."
    )

    parser.add_argument(
        "--pdb_file",
        required=True,
        help="Path to the PDB file."
    )
    
    args = parser.parse_args()

    main(args.pdb_file)
