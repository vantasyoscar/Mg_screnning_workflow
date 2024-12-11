from ase import Atoms
from ase.io import read, write
from ase.build import sort
import os
import shutil
import json
import time
from ase.data import atomic_masses, chemical_symbols
import random
import numpy as np
import itertools

def scale_structure(structure, min_atoms=256, max_atoms=1024):
    # Get original lattice parameters and number of atoms
    cell = structure.get_cell()
    num_atoms = structure.get_number_of_atoms()

    # Calculate the minimum expansion factor required
    min_factor = int(np.cbrt(min_atoms / num_atoms))

    # Find an expansion factor that brings the number of atoms within the specified range and makes a, b, c as close as possible
    best_structure = None
    min_difference = float('inf')
    for factor in itertools.product(range(min_factor, int(np.cbrt(max_atoms / num_atoms)) + 1), repeat=3):
        new_structure = structure * factor
        if min_atoms <= new_structure.get_number_of_atoms() <= max_atoms:
            new_cell = new_structure.get_cell()
            difference = max(new_cell.diagonal()) - min(new_cell.diagonal())
            if difference < min_difference:
                min_difference = difference
                best_structure = new_structure

    if best_structure is not None:
        return sort(best_structure)
    else:
        print("No suitable expansion factor found")
        return structure

# Generate LAMMPS script for calculating this structure
def generate_lmp(atoms, output_path, model_path="../model.pth"):
    
    current_symbols = ""
    for i in atoms.get_chemical_symbols():
        if not i in current_symbols:
            current_symbols += i + " "
    current_symbols = current_symbols[:-1]

    cmds = ["units           metal",
            "atom_style      atomic",
            "newton          off",
            "read_data       lammps.data", 
            "pair_style      nequip",
            "pair_coeff      * * " + model_path + " " + current_symbols,
            ]
    
    # Add the mass for each atom type
    for j, atom_type in enumerate(current_symbols.split(" ")):
        mass = atomic_masses[chemical_symbols.index(atom_type)]
        cmds.append(f"mass {str(j+1)} {mass}")

    cmds.extend(["thermo_style    custom step pe ke etotal temp press vol",
                "thermo          100",
                "timestep        0.001",
                "velocity all create 1200.0 223323333",
                "fix             1 all npt temp 1200 1200 0.1 iso 1.0 1.0 0.5",
                "dump            1 all custom 10 warmup.lammpstrj id type element xu yu zu",
                f"dump_modify     1 element " + current_symbols,
                "run             10000",
                "unfix 1",
                "undump 1",
                "fix             2 all nvt temp 1200 1200 0.1",
                "dump            2 all custom 10 dump.lammpstrj id type element xu yu zu",
                f"dump_modify     2 element " + current_symbols,
                "run             50000",
                "write_data finished"]
                )

    if not os.path.exists(output_path):
        # Create directory if it does not exist
        os.makedirs(output_path)
    # Write LAMMPS input script to file
    with open(output_path + '/in.lammps', 'w') as f:
        for cmd in cmds:
            f.write(cmd + '\n')

# Each loop is a calculation for one file
# Read files from raw_structures and store them in raw_struc variable
raw_structures = os.listdir("high_V")

for i in raw_structures:
    output_path = "./MLMD/" + i[:-4] + "calc"
    # Place structure data in the calculation folder

    current_struc = sort(read('high_V/' + i))
    current_struc = scale_structure(current_struc)
    mg_indices = [atom.index for atom in current_struc if atom.symbol == 'Mg']

    # Remove 1/5 of Mg
    for j in range(len(mg_indices) // 5):
        mg_indices = [atom.index for atom in current_struc if atom.symbol == 'Mg']

        if len(mg_indices) > 1:
            # Randomly select an Mg atom
            mg_to_remove = random.choice(mg_indices)
            # Remove the selected Mg atom
            del current_struc[mg_to_remove]

    print("Currently processing structure: ", str(current_struc.symbols))

    generate_lmp(current_struc, output_path, "../../models/model_Nequip.pth")

    write(output_path + '/lammps.data', current_struc, format='lammps-data')
