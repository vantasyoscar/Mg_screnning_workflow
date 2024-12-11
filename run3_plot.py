import os
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.msd import EinsteinMSD

def calculate_msd(u, atom_type='Mg'):
    # Select the atom type

    # Calculate MSD
    msd = EinsteinMSD(u, select=f'type {atom_type}', fft=False)

    # Run analysis
    msd.run()

    return msd.results.timeseries

def find_mg_position(filename):
    with open(filename, 'r') as file:
        for line in file:
            words = line.split()
            if len(words) > 0 and words[0] == 'pair_coeff':
                if "Mg" in words:
                    return words.index("Mg") + 1
    return None

# Traverse all subfolders in the com1 directory
msd_dict = {}
plt.figure(dpi=160)

for root, dirs, files in os.walk('MLMD'):
    # If there is a "finished" file in the subfolder
    if 'finished' in files:
        # Read dump.lammpstrj trajectory
        position = find_mg_position(os.path.join(root, 'in.lammps')) - 4
        u = mda.Universe(os.path.join(root, 'dump.lammpstrj'), format="LAMMPSDUMP")
        
        # Calculate MSD of Mg
        msd = calculate_msd(u, atom_type=str(position))
        msd_dict[root] = msd
        
        # Plot the results
        plt.plot(np.arange(5001)/100, msd)  # timestep is 10fs
        plt.xlabel('Time (fs)')
        plt.ylabel('MSD of Mg')
        plt.title(f'MSD of all Mg')

import pandas as pd
df = pd.DataFrame(msd_dict)
df.to_csv("msd.csv")

msd = {}
for i, j in msd_dict.items():
    msd[i[4:-4]] = j[-1]

import csv

# The keys of the dictionary as column names of the CSV file
fieldnames = msd.keys()

# Create CSV file
with open('msd50.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write the column names
    writer.writeheader()

    # Write the data
    writer.writerow(msd)
