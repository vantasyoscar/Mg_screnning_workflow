import os
import shutil
import csv
import argparse
from ase.io import read

# Argument parser for command-line options
parser = argparse.ArgumentParser(description="Process voltage and capacity thresholds.")
parser.add_argument('--voltage', type=float, default=3.0, help="Voltage threshold")
parser.add_argument('--capacity', type=float, default=800.0, help="Volumetric capacity threshold")
args = parser.parse_args()

# 1. Run the prediction script
os.system("python CGCNN_for_voltage/predict.py models/model_cgcnn.pth.tar raw_structures/")

# 2. Read test_results.csv and add the first column of rows with voltage > args.voltage to high_V
high_V = []
voltage = []
with open('test_results.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if float(row[2]) > args.voltage:
            high_V.append(row[0])
            voltage.append(float(row[2]))

# 3. Copy corresponding files from raw_structures to high_V if volumetric capacity > args.capacity
if not os.path.exists('high_V'):
    os.makedirs('high_V')

high_V_C = []
voltage_high_C = []
capacity = []
for structure, volt in zip(high_V, voltage):
    src_path = os.path.join('raw_structures', f"{structure}.cif")
    dst_path = os.path.join('high_V', f"{structure}.cif")

    if os.path.exists(src_path):
        atoms = read(src_path)
        volume = atoms.get_volume() * 1e-27  # volume in m^3
        e = 1.602176e-19  # elementary charge in Coulombs
        Mg_num = len([atom for atom in atoms if atom.symbol == 'Mg'])
        charge = 2 * Mg_num * e / 3.6
        volumetric_capacity = charge / (volume * 1000)  # capacity in mAh/L

        if volumetric_capacity > args.capacity:
            shutil.copy2(src_path, dst_path)
            high_V_C.append(structure)
            voltage_high_C.append(volt)
            capacity.append(volumetric_capacity)
    else:
        print(f"Warning: {src_path} does not exist")

# 4. Print information
print("Structures with high voltage and volumetric capacity:")
for structure, volt, cap in zip(high_V_C, voltage_high_C, capacity):
    print(f"{structure} Voltage: {volt:.2f} V Capacity: {cap:.0f} Ah/L")
