#!/home/cwj/opt/anaconda3/bin/python
from ase.io import read, write
from ase.build import sort
import os
import shutil
import json
from scripts.utilities import LDAUini, gen_INCAR, Com_analyzer
import time
import argparse

LDAUref = {"V": 3.25, "Co": 3.32, "Cr": 3.7, "Fe": 5.3, "Mn": 3.9, "Mo": 4.38, "Ni": 6.2, "W": 6.2}

def run_vasp(com_dir, log_file, calc_file, raw_dir):

    current_directory = os.getcwd()

    # Read files in raw_structures and return to raw_structures variable
    raw_structures = os.listdir(raw_dir)

    '''Prepare INCAR, KPOINTS, POTCAR'''

    # Each loop calculates one file
    with open(calc_file, 'r') as f:
        struc_dict = json.load(f)

    for i in raw_structures:
        print(i)
        # Skip already calculated structures, continue calculation
        if i in struc_dict.keys():
            if '0' in struc_dict[i].keys():
                print("This one is already calculated")
                continue
            else:
                # Continue calculation
                current_struc = sort(read(com_dir + '/POSCAR'))
        else:
            current_struc = sort(read(raw_dir + '/' + i))

        write(com_dir + '/POSCAR', current_struc)
        print("Currently calculating structure:", str(current_struc.symbols))

        # Skip this computation if the number of atoms > 40
        if len(current_struc.symbols) > 40:
            write('skipped/' + str(current_struc.symbols) + ".vasp", current_struc)
            continue

        # Count the number of Mg in the structure
        Mg_num = [i for i in current_struc.symbols].count("Mg")
        # Skip this computation if number of Mg > 12
        if Mg_num > 12:
            write('skipped/' + str(current_struc.symbols) + ".vasp", current_struc)
            continue

        with open(log_file, "a") as file:
            file.write(str(current_struc.symbols) + "\n")

        ldauflag = LDAUini(current_struc, LDAUref.keys())

        # Generate INCAR based on the current structure
        gen_INCAR(current_struc, ldauflag, com_dir)

        # Generate KPOINTS, POTCAR and run VASP
        os.chdir(com_dir)
        os.system("rm POTCAR")
        os.system('echo -e "102 \n1 \n0.03" |vaspkit > cache')
        time1 = time.time()
        os.system('mpirun -n 48 ~/opt/vasp.5.4.4/bin/vasp_std >> log_run')
        os.chdir(current_directory)

        # Read energy and volume after calculation
        energy, volume = Com_analyzer("./" + com_dir)
        try:
            print("Structure energy: {:.5g} eV".format(float(energy)))
            print("Structure volume: {:.5g} A^3".format(float(volume)))
        except:
            write('crashed_calc/' + str(current_struc.symbols) + ".vasp", current_struc)

        # Store energy and volume in json
        try:
            struc_dict[i][Mg_num] = {"Structures": str(current_struc.symbols), "Energy": energy, "Volume": volume}
        except:
            struc_dict[i] = {Mg_num: {"Structures": str(current_struc.symbols), "Energy": energy, "Volume": volume}}

        with open(calc_file, 'w') as f:
            json.dump(struc_dict, f)

        time10 = time.time()
        with open(log_file, 'a') as f:
            f.write(f"Energy: {float(energy):.2f} Volume: {float(volume):.2f} Time: {time10 - time1:.2f}\n")

        # Save CONTCAR to opt_structures
        dir_name = i + str(current_struc.symbols)

        if not os.path.exists('./opt_structures/' + dir_name):
            os.mkdir('./opt_structures/' + dir_name)
        shutil.copy2('./' + com_dir + '/CONTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + ".vasp")
        shutil.copy2('./' + com_dir + '/OUTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + "_outcar.vasp")

        # Gradually remove all Mg
        for j in range(Mg_num - 1, -1, -1):
            time3 = time.time()

            # Remove one Mg and rerun the structure optimization
            os.system("./scripts/del_atom.py ./" + com_dir + "/CONTCAR Mg 1 ./" + com_dir + "/POSCAR")

            # Note: read POSCAR again and write to a new file name
            current_struc = read("./" + com_dir + "/POSCAR")
            print("Currently calculating structure:", str(current_struc.symbols))

            # Regenerate INCAR to avoid unreasonable LDAU atoms after removing all Mg
            gen_INCAR(current_struc, ldauflag, com_dir)

            os.chdir(com_dir)
            os.system("rm POTCAR")  # Regenerate POTCAR
            os.system('vaspkit -task 103 > cache')
            os.system('mpirun -n 48 ~/opt/vasp.5.4.4/bin/vasp_std >> log_run')
            os.chdir(current_directory)

            # Read energy and volume from OSZICAR after calculation
            energy, volume = Com_analyzer("./" + com_dir)
            try:
                print("Structure energy: {:.5g} eV".format(float(energy)))
                print("Structure volume: {:.5g} A^3".format(float(volume)))
            except:
                write('crashed_calc/' + str(current_struc.symbols) + ".vasp", current_struc)
                continue

            # Store energy and volume in json
            struc_dict[i][j] = {"Structures": str(current_struc.symbols), "Energy": energy, "Volume": volume}
            with open(calc_file, 'w') as f:
                json.dump(struc_dict, f)

            time4 = time.time()
            with open(log_file, 'a') as f:
                f.write(f"Energy: {float(energy):.2f} Volume: {float(volume):.2f} Time: {time4 - time3:.2f}\n")

            # Save CONTCAR to opt_structures
            shutil.copy2('./' + com_dir + '/CONTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + ".vasp")
            shutil.copy2('./' + com_dir + '/OUTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + "_outcar.vasp")

        # Calculate the time spent on this material
        time2 = time.time()
        elapsed_time = time2 - time1
        print(f"Total time elapsed: {elapsed_time:.2f} seconds")
        with open(log_file, 'a') as f:
            f.write(f"Total time elapsed: {elapsed_time:.2f} seconds\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run VASP calculations with specified parameters.")
    parser.add_argument("--com_dir", type=str, default="set1/com", help="Directory for COM files")
    parser.add_argument("--log_file", type=str, default="set1/log", help="Log file name")
    parser.add_argument("--calc_file", type=str, default="set1/calc.json", help="Calculation file name")
    parser.add_argument("--raw_dir", type=str, default="set1/raw", help="Directory for raw structures")

    args = parser.parse_args()
    run_vasp(args.com_dir, args.log_file, args.calc_file, args.raw_dir)
#!/home/cwj/opt/anaconda3/bin/python
from ase.io import read, write
from ase.build import sort
import os
import shutil
import json
from scripts.utilities import LDAUini, gen_INCAR, Com_analyzer
import time
import argparse

LDAUref = {"V": 3.25, "Co": 3.32, "Cr": 3.7, "Fe": 5.3, "Mn": 3.9, "Mo": 4.38, "Ni": 6.2, "W": 6.2}

def run_vasp(com_dir, log_file, calc_file, raw_dir):

    current_directory = os.getcwd()

    # Read files in raw_structures and return to raw_structures variable
    raw_structures = os.listdir(raw_dir)

    '''Prepare INCAR, KPOINTS, POTCAR'''

    # Each loop calculates one file
    with open(calc_file, 'r') as f:
        struc_dict = json.load(f)

    for i in raw_structures:
        print(i)
        # Skip already calculated structures, continue calculation
        if i in struc_dict.keys():
            if '0' in struc_dict[i].keys():
                print("This one is already calculated")
                continue
            else:
                # Continue calculation
                current_struc = sort(read(com_dir + '/POSCAR'))
        else:
            current_struc = sort(read(raw_dir + '/' + i))

        write(com_dir + '/POSCAR', current_struc)
        print("Currently calculating structure:", str(current_struc.symbols))

        # Skip this computation if the number of atoms > 40
        if len(current_struc.symbols) > 40:
            write('skipped/' + str(current_struc.symbols) + ".vasp", current_struc)
            continue

        # Count the number of Mg in the structure
        Mg_num = [i for i in current_struc.symbols].count("Mg")
        # Skip this computation if number of Mg > 12
        if Mg_num > 12:
            write('skipped/' + str(current_struc.symbols) + ".vasp", current_struc)
            continue

        with open(log_file, "a") as file:
            file.write(str(current_struc.symbols) + "\n")

        ldauflag = LDAUini(current_struc, LDAUref.keys())

        # Generate INCAR based on the current structure
        gen_INCAR(current_struc, ldauflag, com_dir)

        # Generate KPOINTS, POTCAR and run VASP
        os.chdir(com_dir)
        os.system("rm POTCAR")
        os.system('echo -e "102 \n1 \n0.03" |vaspkit > cache')
        time1 = time.time()
        os.system('mpirun -n 48 ~/opt/vasp.5.4.4/bin/vasp_std >> log_run')
        os.chdir(current_directory)

        # Read energy and volume after calculation
        energy, volume = Com_analyzer("./" + com_dir)
        try:
            print("Structure energy: {:.5g} eV".format(float(energy)))
            print("Structure volume: {:.5g} A^3".format(float(volume)))
        except:
            write('crashed_calc/' + str(current_struc.symbols) + ".vasp", current_struc)

        # Store energy and volume in json
        try:
            struc_dict[i][Mg_num] = {"Structures": str(current_struc.symbols), "Energy": energy, "Volume": volume}
        except:
            struc_dict[i] = {Mg_num: {"Structures": str(current_struc.symbols), "Energy": energy, "Volume": volume}}

        with open(calc_file, 'w') as f:
            json.dump(struc_dict, f)

        time10 = time.time()
        with open(log_file, 'a') as f:
            f.write(f"Energy: {float(energy):.2f} Volume: {float(volume):.2f} Time: {time10 - time1:.2f}\n")

        # Save CONTCAR to opt_structures
        dir_name = i + str(current_struc.symbols)

        if not os.path.exists('./opt_structures/' + dir_name):
            os.mkdir('./opt_structures/' + dir_name)
        shutil.copy2('./' + com_dir + '/CONTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + ".vasp")
        shutil.copy2('./' + com_dir + '/OUTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + "_outcar.vasp")

        # Gradually remove all Mg
        for j in range(Mg_num - 1, -1, -1):
            time3 = time.time()

            # Remove one Mg and rerun the structure optimization
            os.system("./scripts/del_atom.py ./" + com_dir + "/CONTCAR Mg 1 ./" + com_dir + "/POSCAR")

            # Note: read POSCAR again and write to a new file name
            current_struc = read("./" + com_dir + "/POSCAR")
            print("Currently calculating structure:", str(current_struc.symbols))

            # Regenerate INCAR to avoid unreasonable LDAU atoms after removing all Mg
            gen_INCAR(current_struc, ldauflag, com_dir)

            os.chdir(com_dir)
            os.system("rm POTCAR")  # Regenerate POTCAR
            os.system('vaspkit -task 103 > cache')
            os.system('mpirun -n 48 ~/opt/vasp.5.4.4/bin/vasp_std >> log_run')
            os.chdir(current_directory)

            # Read energy and volume from OSZICAR after calculation
            energy, volume = Com_analyzer("./" + com_dir)
            try:
                print("Structure energy: {:.5g} eV".format(float(energy)))
                print("Structure volume: {:.5g} A^3".format(float(volume)))
            except:
                write('crashed_calc/' + str(current_struc.symbols) + ".vasp", current_struc)
                continue

            # Store energy and volume in json
            struc_dict[i][j] = {"Structures": str(current_struc.symbols), "Energy": energy, "Volume": volume}
            with open(calc_file, 'w') as f:
                json.dump(struc_dict, f)

            time4 = time.time()
            with open(log_file, 'a') as f:
                f.write(f"Energy: {float(energy):.2f} Volume: {float(volume):.2f} Time: {time4 - time3:.2f}\n")

            # Save CONTCAR to opt_structures
            shutil.copy2('./' + com_dir + '/CONTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + ".vasp")
            shutil.copy2('./' + com_dir + '/OUTCAR', './opt_structures/' + dir_name + "/" + str(current_struc.symbols) + "_outcar.vasp")

        # Calculate the time spent on this material
        time2 = time.time()
        elapsed_time = time2 - time1
        print(f"Total time elapsed: {elapsed_time:.2f} seconds")
        with open(log_file, 'a') as f:
            f.write(f"Total time elapsed: {elapsed_time:.2f} seconds\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run VASP calculations with specified parameters.")
    parser.add_argument("--com_dir", type=str, default="set1/com", help="Directory for COM files")
    parser.add_argument("--log_file", type=str, default="set1/log", help="Log file name")
    parser.add_argument("--calc_file", type=str, default="set1/calc.json", help="Calculation file name")
    parser.add_argument("--raw_dir", type=str, default="set1/raw", help="Directory for raw structures")

    args = parser.parse_args()
    run_vasp(args.com_dir, args.log_file, args.calc_file, args.raw_dir)
