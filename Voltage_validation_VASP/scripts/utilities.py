import ase
import shutil
import numpy as np

global LDAUref

LDAUref = {"V": 3.25, "Co": 3.32, "Cr": 3.7, "Fe": 5.3, "Mn": 3.9, "Mo": 4.38, "Ni": 6.2, "W": 6.2}

def LDAUini(atoms: ase.atoms.Atoms, ldaulist: list) -> bool:
    r"""
    Determine whether to enable LDA+U
    Args:
        atoms: The structure currently being calculated, ase.symbols.Symbols
        ldaulist: List of element symbols that require LDA+U calculation

    Returns:
        ldauflag: Whether to perform LDA+U calculation
    """
    return any([ele in atoms.symbols for ele in ldaulist])

def gen_INCAR(atoms: ase.atoms.Atoms, ldauflag: int, output="com"):
    r"""
    Generate INCAR file for the current structure
    Args:
        atoms: The structure currently being calculated, ase.symbols.Atoms
        ldauflag: Whether to enable LDA+U calculation
    """

    if ldauflag == 1:
        with open('INCAR', 'r') as f:
            lines = f.readlines()

        with open(output + '/INCAR', 'w+') as f:
            for line in lines:
                if 'LDAU' in line:
                    line = line.replace('False', 'True')  # Enable LDA+U
                f.write(line)

            # Set U values based on recommendation
            LDAUU = "LDAUU = "
            LDAUL = "LDAUL = "
            # Set atom_list according to the order in POSCAR
            lst = [i[1] for i in enumerate(atoms.symbols)]
            atom_list = [lst[i] for i in range(len(lst)) if i == 0 or lst[i] != lst[i-1]]
            for i in atom_list:
                if i in LDAUref.keys():
                    LDAUU += str(LDAUref[i]) + " "
                    LDAUL += "2 "
                else:
                    LDAUU += "0 "
                    LDAUL += "0 "
            f.write(LDAUU + "\n")
            f.write(LDAUL + "\n")
    else:
        # If LDA+U is not needed, directly copy the INCAR file
        shutil.copy2('INCAR', './' + output + '/INCAR')

def Com_analyzer(com_folder: str) -> (float, float):
    r"""
    Read energy and volume from OSZICAR and CONTCAR
    Args:
        com_folder: Path to the VASP calculation folder
    Returns:
        energy: Energy of the current structure
        volume: Volume of the current structure
    """
    # Read information from a specific line and return a list
    def wash(line):
        line = line.split(" ")
        while "" in line:
            line.remove("")
        # Remove the last "\n"
        line[-1] = line[-1][:-1]
        while "" in line:
            line.remove("")
        return line

    # Read the last line of energy from the OSZICAR file
    with open(com_folder + "/OSZICAR", 'r') as f:
        energy = f.readlines()[-1].split()[4]

    # Read volume from the CONTCAR file
    with open(com_folder + "/CONTCAR", 'r') as f:
        lines = f.readlines()
        # Read lattice vectors
        a, b, c = [np.array(wash(lines[i]), dtype="double") for i in range(2, 5)]
        volume = np.dot(np.cross(a, b), c)

    return energy, volume
