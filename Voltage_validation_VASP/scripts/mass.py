def Relative_mass(formula):

    molar_mass = 0

    elements = []
    quantities = []
    current_element = ""
    current_quantity = ""
    for char in formula:
        if char.isupper():
        
            if current_element:
                elements.append(current_element)
                quantities.append(int(current_quantity) if current_quantity else 1)
            current_element = char
            current_quantity = ""
        elif char.islower():
           
            current_element += char
        elif char.isdigit():
            
            current_quantity += char
   
    elements.append(current_element)
    quantities.append(int(current_quantity) if current_quantity else 1)

   
    atomic_weights = {
    "H": 1.008,
    "He": 4.003,
    "Li": 6.941,
    "Be": 9.012,
    "B": 10.81,
    "C": 12.01,
    "N": 14.01,
    "O": 16.00,
    "F": 19.00,
    "Ne": 20.18,
    "Na": 22.99,
    "Mg": 24.31,
    "Al": 26.98,
    "Si": 28.09,
    "P": 30.97,
    "S": 32.07,
    "Cl": 35.45,
    "K": 39.10,
    "Ca": 40.08,
    "Sc": 44.96,
    "Ti": 47.87,
    "V": 50.94,
    "Cr": 52.00,
    "Mn": 54.94,
    "Fe": 55.85,
    "Ni": 58.69,
    "Co": 58.93,
    "Cu": 63.55,
    "Zn": 65.38,
    "Ga": 69.72,
    "Ge": 72.63,
    "As": 74.92,
    "Se": 78.96,
    "Br": 79.90,
    "Kr": 83.80,
    "Rb": 85.47,
    "Sr": 87.62,
    "Y": 88.91,
    "Zr": 91.22,
    "Nb": 92.91,
    "Mo": 95.94,
    "Tc": 98.00,
    "Ru": 101.1,
    "Rh": 102.9,
    "Pd": 106.4,
    "Ag": 107.9,
    "Cd": 112.4,
    "In": 114.8,
    "Sn": 118.7,
    "Sb": 121.8,
    "I": 126.9,
    "Te": 127.6,
    "Xe": 131.3,
    "Cs": 132.9,
    "Ba": 137.3,
    "La": 138.9,
    "Ce": 140.1,
    "Pr": 140.9,
    "Nd": 144.2,
    "Pm": 145.0,
    "Sm": 150.4,
    "Eu": 152.0,
    "Gd": 157.3,
    "Tb": 158.9,
    "Dy": 162.5,
    "Ho": 164.9,
    "Er": 167.3,
    "Tm": 168.9,
    "Yb": 173.0,
    "Lu": 175.0,
    "Hf": 178.5,
    "Ta": 180.9,
    "W": 183.8,
    "Re": 186.2,
    "Os": 190.2,
    "Ir": 192.2,
    "Pt": 195.1,
    "Au": 197.0,
    "Hg": 200.6,
    "Tl": 204.4,
    "Pb": 207.2,
    "Bi": 208.9,
    "Th": 232.0,
    "Pa": 231.0,
    "U": 238.0,
    "Np": 237.0,
    "Pu": 244.0,
    "Am": 243.0,
    "Cm": 247.0,
    "Bk": 247.0,
    "Cf": 251.0,
    "Es": 252.0,
    "Fm": 257.0,
    "Md": 258.0,
    "No": 259.0,
    "Lr": 262.0,
    "Rf": 267.0,
    "Db": 270.0,
    "Sg": 271.0,
    "Bh": 270.0,
    "Hs": 277.0,
    "Mt": 276.0,
    "Ds": 281.0,
    "Rg": 280.0,
    "Cn": 285.0,
    "Nh": 284.0,
    "Fl": 289.0,
    "Mc": 288.0,
    "Lv": 293.0,
    "Ts": 294.0,
    "Og": 294.0,
    }


    for i in range(len(elements)):
        element_mass = atomic_weights[elements[i]]
        molar_mass += element_mass * quantities[i]

    return molar_mass


