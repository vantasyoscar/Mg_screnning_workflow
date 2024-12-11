# Mg Cathode Analysis Workflow

## Overview

This workflow provides an end-to-end solution to analyze magnesium (Mg) cathode materials. By placing your Mg cathode files in the `raw_structures` directory, the workflow will calculate both the electrode voltage and the ionic conductivity of the materials.

## Features

1. **Electrode Voltage Prediction**: The prediction part of the code is a modified version of the [CGCNN](https://github.com/txie-93/cgcnn) model.
2. **Ionic Conductivity Simulation**: Molecular dynamics simulations for ionic conductivity are conducted using [NequIP](https://github.com/mir-group/nequip).

## Usage

To run the workflow, follow these steps:

1. **Predict Voltage**:
   ```bash
   python run1_voltage_select.py
   ```

2. **Generate Structures for NequIP**:
   ```bash
   python run2_Nequip_generate.py
   ```

3. **Run NequIP Simulations** (requires LAMMPS with NequIP module):
   ```bash
   sh ./run3_Nequip_run.sh
   ```

Make sure you have the necessary dependencies installed and the NequIP module integrated with LAMMPS.

## Dependencies

- [CGCNN](https://github.com/txie-93/cgcnn)
- [NequIP](https://github.com/mir-group/nequip)
- LAMMPS with NequIP module

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

We would like to thank the developers of [CGCNN](https://github.com/txie-93/cgcnn) and [NequIP](https://github.com/mir-group/nequip) for their outstanding contributions to the field of material science.

---
