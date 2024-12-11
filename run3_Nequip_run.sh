#!/bin/bash
cd MLMD

for dir in */
do
  cd "$dir"
  
  if [ ! -f "finished" ]
  then
    ~/opt/lammps-nequip/new/lammps-17Apr2024/build/lmp < in.lammps
  fi
  cd ..
done
