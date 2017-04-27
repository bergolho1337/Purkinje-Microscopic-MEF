#!/bin/bash

PROGRAM_NAME="purkinjeFEM"
ARGS="0.1 1000"
C_TEST=1
MAX_TEST=10

echo "======= RUNNING BIFURCATION SIMULATION ======="

# Compile the source code if doesn't exist
if [ ! -f $PROGRAM_NAME ]; then
    echo "-----------------------------------------------------------------------------------------------------"
    echo "---> Compile source code ..."
    make clean
    make
    echo "-----------------------------------------------------------------------------------------------------"
fi

# Check if the Results directory already exist, if not make the directory
if [ ! -d Resultados ]; then
    echo "Making Resultados folder ..."
    mkdir Resultados
fi

# Run the simulation for each mesh
echo "-----------------------------------------------------------------------------------------------------"
while [ $C_TEST -le $MAX_TEST ]; do
    echo "------ Simulation $C_TEST ---------"
    ./$PROGRAM_NAME $ARGS Malhas/NElem_300/test$C_TEST.msh
    # Copy the results in folder VTK to correct one in the Results folder
    cp -r VTK Resultados/
    # Rename it to the appropriate mesh name
    mv Resultados/VTK Resultados/Mesh_$C_TEST     
    let C_TEST=C_TEST+1
    echo
done
echo "-----------------------------------------------------------------------------------------------------"
    




