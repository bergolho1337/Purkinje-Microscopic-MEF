# ==========================================================================================================
# Script para execucao automatica do SteadyState de um conjunto de malhas.
# ==========================================================================================================

#!/bin/bash

# Variables
PROGRAM_NAME="steadyState"                  # Nome do programa
ARGS="0.1 5000"                             # Argumentos do programa
MIN_ELEM=150                                # Numero minimo de elementos
MAX_ELEM=300                                # Numero maximo de elementos
MIN_MESH=1                                  # Numero minimo de malhas
MAX_MESH=3                                  # Numero maximo de malhas
FIBER_SIZE=1                                # Tamanho da fibra

echo "======= RUNNING STEADY STATE SIMULATION ======="
# Compile the source code if doesn't exist
if [ ! -f $PROGRAM_NAME ]; then
    echo "-----------------------------------------------------------------------------------------------------"
    echo "---> Compile source code ..."
    make clean
    make
    echo "-----------------------------------------------------------------------------------------------------"
fi

# Check if the SteadyState directory already exist, if not make the directory
if [ ! -d SteadyState ]; then
    echo "Making SteadyState folder ..."
    mkdir SteadyState
fi

echo "-----------------------------------------------------------------------------------------------------"
# For each NElem
for i in $(seq $MIN_ELEM 150 $MAX_ELEM); do
    echo ">>>>>>>>>>>>>>>>>>>>>>>>> Number of elements $i <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    mkdir ./SteadyState/e-$i
    for j in $(seq $MIN_MESH $MAX_MESH); do
        echo "------ SteadyState - Mesh $j ---------"
        ./$PROGRAM_NAME $ARGS Malhas/$FIBER_SIZE/E_$i/test$j.msh $j
    done
    # Move the data to the SteadyState folder
    mv ./steadystate* ./SteadyState/e-$i/
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
done
echo "-----------------------------------------------------------------------------------------------------"


exit 0    
