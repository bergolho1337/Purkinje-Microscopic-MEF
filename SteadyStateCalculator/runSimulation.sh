# ==========================================================================================================
# Script para execucao automatica do SteadyState de um conjunto de malhas.
# ==========================================================================================================

#!/bin/bash

# Constants
PROGRAM_NAME="steadyState"                  # Nome do programa
ARGS="0.1 5000"                             # Argumentos do programa
BCL=250                                     # Basic cycle length = Ciclo basico do pacing
NELEM=50                                    # Numero de elementos a ser utilizado

# Variables
C_MESH=1                                    # Contador da malha a ser resolvida
C_BCL=1                                     # Contador do BCL a ser executado

MAX_NELEM=100                                # Numero maximo de elementos que serao testados
MAX_BCL=500                                 # Numero maximo de BCL que serao testados
MAX_MESH=5                                 # Numero maximo de malhas que serao testados


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
while [ $NELEM -le $MAX_NELEM ]; do
    echo ">>>>>>>>>>>>>>>>>>>>>>>>> Number of elements $NELEM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    let C_BCL=BCL
    # For each BCL
    while [ $C_BCL -le $MAX_BCL ]; do
        echo "................................ Cycle length $C_BCL ......................................."
        mkdir ./SteadyState/e-$NELEM-c-$C_BCL
        # For each mesh
        let C_MESH=1
        while [ $C_MESH -le $MAX_MESH ]; do
            echo "------ SteadyState - Mesh $C_MESH ---------"
            ./$PROGRAM_NAME $ARGS Malhas/NElem_$NELEM/test$C_MESH.msh $C_MESH $BCL
            let C_MESH=C_MESH+1
        echo
        # Copy the data to the SteadyState folder
        mv ./steadystate* ./SteadyState/e-$NELEM-c-$C_BCL/
        done
        let C_BCL=C_BCL*2
    done
    let NELEM=NELEM*2
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
done
echo "-----------------------------------------------------------------------------------------------------"


exit 0    
