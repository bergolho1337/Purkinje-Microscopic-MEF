#!/bin/bash

# Constants
PROGRAM_NAME="purkinjeFEM"
ARGS="0.1 5000"
BCL=250
NELEM=50

C_MESH=1
C_BCL=1

MAX_NELEM=50
MAX_BCL=250
MAX_MESH=10


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