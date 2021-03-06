# ================================================================================================================
# Script para solucao automatica de um conjunto de malhas
# Cada simulacao eh feita para um tamanho fixo de fibra (FIBER_SIZE) e de elementos (NELEM)
# ================================================================================================================
#!/bin/bash

PROGRAM_NAME="purkinjeFEM"
ARGS="0.1 1000"
MIN_MESH=1                                  # Numero minimo de malhas
MAX_MESH=5                                  # Numero maximo de malhas
FIBER_SIZE=0.50                                # Tamanho da fibra
NELEM=75                                   # Numero de elementos

echo "======= RUNNING BIFURCATION SIMULATION ======="
echo "FIBER SIZE = $FIBER_SIZE cm"

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
for i in $(seq $MIN_MESH $MAX_MESH); do
    echo "------ Simulation $i ---------"
    ./$PROGRAM_NAME $ARGS Malhas/$FIBER_SIZE/E_$NELEM/test$i.msh SteadyState/e-$NELEM/steadystate$i.dat
    # Copy the results in folder VTK to correct one in the Results folder
    cp -r VTK Resultados/
    cp ./velocity.txt Resultados/VTK/velocity$i.txt
    cp ./block.txt Resultados/VTK/block$i.txt
    # Rename it to the appropriate mesh name
    mv Resultados/VTK Resultados/Mesh_$i
    echo 
done
echo "-----------------------------------------------------------------------------------------------------"
    




