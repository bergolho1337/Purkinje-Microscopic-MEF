#!/bin/bash

# Variables
MIN_SIZE=1
MAX_SIZE=2
C_MESH=1
MAX_MESH=10
MIN_ELEM=150
MAX_ELEM=300

echo "========================================================================================"
echo "[!] Generating meshes ..."
for i in $(seq $MIN_SIZE $MAX_SIZE); do
    #echo "----- Size = $i -----"
    mkdir ./Mesh-Networks/$i
    for j in $(seq $MIN_ELEM 50 $MAX_ELEM); do
        #echo "--------------------- Elements $j ------------------------"
        mkdir ./Mesh-Networks/$i/E_$j
        for k in $(seq $MIN_MESH $MAX_MESH); do
            ./meshGenerator $i $j ./Skeleton-Networks/$i/test$k.vtk ./Mesh-Networks/$i/E_$j/test$k.msh
            #echo "========================================= Mesh $k =========================================="
        done
    done
done

echo "========================================================================================"
exit 0