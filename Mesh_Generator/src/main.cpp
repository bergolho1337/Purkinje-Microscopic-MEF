#include <cstdio>
#include <cstdlib>
#include "../include/mesh.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 < 3)
    {
        printf("============================================================================================\n");
        printf("Gera uma malha de elementos finitos a partir em um arquivo .vtk com o esqueleto da mesma.\n");
        printf("--------------------------------------------------------------------------------------------\n");
        printf("Usage:> %s <Xmax> <nElem> <in_VTK_file>\n",argv[0]);
        printf("<Xmax> = Tamanho de uma fibra\n");
        printf("<nElem> = Numero de elementos em uma fibra\n");
        printf("<in_VTK_file> = Nome do arquivo de entrada .vtk (saida --> <in_VTK_file>.txt)\n");
        printf("Try for example:> %s 2.0 50.0 ./Purkinje_Networks/tree1.vtk\n",argv[0]);
        printf("============================================================================================\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("============================================================================================\n");
        Mesh *mesh = newMesh(argc,argv);
        writeMeshToFile(mesh,argv[3]);
        printf("[+] Done\n");
        printf("============================================================================================\n");
        return 0;
    }
}