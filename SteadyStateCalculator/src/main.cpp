/*
******************************************
    Autor: Lucas Berg
    Last Update: 27/04/17
******************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "../include/monodomainFEM.h"

int main (int argc, char *argv[])
{
  /*
  printf("========= METODO DOS ELEMENTOS FINITOS ----- EQUACAO DO MONODOMINIO (CABO) ===== v.2.0 ========\n");
  printf("Problema: resolver a equacao do cabo utilizando a equacao do monodominio aplicando o MEF.\n");
  printf("\t{ BETA*Cm*V_t = SIGMA*V_xx\n");
  printf("\t{ V'(0,t) = V'(1,t) = 0\n");
  printf("\t{ V(x,0) = 0\n");
  printf("***********************************************************************************************\n");
  printf("BETA = Razao entre area superficial por volume da celula (cm^-1)\n");
  printf("Cm = Capacitancia da membrana celular (uF/cm^2)\n");
  printf("SIGMA = Condutividade da membrana celular (mS/cm^2)\n");
  printf("***********************************************************************************************\n");
  printf("!!!!       Este exemplo utiliza o modelo celular de Noble      !!!!\n");
  printf("!!!!             Utilizando elementos de Hermite               !!!!\n");
  printf("-----------------------------------------------------------------------------------------------\n");
  */
  if (argc-1 < 5)
  {
    printf("Usage:> %s <dt> <t_max> <mesh_file> <mesh_id> <BCL>\n",argv[0]);
    printf("---------------------------------------------------------------------------------------------\n");
    printf("<dt> = Tamanho da discretizacao no tempo\n");
    printf("<t_max> = Tempo maximo de simulacao\n");
    printf("<mesh_file> = Arquivo contendo a malha com os elementos e os pontos definidos\n");
    printf("<mesh_id> = Identiicador da malha\n");
    printf("<BCL> = Ciclo do pacing\n");
    printf("[!] DEBUGACAO EH ATIVADA POR FLAG NO ARQUIVO \"monodomainFEM.h\".\n");
    printf("---------------------------------------------------------------------------------------------\n");
    printf("Exemplo: ./purkinjeFEM 0.1 10000 Malhas/test1.msh 1 500\n");
    printf("=============================================================================================\n");
    return 1;
  }
  else
  {
    double start, finish, elapsed;
    MonodomainFEM *monoFEM = newMonodomainFEM(argc,argv);
    
    start = clock();
    solveMonodomain(monoFEM);
    finish = clock();
    elapsed = (double)(finish-start) / (double)CLOCKS_PER_SEC;

    freeMonodomain(monoFEM);
    printf("==========================================================\n");
    printf("[!] Time elapsed = %.10lf\n",elapsed);
    printf("==========================================================\n");

    return 0;
  }
}
