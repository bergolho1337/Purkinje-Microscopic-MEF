#ifndef MONODOMAINFEM_H_
#define MONODOMAINFEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "../include/noble.h"
//#include "../include/fitz.h"
#include "../include/plot.h"
#include "../include/linSolver.h"

using namespace std;

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
//#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#define LU                                    // Metodo resolvedor do sistema linear

// Function pointer
typedef double (*Func) (int point, double t, double vm, double m, double h, double n);     

/* ============================== CONSTANTES ================================================== */
const double BETA = 0.14;                     // Razao area superficial por volume (cm^-1) (Bruno)
const double Cm = 1.0;                        // Capacitancia da membrana (uF/cm^2) (Bruno)
const double SIGMA = 0.004;                  // Condutividade citoplasmatica da celula (mS/cm) (default=0.004)
/* ============================================================================================ */

struct MonodomainFEM;
struct Point;
struct Bifurcation;

// Estrutura do resolvedor da equacao do Monodominio
struct MonodomainFEM
{
  int id;                                     // Identificador da malha sendo utilizada
  int nElem;                                  // Numero de elementos a serem utilizados
  int nPoints;                                // Numero de pontos do dominio
  int M;                                      // Numero de subintervalos no tempo
  int *map;                                   // Mapeamento dos elementos locais para os globais
  double dx;                                  // Tamanho da discretizacao no espaco (h)
  double dt;                                  // Tamanho da discretizacao no tempo (k)
  double t_max;                               // Tempo maximo de simulacao
  double x_max;                               // Tamanho maximo do dominio
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  Point *points;                              // Vetor com os pontos do dominio
  double *K;                                  // Matriz global do sistema linear
  double *A;                                  // Matriz massa global
  double *B;                                  // Matriz rigidez global
  double *F;                                  // Vetor de carga global
  double *VNew;                               // Vetor com o valor do potencial transmembranico de cada ponto no tempo n
  double *VOld;                               // Vetor com o valor do potencial transmembranico de cada ponto no tempo n-1
  double *Vstar;                              // Vetor com o valor do potencial transmembranico de cada ponto no tempo intermediario *
  double *mNew;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  double *mOld;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  double *hNew;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  double *hOld;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  double *nNew;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  double *nOld;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  vector<Bifurcation> bif;                    // Vetor das bifurcacoes
  char filename[30];                          // Nome do arquivo do steady-state que depende do id da malha
}typedef MonodomainFEM;

// Estrutura de um ponto do dominio
struct Point
{
  double x, y, z;
}typedef Point; 

// Estrutura de uma bifurcacao
struct Bifurcation
{
  int id;                                     // Identificador do Node da bifurcacao
  vector<int> links;                          // Vetor com os indices dos Node que estao ligados
}typedef Bifurcation;

/* ================================= FUNCTIONS ======================================================= */
MonodomainFEM* newMonodomainFEM (int argc, char *argv[]);
void freeMonodomain (MonodomainFEM *monoFEM);
Func* buildFunctions ();
void printInfoModel (MonodomainFEM *monoFEM);
void assembleMatrix (MonodomainFEM *monoFEM);
double* buildLocalMassMatrix (double h);
double* buildLocalStiffMatrix (double h);
double* buildGlobalMatrixFromLocal (double *local_A, int *map, int np, int ne);
double* buildGlobalMatrix (double *A, double *B, double dt, int np);
void setBoundaryConditions (double *K, int np);
void scaleFactor (double *V, double scale, int np);
void setInitialConditionsModel (MonodomainFEM *monoFEM, int argc, char *argv[]);
void setInitialConditionsModel_FromFile (MonodomainFEM *monoFEM, char *filename);
void calcPropagationVelocity (double *V, double t);
void findBifurcation (MonodomainFEM *monoFEM);
void kirchoffCondition_Matrix (MonodomainFEM *monoFEM);
void kirchoffCondition_Vector (MonodomainFEM *monoFEM);
void setVelocityPoints (double dx, int p1, int p2);

void solveMonodomain (MonodomainFEM *monoFEM);
void assembleLoadVector (MonodomainFEM *monoFEM);
void solveEDO (MonodomainFEM *monoFEM, double t);

void writeVTKFile (double *Vm, Point *points, int *map, int np, int ne, int k);
void writeSteadyStateFile (FILE *steadyFile, int nPoints, double vm[], double m[], double h[], double n[]);
void printError (char *msg);


/* =================================================================================================== */

#endif