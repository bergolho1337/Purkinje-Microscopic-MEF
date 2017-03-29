#include "../include/monodomainFEM.h"

bool stage1 = true;
bool stage2 = false;
int id_1, id_2;
double t1, t2;
double delta_x;

// Construtor da estrutura do resolvedor da equacao do monodominio
MonodomainFEM* newMonodomainFEM (int argc, char *argv[])
{
    MonodomainFEM *monoFEM = (MonodomainFEM*)malloc(sizeof(MonodomainFEM));
    monoFEM->dt = atof(argv[1]);
    monoFEM->t_max = atof(argv[2]);
    monoFEM->M = nearbyint(monoFEM->t_max / monoFEM->dt);
    // Le o arquivo de malha e retirar as informacoes
    int left, right;
    FILE *file = fopen(argv[3],"r");
    if (!fscanf(file,"%d %d %lf",&monoFEM->nElem,&monoFEM->nPoints,&monoFEM->dx)) printError("Reading file");
    // Lendo os pontos
    monoFEM->points = (Point*)malloc(sizeof(Point)*monoFEM->nPoints);
    for (int i = 0; i < monoFEM->nPoints; i++)
      if (!fscanf(file,"%lf %lf %lf",&monoFEM->points[i].x,&monoFEM->points[i].y,&monoFEM->points[i].z)) printError("Reading file");
    // Lendo os elementos
    monoFEM->map = (int*)calloc(monoFEM->nElem*4,sizeof(int));
    for (int i = 0; i < monoFEM->nElem; i++)
    {
      if (!fscanf(file,"%d %d",&left,&right)) printError("Reading file");
      monoFEM->map[i*4] = left;
      monoFEM->map[i*4+1] = right;
      monoFEM->map[i*4+2] = left + monoFEM->nPoints;
      monoFEM->map[i*4+3] = right + monoFEM->nPoints;
    }
    monoFEM->functions = buildFunctions();
    
    // Alocar memoria e setar as condicoes iniciais (para o Fitz-Hugh Nagumo eh tudo 0)
    // ** Lembrando que o elemento de Hermite possui 2*N o numero de pontos
    // Ou seja, o potencial Vm --> [0,N-1], I_L = [N,2N-1]
    monoFEM->VOld = (double*)calloc(monoFEM->nPoints*2,sizeof(double));
    monoFEM->Vstar = (double*)calloc(monoFEM->nPoints*2,sizeof(double));
    monoFEM->VNew = (double*)calloc(monoFEM->nPoints*2,sizeof(double));
    monoFEM->wOld = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->wNew = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->F = (double*)calloc(monoFEM->nPoints*2,sizeof(double));

    // Construir a matriz global do sistema linear ligado a solucao da EDP
    assembleMatrix(monoFEM);

    // Tratar a condicao de conservacao da corrente nas bifurcacoes na matriz K
    kirchoffCondition_Matrix(monoFEM);

    // Decompor a matriz em LU
    LUDecomposition(monoFEM->K,monoFEM->nPoints*2);

    #ifdef DEBUG
    printInfoModel(monoFEM);
    #endif

    return monoFEM;
}

// Imprime informacoes sobre o modelo construido
void printInfoModel (MonodomainFEM *monoFEM)
{
    int i;
    printf("======================== INFO MODEL ============================\n");
    printf("Number of elements = %d\n",monoFEM->nElem);
    printf("Number of points = %d\n",monoFEM->nPoints);
    printf("dt = %e\n",monoFEM->dt);
    printf("t_max = %e\n",monoFEM->t_max);
    printf("dx = %e\n",monoFEM->dx);
    printf("M = %d\n",monoFEM->M);
    printf("------------------------- MAPPING ------------------------------\n");
    for (i = 0; i < monoFEM->nElem; i++)
    printf("%d %d %d %d\n",monoFEM->map[i*4],monoFEM->map[i*4+1],monoFEM->map[i*4+2],monoFEM->map[i*4+3]);
    printf("------------------------- POINTS -------------------------------\n");
    for (i = 0; i < monoFEM->nPoints; i++)
    printf("%d - %lf %lf %lf\n",i,monoFEM->points[i].x,monoFEM->points[i].y,monoFEM->points[i].z);
    printf("================================================================\n");
}

void setInitialConditionsModel (MonodomainFEM *monoFEM, double v0, double w0)
{
  int i, np;
  np = monoFEM->nPoints;
  for (i = 0; i < np; i++)
  {
    monoFEM->VOld[i] = v0;
    monoFEM->wOld[i] = w0;
  }
}

// Constroi o vetor de funcoes que serao utilizados no programa
// 0 = Funcao do potencial transmembranico (Vm)
// 1 = Funcao da variavel de estado (w)
Func* buildFunctions ()
{
  Func *func = (Func*)malloc(sizeof(Func)*num_eq);
  func[0] = dvdt__Fitz;
  func[1] = dwdt__Fitz;
  return func;
}

// Constroi a matriz global K do metodo implicito da solucao da EDP pelo MEF
// K = (BETA*Cm*A + SIGMA*dt*B)
void assembleMatrix (MonodomainFEM *monoFEM)
{
  printf("[!] Construindo matriz ... ");
  fflush(stdout);
  int np, ne;
  double *local_A, *local_B;
  ne = monoFEM->nElem;
  np = monoFEM->nPoints;

  // Construir as matrizes locais de massa e de rigidez 
  local_A = buildLocalMassMatrix(monoFEM->dx);
  local_B = buildLocalStiffMatrix(monoFEM->dx);

  // Construir as matrizes globais A (massa) e B (rigidez)
  monoFEM->A = buildGlobalMatrixFromLocal(local_A,monoFEM->map,np,ne);
  monoFEM->B = buildGlobalMatrixFromLocal(local_B,monoFEM->map,np,ne);

  // Construir a matriz global do sistema linear: 
  monoFEM->K = buildGlobalMatrix(monoFEM->A,monoFEM->B,monoFEM->dt,np);

  // Setar as condicoes de contorno, a corrente i_L deve ser 0 em x=0 e x=L ?????
  //setBoundaryConditions(monoFEM->K,np*2);

  printf("ok\n");  
  #ifdef DEBUG
  printMatrix("Global mass matrix A",monoFEM->A,np*2);
  printMatrix("Global stiff matrix B",monoFEM->B,np*2);
  printMatrix("Global matrix K",monoFEM->K,np*2);
  #endif

  // !!! As condicoes de contorno saem naturalmente da formulacao variacional !!!
  // Decompor a matriz K em LU
  //LUDecomposition(monoFEM->K,monoFEM->nPoints);

  //#ifdef DEBUG
  //printMatrix("Global matrix K after LU decomposition",monoFEM->K,np*2);
  //#endif

  // Liberar memoria
  free(local_A);
  free(local_B);

}

// Constroi a matriz local de massa de cada elemento (4x4)
// !!! USANDO ELEMENTOS HERMITE !!!
// || phi_1 = 2*(x/h)**3 - 3*(x/h)**2 + 1 || phi_2 = 3*(x/h)**2 - 2*(x/h)**3 ||
// || phi_3 = (x/h)**3 - 2*(x/h)**2 + (x/h) || phi_4 = (x/h)**3 - (x/h)**2 
// --> integral_(0,h) {phi_i . phi_j dx}
double* buildLocalMassMatrix (double h)
{
  double *A = (double*)calloc(16,sizeof(double));
  double r = h / 35.0;
  A[0*4] = r*13.0;
  A[1*4+1] = r*13.0;
  A[2*4+2] = r*(1.0/3.0);
  A[3*4+3] = r*(1.0/3.0);
  A[0*4+1] = A[1*4] = r*(9.0/2.0);
  A[0*4+2] = A[2*4] = r*(11.0/6.0);
  A[0*4+3] = A[3*4] = r*(-13.0/12.0);
  A[1*4+2] = A[2*4+1] = r*(13.0/12.0);
  A[1*4+3] = A[3*4+1] = r*(-11.0/6.0);
  A[2*4+3] = A[3*4+2] = r*(-1.0/4.0);

  return A;
}

// Constroi a matriz local de rigidez de cada elemento (4x4)
// !!! USANDO ELEMENTOS HERMITE !!!
// || phi_1' = 6*(x**2/h**3) - 6*(x/h**2) || phi_2' = 6*(x/h**2) - 6*(x**2/h**3) ||
// || phi_3' = 3*(x**2/h**3) - 4*(x/h**2) + (1/h) || phi_4' = 3*(x**2/h**3) - 2*(x/h**2) 
// --> integral_(0,h) {phi_i' . phi_j' dx}
double* buildLocalStiffMatrix (double h)
{
  double *B = (double*)calloc(16,sizeof(double));
  double r = 1.0 / (30.0*h);
  B[0*4] = r*36.0;
  B[1*4+1] = r*36.0;
  B[2*4+2] = r*4.0;
  B[3*4+3] = r*4.0;
  B[0*4+1] = B[1*4] = -36.0*r;
  B[0*4+2] = B[2*4] = r*3.0;
  B[0*4+3] = B[3*4] = r*3.0;
  B[1*4+2] = B[2*4+1] = -3.0*r;
  B[1*4+3] = B[3*4+1] = -3.0*r;
  B[2*4+3] = B[3*4+2] = -1.0*r;
  return B;
}

// Constroi uma matriz global a partir da matriz local de cada elemento
double* buildGlobalMatrixFromLocal (double *local_M, int *map, int np, int ne)
{
  int i, j, k;
  int ig, jg;
  // Hermite (2N)x(2N)
  np *= 2;
  double *M = (double*)calloc(np*np,sizeof(double));
  // Para cada elemento
  for (k = 0; k < ne; k++)
  {
    // Percorrer a matriz local
    for (i = 0; i < 4; i++)
    {
      ig = map[k*4+i];
      for (j = 0; j < 4; j++)
      {
        jg = map[k*4+j];
        // Juntar a contribuicao local do elemento na matriz global
        M[ig*np+jg] += local_M[i*4+j];
      }
    }
  }
  return M;
}

// Constroi a matriz global do sistema linear da EDP
// K = (BETA*Cm*A + SIGMA*dt*B)
double* buildGlobalMatrix (double *A, double *B, double dt, int np)
{
  int i, j;
  // Elemento de Hermite
  np *= 2;
  double *K = (double*)calloc(np*np,sizeof(double));
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
      K[i*np+j] = (BETA*Cm*A[i*np+j]) + (SIGMA*dt*B[i*np+j]);
  }
  return K;
}

// Setar as condicoes de contorno de Neumann na matriz
void setBoundaryConditions (double *K, int np)
{
  // Pelo metodo da penalizacao
  double epsilon = 1.0e+08;
  K[(np/2)*np+(np/2)] += epsilon;
  K[(np-1)*np+(np-1)] += epsilon;
}

// Multiplicar os valores da derivada pelo fator de escala para corrigi-los
void scaleFactor (double *V, double scale, int np)
{
  int i;
  for (i = np; i < 2*np; i++)
    V[i] *= scale;
}

// Constroi o vetor de carga do sistema linear relacionado a EDP
// F = (BETA*Cm*A*VOld)
void assembleLoadVector (MonodomainFEM *monoFEM)
{
  int i, j;
  int np;
  // Elemento de Hermite possui 2*N
  np = monoFEM->nPoints*2;
  // Iniciliza o vetor global com 0's
  memset(monoFEM->F,0,sizeof(double)*np);
  // Montar o vetor global a partir da multiplicacao A*VOld
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
      monoFEM->F[i] += monoFEM->A[i*np+j]*monoFEM->VOld[j];
    // Multiplicar pelo coeficiente BETA*Cm
    monoFEM->F[i] *= BETA*Cm;
  }

  #ifdef DEBUG
  printVector("Load vector F",monoFEM->F,np);
  #endif
}

// Resolver o sistema nao linear de EDOs no tempo atual
void solveEDO (MonodomainFEM *monoFEM, double t)
{
  int np;
  int i, point;
  double f, dt;
  np = monoFEM->nPoints;
  dt = monoFEM->dt;
  // Resolver o sistema de EDO para cada ponto (Potencial e as variaveis de estado)
  for (point = 0; point < np; point++)
  {
    for (i = 0; i < num_eq; i++)
    {
      // V^{n+1} = V^{*} + f*dt
      if (i == 0)
      {
        f = monoFEM->functions[i](point,t,monoFEM->Vstar[point],monoFEM->wOld[point]);
        monoFEM->VNew[point] = monoFEM->Vstar[point] + f*dt;
      }
      // w^{n+1} = w^{n} + f*dt
      else
      {
        f = monoFEM->functions[i](point,t,monoFEM->VOld[point],monoFEM->wOld[point]);
        monoFEM->wNew[point] = monoFEM->wOld[point] + f*dt;
      }     
    }
  }
  // Resolve agora para a corrente i_L^(n-1) = i_L^(n)  
  for (point = np; point < 2*np; point++)
    monoFEM->VNew[point] = monoFEM->Vstar[point]; 
  
}


// Resolve a equacao do monodominio
void solveMonodomain (MonodomainFEM *monoFEM)
{
  double t;
  // A matriz global do problema jah se encontra como LU
  printf("[!] Resolvendo o problema transiente ... ");
  fflush(stdout);
  int i, M;
  M = monoFEM->M;
  // Iterar o metodo a cada passo de tempo
  for (i = 0; i < M; i++)
  {
    t = i*monoFEM->dt;

    // Escrever em arquivo .vtk
    writeVTKFile(monoFEM->VOld,monoFEM->points,monoFEM->map,monoFEM->nPoints,monoFEM->nElem,i);

    // Resolver a EDP (parte difusiva)
    assembleLoadVector(monoFEM);
    kirchoffCondition_Vector(monoFEM);
    solveLinearSystem_LU(monoFEM->K,monoFEM->F,monoFEM->Vstar,monoFEM->nPoints*2);
    //solveLinearSystem_CG(monoFEM->K,monoFEM->F,monoFEM->Vstar,monoFEM->nPoints*2);

    // Multiplicar os termos da derivada pelo fator de escala
    scaleFactor(monoFEM->Vstar,SIGMA/monoFEM->dx,monoFEM->nPoints);

    // Resolver as EDOs (parte reativa)
    solveEDO(monoFEM,t);

    calcPropagationVelocity(monoFEM->VNew,t);

    #ifdef DEBUG
    printVector("Vstar",monoFEM->Vstar,monoFEM->nPoints*2);
    printVector("VNew",monoFEM->VNew,monoFEM->nPoints*2);
    #endif

    // Passa para a proxima iteracao
    memcpy(monoFEM->VOld,monoFEM->VNew,sizeof(double)*monoFEM->nPoints*2);
    memcpy(monoFEM->wOld,monoFEM->wNew,sizeof(double)*monoFEM->nPoints);
  } 
  //printf("ok\n");
}

// Calcula a velocidade de propagacao entre dois pontos
void calcPropagationVelocity (double *V, double t)
{
  if (stage1)
  {
    // Potencial do ponto 1 chegou no ponto minimo ?
    if (V[id_1] < -1.90)
    {
      t1 = t;
      stage2 = true;
      stage1 = false;
      printf("\n\n[!] Propagation velocity! Stage 1 clear!\n");
      printf("t1 = %e\n",t1);
      printf("V[%d] = %e\n",id_1,V[id_1]);
    }
  }
  else if (stage2)
  {
    // Potencial do ponto 2 chegou no ponto minimo ?
    if (V[id_2] < -1.90)
    {
      double velocity;
      t2 = t;
      stage2 = false;
      // Calcula velocidade: v = dx/dt
      velocity = delta_x / (t2-t1);
      printf("\n[!] Propagation velocity! Stage 2 clear!\n");
      printf("t2 = %e\n",t2);
      printf("V[%d] = %e\n",id_2,V[id_2]);
      printf("\n!!!!!!!! Propagation velocity = %e cm/s !!!!!!!!!!\n",velocity);
    }
  }
}

// Aplica a condicao de Kirchoff (conservacao de corrente) nas bifurcacoes na parte da matriz
void kirchoffCondition_Matrix (MonodomainFEM *monoFEM)
{
  int np = monoFEM->nPoints*2;
  findBifurcation(monoFEM);
  // Percorrer as bifurcacoes colocando a condicao: i_A - i_B - i_C = 0 (ver Vigmond)
  for (int i = 0; i < (int)monoFEM->bif.size(); i++)
  {
    int lin = monoFEM->bif[i].id + monoFEM->nPoints;
    // Zerar a linha
    for (int j = 0; j < np; j++) monoFEM->K[lin*np+j] = 0;
    // Aplicar a condicao de Kirchoff sobre esse Node
    monoFEM->K[lin*np+lin] = 1;
    for (int j = 0; j < (int)monoFEM->bif[i].links.size(); j++)
    {
      int col = monoFEM->bif[i].links[j] + monoFEM->nPoints;
      monoFEM->K[lin*np+col] = -1;
    }
  }
}

// Descobre aonde estao as bifurcacoes
void findBifurcation (MonodomainFEM *monoFEM)
{
  int cont;
  // Descobrir quais Nodes sao bifurcacoes, verifica pelo numero de vizinhos (> 2)
  for (int i = 0; i < monoFEM->nPoints; i++)
  {
    cont = 0;
    for (int j = 0; j < monoFEM->nElem; j++)
    {
      if (i == monoFEM->map[j*4] || i == monoFEM->map[j*4+1])
        cont++;
    }
    // Eh bifurcacao ?
    if (cont > 2)
    {
      Bifurcation new_bif;
      new_bif.id = i;
      // Descobrir quais os Nodes envolvidos na bifurcacao
      for (int j = 0; j < monoFEM->nElem; j++)
      {
        if (i == monoFEM->map[j*4])
          new_bif.links.push_back(monoFEM->map[j*4+1]);
      }
      monoFEM->bif.push_back(new_bif);
    }
  }
}

// Aplica a condicao de Kirchoff (conservacao de corrente) nas bifurcacoes na parte do vetor
void kirchoffCondition_Vector (MonodomainFEM *monoFEM)
{
  for (int i = 0; i < (int)monoFEM->bif.size(); i++)
  {
    int lin = monoFEM->bif[i].id + monoFEM->nPoints;
    monoFEM->F[lin] = 0;
  }
}

void freeMonodomain (MonodomainFEM *monoFEM)
{
  printf("[!] Liberando memoria ... ");
  fflush(stdout);
  free(monoFEM->map);
  free(monoFEM->functions);
  free(monoFEM->points);
  free(monoFEM->K);
  free(monoFEM->A);
  free(monoFEM->B);
  free(monoFEM->F);
  free(monoFEM->VNew);
  free(monoFEM->VOld);
  free(monoFEM->Vstar);
  free(monoFEM->wOld);
  free(monoFEM->wNew);
  free(monoFEM);
  printf("ok\n");
}

void writeVTKFile (double *Vm, Point *points, int *map, int np, int ne, int k)
{
  FILE *file;
  int i;
  char filename[50];
  // Escrever o potencial transmembranico
  sprintf(filename,"VTK/Potential/potential%d.vtk",k);
  file = fopen(filename,"w+");
  fprintf(file,"# vtk DataFile Version 3.0\n");
  fprintf(file,"Monodomain FEM\n");
  fprintf(file,"ASCII\n");
  fprintf(file,"DATASET POLYDATA\n");
  fprintf(file,"POINTS %d float\n",np);
  for (i = 0; i < np; i++)
    fprintf(file,"%e %e %e\n",points[i].x,points[i].y,points[i].z);
  fprintf(file,"LINES %d %d\n",np,ne*3);
  for (i = 0; i < ne; i++)
  {
    fprintf(file,"2 %d %d\n",map[i*4],map[i*4+1]);
  }
    
  fprintf(file,"POINT_DATA %d\n",np);
  fprintf(file,"SCALARS vm float 1\n");
  fprintf(file,"LOOKUP_TABLE default\n");
  // Vm --> [0,N-1]
  for (i = 0; i < np; i++)
    fprintf(file,"%e\n",Vm[i]);
  fclose(file);

  // Escrever a corrente (derivada do potencial)
  sprintf(filename,"VTK/Current/current%d.vtk",k);
  file = fopen(filename,"w+");
  fprintf(file,"# vtk DataFile Version 3.0\n");
  fprintf(file,"Monodomain FEM\n");
  fprintf(file,"ASCII\n");
  fprintf(file,"DATASET POLYDATA\n");
  fprintf(file,"POINTS %d float\n",np);
  for (i = 0; i < np; i++)
    fprintf(file,"%e %e %e\n",points[i].x,points[i].y,points[i].z);
  fprintf(file,"LINES %d %d\n",np,ne*3);
  for (i = 0; i < ne; i++)
  {
    fprintf(file,"2 %d %d\n",map[i*4],map[i*4+1]);
  }
    
  fprintf(file,"POINT_DATA %d\n",np);
  fprintf(file,"SCALARS vm float 1\n");
  fprintf(file,"LOOKUP_TABLE default\n");
  // i_L --> [N,2N-1]
  for (i = np; i < 2*np; i++)
    fprintf(file,"%e\n",Vm[i]);
  fclose(file);
}

// Imprime uma mensagem de erro e sai do programa
void printError (char *msg)
{
  printf("[-] ERROR ! %s\n",msg);
  exit(1);
}