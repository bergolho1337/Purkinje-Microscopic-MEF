#include "../include/plot.h"

// Imprime uma matriz NxN
void printMatrix (char *str, double *A, int N)
{
  int i, j;
  printf("\n%s\n",str);
  for (i = 0; i < N; i++)
  {
    printf("\n");
    for (j = 0; j < N; j++)
    {
      printf("%e ",A[i * N + j]);
    }
  }
  printf("\n");
}

// Imprime um vetor Nx1
void printVector (char *str, double *b, int N)
{
  int i;
  printf("\n%s\n",str);
  for (i = 0; i < N; i++)
    printf("%e\n",b[i]);
  printf("\n");
}
