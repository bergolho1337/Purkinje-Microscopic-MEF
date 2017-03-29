#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*Func) (int elem, double t, double vm, double w);

void printMatrix (char *str, double *A, int N);
void printVector (char *str, double *b, int N);


#endif