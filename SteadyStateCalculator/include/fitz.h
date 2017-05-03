#ifndef FITZ_H_
#define FITZ_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Maximum time = 200s
// dt = 0.1

// Function pointer
typedef double (*Func) (int elem, double t, double vm, double w);

/* ============================== CONSTANTS ========================================== */
const int num_eq = 2;
//const double v0__Fitz = 0.0;
//const double w0__Fitz = 0.0;
const double v0__Fitz = -2.0;                       // Equilibrio
const double w0__Fitz = -2.0;                       // Equilibrio
const double v_stim__Fitz = 0.1;
const double a__Fitz = 0.0;                         // Auto-oscillatory
const double b__Fitz = 0.0;                         // Auto-oscillatory
const double tal__Fitz = 12.5;
/* ============================ FUNCTIONS ============================================= */
double dvdt__Fitz (int point, double t, double vm, double w);
double dwdt__Fitz (int point, double t, double vm, double w);
double I_stim__Fitz (int point, double t);
/* ==================================================================================== */

void setInitialConditions__Fitz (double *y, int num_eq);
void setFunctions__Fitz (Func *f, int num_eq);

#endif