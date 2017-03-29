#include "../include/fitz.h"

void setInitialConditions__Fitz (double *y, int num_eq)
{
    y[0] = v0__Fitz;
    y[1] = w0__Fitz;
}

void setFunctions__Fitz (Func *f, int num_eq)
{
    f[0] = dvdt__Fitz;
    f[1] = dwdt__Fitz;
}

double I_stim__Fitz (int point, double t)
{
    // Somente os 3 primeiros pontos ficam como celulas de estimulo
    if (point <= 1)
    {
        // Somente nesse periodo de tempo que o estimulo ira ocorrer
        if (t >= 0 && t < 2)
            return v_stim__Fitz;
        else
            return 0;
    }
    else
        return 0;
}

// Recebe o potencial 'vm' e a variavel de estado 'w' de cada ponto
double dvdt__Fitz (int point, double t, double vm, double w)
{
    return vm - pow(vm,3)/3.0 - w + I_stim__Fitz(point,t);
}

// Recebe o potencial 'vm' e a variavel de estado 'w' de cada ponto
double dwdt__Fitz (int point, double t, double vm, double w)
{
    return (vm + a__Fitz - b__Fitz*w) / tal__Fitz;
}