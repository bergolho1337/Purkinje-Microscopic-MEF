#include "../include/noble.h"

void setInitialConditions__Nob (double *y, int num_eq)
{
    y[0] = v0__Nob;
    y[1] = m0__Nob;
    y[2] = h0__Nob;
    y[3] = n0__Nob;
}

void setFunctions__Nob (Func *f, int num_eq)
{
    f[0] = dvdt__Nob;
    f[1] = dmdt__Nob;
    f[2] = dhdt__Nob;
    f[3] = dndt__Nob;
}

double I_Stim__Nob (int point, double t)
{
    // Somente os 2 primeiros pontos ficam como celulas de estimulo
    if (point <= 5)
    {
        double min_time, max_time;
        for (int k = 0; k < 20; k++)
        {
            min_time = k*cycle_length;
            max_time = k*cycle_length + 2.0;
            // Somente nesse periodo de tempo que o estimulo ira ocorrer
            if (t >= min_time && t < max_time)
                return v_stim__Nob;
        }
        return 0;
    }
    else
        return 0;
}

double I_Leak__Nob (double t, double vm, double m, double h, double n)
{
	return ((g_L*(vm-E_L)));
}

// Oscilacoes ainda presentes podem ser eliminadas somando um termo 1.0e+02
double g_K2__Nob (double t, double vm, double m, double h, double n)
{
    return ((1.2*pow(n,4.0e+00)));
	//return ((1.2e+03*pow(y[3],4.0e+00))+1.0e-01);          // Sem oscilacoes
}

// Trocar o primeiro coeficiente de 1.2e+03 para 1.3e+03 tambem tira o caracter auto-oscilatorio
double g_K1__Nob (double t, double vm, double m, double h, double n)
{
	return (((1.2*exp((((-vm)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((vm+9.0e+01)/6.0e+01)))));
  //return (((1.3*exp((((-y[0])-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((y[0]+9.0e+01)/6.0e+01)))));    // Sem oscilacoes
}

double I_K__Nob (double t, double vm, double m, double h, double n)
{
    return (((g_K1__Nob(t,vm,m,h,n)+g_K2__Nob(t,vm,m,h,n))*(vm+1.0e+02)));
}

double g_Na__Nob (double t, double vm, double m, double h, double n)
{
    return ((pow(m,3.0e+00)*h*g_Na_max));
}

// Mudando o coeficiente de 1.4e+02 para 1.225e+02 pode-se eliminar o comportamento auto-oscilatorio
double I_Na__Nob (double t, double vm, double m, double h, double n)
{
    return ((g_Na__Nob(t,vm,m,h,n)+1.4e-01)*(vm-E_Na));
    //return ((g_Na__Nob(t,vm,m,h,n)+1.220e-01)*(vm-E_Na));      // Sem oscilacoes
}

double dvdt__Nob (int point, double t, double vm, double m, double h, double n)
{
	//return ((-(I_Na(t,y)+I_K(t,y)+I_Leak(t,y))/CM));
    return ((-(I_Na__Nob(t,vm,m,h,n)+I_K__Nob(t,vm,m,h,n)+I_Leak__Nob(t,vm,m,h,n))+I_Stim__Nob(point,t))/CM);
}

/* ====================================================================================================== */
// dm/dt

double beta_m__Nob (double t, double vm, double m, double h, double n)
{
    return (((1.2e-01*(vm+8.0e+00))/(exp(((vm+8.0e+00)/5.0e+00))-1.0e+00)));
}

double alpha_m__Nob (double t, double vm, double m, double h, double n)
{
	return (((1.0e-01*((-vm)-4.8e+01))/(exp((((-vm)-4.8e+01)/1.5e+01))-1.0e+00)));
}

double dmdt__Nob (int point, double t, double vm, double m, double h, double n)
{
	return ((alpha_m__Nob(t,vm,m,h,n)*(1.0e+00-m))-(beta_m__Nob(t,vm,m,h,n)*m));
}

/* ====================================================================================================== */
// dh/dt

double beta_h__Nob (double t, double vm, double m, double h, double n)
{
	return ((1.0/(1.0e+00+exp((((-vm)-4.2e+01)/1.0e+01)))));
}

double alpha_h__Nob (double t, double vm, double m, double h, double n)
{
	return ((1.7e-01*exp((((-vm)-9.0e+01)/2.0e+01))));
}

double dhdt__Nob (int point, double t, double vm, double m, double h, double n)
{
	return ((alpha_h__Nob(t,vm,m,h,n)*(1.0e+00-h))-(beta_h__Nob(t,vm,m,h,n)*h));
}

/* ====================================================================================================== */
// dn/dt

double beta_n__Nob (double t, double vm, double m, double h, double n)
{
	return ((2.0e-03*exp((((-vm)-9.0e+01)/8.0e+01))));
}

double alpha_n__Nob (double t, double vm, double m, double h, double n)
{
	return (((1.0e-04*((-vm)-5.0e+01))/(exp((((-vm)-5.0e+01)/1.0e+01))-1.0e+00)));
}

double dndt__Nob (int point, double t, double vm, double m, double h, double n)
{
	return ((alpha_n__Nob(t,vm,m,h,n)*(1.0e+00-n))-(beta_n__Nob(t,vm,m,h,n)*n));
}

/* ====================================================================================================== */