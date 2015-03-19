//Multiband Frequency Shifter
//written by Matthias Wagner March 2015, for MAT 240B

#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265
#endif

typedef struct
{
    double xnm1[12];
    double ynm1[12];
    double coef[12];

} H_coefficients;

void Hil_coef(int sr, H_coefficients *hc);

float Hil_real(float in, H_coefficients *hc);

float Hil_imaginary(float in, H_coefficients *hc);

