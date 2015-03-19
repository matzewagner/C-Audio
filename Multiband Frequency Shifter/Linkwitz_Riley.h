//Multiband Frequency Shifter
//written by Matthias Wagner March 2015, for MAT 240B

#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265
#endif

typedef struct
{
    double x[2];
    double y[2];

} LR_state;

typedef struct
{
    double a[3];
    double b[3];

} LR_coefficients;


void LP_coef(int sr, double freq, LR_coefficients *coef);

void HP_coef(int sr, double freq, LR_coefficients *coef);

float LP_LRFilter(float in, LR_coefficients *coef, LR_state *s1, LR_state *s2);

float HP_LRFilter(float in, LR_coefficients *coef, LR_state *s1, LR_state *s2);

