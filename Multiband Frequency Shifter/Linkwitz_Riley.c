//Multiband Frequency Shifter
//written by Matthias Wagner March 2015, for MAT 240B

#include <stdio.h>
#include <math.h>

#include "Linkwitz_Riley.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif

void LP_coef(int sr, double freq, LR_coefficients *coef) {
        double q;
        double w0;
        double alpha;

        w0 = 2.0*M_PI*freq/sr;
        double sn = sin(w0);
        double cs = cos(w0);
        q = (double)1.0/sqrt(2.0);
        alpha = sn/(2.0*q);

        coef->b[0] = (1-cs) / 2.0;
        coef->b[1] = (1-cs);
        coef->b[2] = (1-cs) / 2.0;
        coef->a[0] = 1+alpha;
        coef->a[1] = -2.0*cs;
        coef->a[2] = 1-alpha;

        coef->b[0] /= coef->a[0];
        coef->b[1] /= coef->a[0];
        coef->b[2] /= coef->a[0];
        coef->a[1] /= coef->a[0];
        coef->a[2] /= coef->a[0];
}

float LP_LRFilter(float in, LR_coefficients *coef, LR_state *s1, LR_state *s2) {
    float stage1, stage2;

    stage1 = coef->b[0]*(double)in + coef->b[1]*s1->x[0] + coef->b[2]*s1->x[1]
                - coef->a[1]*s1->y[0] - coef->a[2]*s1->y[1];

    s1->x[1] = s1->x[0];
    s1->x[0] = (double)in;
    s1->y[1] = s1->y[0];
    s1->y[0] = stage1;

    stage2 = coef->b[0]*(double)stage1 + coef->b[1]*s1->x[0] + coef->b[2]*s1->x[1]
                - coef->a[1]*s1->y[0] - coef->a[2]*s1->y[1];

    s1->x[1] = s1->x[0];
    s1->x[0] = (double)stage1;
    s1->y[1] = s1->y[0];
    s1->y[0] = stage2;

    return stage2;
}

void HP_coef(int sr, double freq, LR_coefficients *coef) {
        double q;
        double w0;
        double alpha;

        w0 = 2.0*M_PI*freq/sr;
        double sn = sin(w0);
        double cs = cos(w0);
        q = (double)1.0/sqrt(2.0);
        alpha = sn/(2.0*q);

        coef->b[0] = (1+cs) / 2.0;
        coef->b[1] = -(1+cs);
        coef->b[2] = (1+cs) / 2.0;
        coef->a[0] = 1+alpha;
        coef->a[1] = -2.0*cs;
        coef->a[2] = 1-alpha;

        coef->b[0] /= coef->a[0];
        coef->b[1] /= coef->a[0];
        coef->b[2] /= coef->a[0];
        coef->a[1] /= coef->a[0];
        coef->a[2] /= coef->a[0];
}

float HP_LRFilter(float in, LR_coefficients *coef, LR_state *s1, LR_state *s2) {
    float stage1, stage2;

    stage1 = coef->b[0]*(double)in + coef->b[1]*s1->x[0] + coef->b[2]*s1->x[1]
                - coef->a[1]*s1->y[0] - coef->a[2]*s1->y[1];

    s1->x[1] = s1->x[0];
    s1->x[0] = (double)in;
    s1->y[1] = s1->y[0];
    s1->y[0] = stage1;

    stage2 = coef->b[0]*(double)stage1 + coef->b[1]*s1->x[0] + coef->b[2]*s1->x[1]
                - coef->a[1]*s1->y[0] - coef->a[2]*s1->y[1];

    s1->x[1] = s1->x[0];
    s1->x[0] = (double)stage1;
    s1->y[1] = s1->y[0];
    s1->y[0] = stage2;

    return stage2;
}

