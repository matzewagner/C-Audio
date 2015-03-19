//Multiband Frequency Shifter
//written by Matthias Wagner March 2015, for MAT 240B

#include <stdio.h>
#include <math.h>

#include "Hilbert.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif


void Hil_coef(int sr, H_coefficients *hc) {
    double poles[12] = {0.3609, 2.7412, 11.1573, 44.7581, 179.6242, 798.4578,
                        1.2524, 5.5671, 22.3423, 89.6271, 364.7914, 2770.1114};
    double poleFreq, rc, alpha, beta;
    int i;
    for (i=0; i<12; i++) {
        poleFreq = poles[i] * 15.0;
        rc = (double)(1.0 / (2.0*M_PI*poleFreq));
        alpha = (double)(1.0/ rc);
        alpha = alpha * 0.5 * (double)(1.0/sr);
        beta = (double)((1.0 - alpha) / (1.0 + alpha));
        hc->coef[i] = -(double)(beta);
    }
}

float Hil_real(float in, H_coefficients *hc) {
        double xn1, yn1;
        xn1 = in;
        int i;
        // compute the real (sine) output with 6 allpass filters
        for (i=0; i<6; i++) {
            yn1 = hc->coef[i] * (xn1 - hc->ynm1[i]) + hc->xnm1[i];
            hc->xnm1[i] = xn1;
            hc->ynm1[i] = yn1;
            xn1 = yn1;
        }
        float result = yn1;
        return result * 0.99;
}

float Hil_imaginary(float in, H_coefficients *hc) {
        double xn2, yn2;
        xn2 = in;
        int i;
        // compute the imaginary (cosine) output with 6 allpass filters
        for (i=6; i<12; i++) {
            yn2 = hc->coef[i] * (xn2 - hc->ynm1[i]) + hc->xnm1[i];
            hc->xnm1[i] = xn2;
            hc->ynm1[i] = yn2;
            xn2 = yn2;
        }
        float result = yn2;
        return result * 0.99;
}
