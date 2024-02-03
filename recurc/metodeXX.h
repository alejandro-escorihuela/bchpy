/* 03-02-2024 */
/* alex */
/* metodeXX.h */
#ifndef _METODEXX_H
#define _METODEXX_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "recurXX.h"
#include "recurXXc.h"

void metode_setXX(int tam, double * cofs, double * res, int order);
void metode_setXX_c(int tam, double complex * cofs, double complex * res, int order);

#endif
