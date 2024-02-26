/* 28-02-2022 */
/* alex */
/* recurABc.h */
#ifndef _RECURABC_H
#define _RECURABC_H
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include "esq.h"

void recAsim_c(double complex * res, double complex x, int order, int rkn);
void recBsim_c(double complex * res, double complex x, int order, int rkn);
void recA_c(double complex * res, double complex x, int order, int rkn);
void recB_c(double complex * res, double complex x, int order, int rkn);

#endif
