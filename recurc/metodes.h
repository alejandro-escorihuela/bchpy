/* 29-01-2024 */
/* alex */
/* metodes.h */
#ifndef _METODES_H
#define _METODES_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "recurAB.h"
#include "recurABc.h"
#include "recurXX.h"
#include "recurXXc.h"
#include "recurS2.h"
#include "recurS2c.h"
#include "recurS4.h"
#include "recurS4c.h"

void metode_setABAsim(int tam, double * cofs, double * res, int order, int rkn);
void metode_setBABsim(int tam, double * cofs, double * res, int order, int rkn);
void metode_setABA(int tam, double * cofs, double * res, int order, int rkn);
void metode_setBAB(int tam, double * cofs, double * res, int order, int rkn);
void metode_setXXsim(int tam, double * cofs, double * res, int order);
void metode_setXX(int tam, double * cofs, double * res, int order);
void metode_setS2sim(int tam, double * cofs, double * res, int order);
void metode_setS2(int tam, double * cofs, double * res, int order);
void metode_setS4sim(int tam, double * cofs, double * res, int order);
void metode_setS4(int tam, double * cofs, double * res, int order);


void metode_setABAsim_c(int tam, double complex * cofs, double complex * res, int order, int rkn);
void metode_setBABsim_c(int tam, double complex * cofs, double complex * res, int order, int rkn);
void metode_setABA_c(int tam, double complex * cofs, double complex * res, int order, int rkn);
void metode_setBAB_c(int tam, double complex * cofs, double complex * res, int order, int rkn);
void metode_setXXsim_c(int tam, double complex * cofs, double complex * res, int order);
void metode_setXX_c(int tam, double complex * cofs, double complex * res, int order);
void metode_setS2sim_c(int tam, double complex * cofs, double complex * res, int order);
void metode_setS2_c(int tam, double complex * cofs, double complex * res, int order);
void metode_setS4sim_c(int tam, double complex * cofs, double complex * res, int order);
void metode_setS4_c(int tam, double complex * cofs, double complex * res, int order);

#endif
