/* 29-01-2024 */
/* alex */
/* metodeABA.h */
#ifndef _METODEABA_H
#define _METODEABA_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "recurABc.h"

void metode_setABAsim(int tam, double * cofs, double * res, int order, int rkn);
void metode_setBABsim(int tam, double * cofs, double * res, int order, int rkn);
void metode_setABA(int s, double * cofs, double * res, int order, int rkn);
void metode_setBAB(int s, double * cofs, double * res, int order, int rkn);
void metode_setABAcomp(int s, double complex * cofs, double complex * res, int order, int rkn);
void metode_setBABcomp(int s, double complex * cofs, double complex * res, int order, int rkn);
#endif
