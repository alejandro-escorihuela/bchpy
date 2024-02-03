/* 03-02-2024 */
/* alex */
/* esq.h */
#ifndef _ESQ_H
#define _ESQ_H
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#ifndef TAMBCH
#define TAMBCH 128
#endif

void copyesq(double * dest, double * orig);
void copyesq_c(double complex * dest, double complex * orig);

#endif
