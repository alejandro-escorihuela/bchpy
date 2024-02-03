/* 03-02-2024 */
/* alex */
/* copyesq.c */
#include <stdio.h>
#include <stdlib.h>
#include "esq.h"

void copyesq(double * dest, double * orig) {
  int i = 0;
  for (i = 0; i < TAMBCH; i++)
    dest[i] = orig[i];
}

void copyesq_c(double complex * dest, double complex * orig) {
  int i = 0;
  for (i = 0; i < TAMBCH; i++)
    dest[i] = orig[i];
}
