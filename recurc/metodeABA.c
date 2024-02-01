/* 29-01-2024 */
/* alex */
/* metodeABA.c */
#include <stdio.h>
#include <stdlib.h>
#include "metodeABA.h"

void metode_setABAsim(int tam, double * cofs, double * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  if (s % 2 == 0)
    res[0] = cofs[s];
  else
    res[1] = cofs[s];
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recAsim(res, cofs[i], order, rkn);
    else
      recBsim(res, cofs[i], order, rkn);
}

void metode_setBABsim(int tam, double * cofs, double * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  if (s % 2 == 0)
    res[1] = cofs[s];
  else
    res[0] = cofs[s];
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recBsim(res, cofs[i], order, rkn);
    else
      recAsim(res, cofs[i], order, rkn);
}

void metode_setABA(int tam, double * cofs, double * res, int order, int rkn) {
  int i;
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  if ((tam - 1) % 2 == 0)
    res[0] = cofs[tam - 1];
  else
    res[1] = cofs[tam - 1];
  for (i = tam - 2; i >= 0; i--)
    if (i % 2 == 0)
      recA(res, cofs[i], order, rkn);
    else
      recB(res, cofs[i], order, rkn);
}

/* void metode_setABA(int tam, double * cofs, double * res, int order, int rkn) { */
/*   int i; */
/*   double complex resc[TAMBCH], cofsc[TAMBCH]; */
  
/*   for (i = 0; i < TAMBCH; i++) */
/*     resc[i] = 0.0 + 0.0*I; */
/*   for (i = 0; i < tam; i++) */
/*     cofsc[i] = cofs[i] + 0.0*I; */
/*   if ((tam - 1) % 2 == 0) */
/*     resc[0] = cofsc[tam - 1]; */
/*   else */
/*     resc[1] = cofsc[tam - 1];   */
/*   for (i = tam - 2; i >= 0; i--) */
/*     if (i % 2 == 0) */
/*       recA(resc, cofsc[i], order, rkn); */
/*     else */
/*       recB(resc, cofsc[i], order, rkn); */
/*   for (i = 0; i < TAMBCH; i++) */
/*     res[i] = creal(resc[i]); */
/* } */

void metode_setBAB(int tam, double * cofs, double * res, int order, int rkn) {
  int i;
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  if ((tam - 1) % 2 == 0)
    res[1] = cofs[tam - 1];
  else
    res[0] = cofs[tam - 1];
  for (i = tam - 2; i >= 0; i--)
    if (i % 2 == 0)
      recB(res, cofs[i], order, rkn);
    else
      recA(res, cofs[i], order, rkn);
}
