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
      recAsim_c(res, cofs[i], order, rkn);
    else
      recBsim_c(res, cofs[i], order, rkn);
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
      recBsim_c(res, cofs[i], order, rkn);
    else
      recAsim_c(res, cofs[i], order, rkn);
}
