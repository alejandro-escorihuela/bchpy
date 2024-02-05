/* 03-02-2024 */
/* alex */
/* metodeXX.c */
#include <stdio.h>
#include <stdlib.h>
#include "metodeXX.h"

void metode_setXX(int tam, double * cofs, double * res, int order) {
  int i;
  int vp[9] = {0, 1, 2, 4, 7, 13, 22, 40, 70};

  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[0];
  for (i = 1; i < order; i++)
    res[vp[i]] = res[vp[i - 1]]*cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recXa(res, cofs[i], order);
    else
      recXb(res, cofs[i], order); 
}

void metode_setXX_c(int tam, double complex * cofs, double complex * res, int order) {
  int i;
  int vp[9] = {0, 1, 2, 4, 7, 13, 22, 40, 70};
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[0];
  for (i = 1; i < order; i++)
    res[vp[i]] = res[vp[i - 1]]*cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recXa_c(res, cofs[i], order);
    else
      recXb_c(res, cofs[i], order); 
}
