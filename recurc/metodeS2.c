/* 05-02-2024 */
/* alex */
/* metodeS2.c */
#include <stdio.h>
#include <stdlib.h>
#include "metodeS2.h"

void metode_setS2(int tam, double * cofs, double * res, int order) {
  int i;
  int vp[13] = {0, -1, 1, 2, 3, 5, 7, 11, 16, 24, 35, 53, 78};
  double auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[0];
  auxp = cofs[0]*cofs[0];
  for (i = 2; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS2(res, cofs[i], order);  
}
