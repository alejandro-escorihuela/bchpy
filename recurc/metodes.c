/* 29-01-2024 */
/* alex */
/* metodes.c */
#include <stdio.h>
#include <stdlib.h>
#include "metodes.h"

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
  res[0] = cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recA(res, cofs[i], order, rkn);
    else
      recB(res, cofs[i], order, rkn);
}

void metode_setBAB(int tam, double * cofs, double * res, int order, int rkn) {
  int i;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[1] = cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recB(res, cofs[i], order, rkn);
    else
      recA(res, cofs[i], order, rkn);
}

void metode_setXXsim(int tam, double * cofs, double * res, int order) {
  int i, s = tam/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recXasim(res, cofs[i], order);
    else
      recXbsim(res, cofs[i], order);
}

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

void metode_setS2sim(int tam, double * cofs, double * res, int order) {
}


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

void metode_setS4sim(int tam, double * cofs, double * res, int order) {
}


void metode_setS4(int tam, double * cofs, double * res, int order) {
  int i;
  int vp[15] = {0, -1, -1, -1, 1, 2, 3, 5, 7, 10, 13, 18, 24, 33, 44};
  double auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[0];
  res[vp[4]] = pow(cofs[0], 5);
  auxp = cofs[0]*cofs[0];
  for (i = 6; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS4(res, cofs[i], order);  
}

void metode_setABAsim_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
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

void metode_setBABsim_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
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

void metode_setABA_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[0] = cofs[0];
  for (i = 1; i <= tam; i++)
    if (i % 2 == 0)
      recA_c(res, cofs[i], order, rkn);
    else
      recB_c(res, cofs[i], order, rkn);
}

void metode_setBAB_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[1] = cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recB_c(res, cofs[i], order, rkn);
    else
      recA_c(res, cofs[i], order, rkn);
}

void metode_setXXsim_c(int tam, double complex * cofs, double complex * res, int order) {
  int i, s = tam/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recXasim_c(res, cofs[i], order);
    else
      recXbsim_c(res, cofs[i], order);
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

void metode_setS2sim_c(int tam, double complex * cofs, double complex * res, int order) {

}

void metode_setS2_c(int tam, double complex * cofs, double complex * res, int order) {
  int i;
  int vp[13] = {0, -1, 1, 2, 3, 5, 7, 11, 16, 24, 35, 53, 78};
  double complex auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[0];
  auxp = cofs[0]*cofs[0];
  for (i = 2; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS2_c(res, cofs[i], order);  
}

void metode_setS4sim_c(int tam, double complex * cofs, double complex * res, int order) {
}


void metode_setS4_c(int tam, double complex * cofs, double complex * res, int order) {
  int i;
  int vp[15] = {0, -1, -1, -1, 1, 2, 3, 5, 7, 10, 13, 18, 24, 33, 44};
  double complex auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[0];
  res[vp[4]] = cpow(cofs[0], 5);
  auxp = cofs[0]*cofs[0];
  for (i = 6; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS4_c(res, cofs[i], order);  
}
